
#include <algorithm>
#include <cerrno>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_set>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>

#include "Version.hpp"

namespace bio = boost::iostreams;


class FileException: public std::runtime_error {
public:
    FileException() : std::runtime_error("") { }
    explicit FileException(std::string msg) : std::runtime_error(msg) { }
};

class FASTQRecord {
public:
    std::string name;
    std::string comment;
    std::string sequence;
    std::string quality;

    void append_barcode(const std::string& barcode) {
        if (name.find(' ')!=std::string::npos) {
                name.insert(name.find(' '), "_" + barcode);
        } else {
                name += "_" + barcode;
        }
    }

    void strip_description() {
        // the description is the part of the name after the first space
        if (name.find(' ')!=std::string::npos) {
            name = name.substr(0, name.find(' '));
        }
    }

    std::string get_description() {
        // the description is the part of the name after the first space
        if (name.find(' ')!=std::string::npos) {
            return name.substr(name.find(' '));
        } else {
            return "";
        }
    }


};

std::ostream& operator<<(std::ostream& os, const FASTQRecord& record) {
    os << record.name << '\n'
       << record.sequence << '\n'
       << record.comment << '\n'
       << record.quality << '\n';
    return os;
}

std::istream& operator>>(std::istream& is, FASTQRecord& record) {
    getline(is, record.name);
    getline(is, record.sequence);
    getline(is, record.comment);
    getline(is, record.quality);
    return is;
}

///
/// Check for the GZIP header in a file
///
bool is_gzipped_file(const std::string& filename) {
    bool gzipped = false;
    FILE* f = fopen(filename.c_str(), "rb");

    if (f == NULL) {
        throw FileException("Could not open file \"" + filename + "\": " + strerror(errno));
    } else {
        if (fgetc(f) == 0x1f && fgetc(f) == 0x8b) {
            gzipped = true;
        }
    }
    fclose(f);
    return gzipped;
}

bool is_gzipped_filename(const std::string& filename) {
    std::string ext = ".gz";
    return std::equal(ext.rbegin(), ext.rend(), filename.rbegin());
}

std::regex fastq_filename_re("\\.(fq|fastq)?(\\.gz)?$");
std::string make_trimmed_filename(const std::string& filename) {
    return std::regex_replace(filename, fastq_filename_re, ".trimmed.fq.gz");
}

std::map<char, char> nucleotide_complements = {
    {'A', 'T'},
    {'T', 'A'},
    {'C', 'G'},
    {'G', 'C'},
    {'N', 'N'},
    {'a', 't'},
    {'t', 'a'},
    {'c', 'g'},
    {'g', 'c'},
    {'n', 'n'},
};

std::string complement(const std::string& seq) {
    std::string complemented;
    for (const auto it : seq) {
        complemented += nucleotide_complements.at(it);
    }
    return complemented;
}

std::string reverse_complement(const std::string& seq) {
    std::string rc = complement(seq);
    std::reverse(rc.begin(), rc.end());
    return rc;
}

bool shares_kmer(const std::string& s1, const std::string& s2, int k) {
    // return true if any kmer from s1 is in s2, else false

    // get all s1 kmers
    std::unordered_set<std::string> kmers;
    for (int i = 0; i <= s1.size() - k; i++) {
        kmers.insert(s1.substr(i, k));
    }


    // check if any kmer from s2 is in the set
    for (int i = 0; i <= s2.size() - k; i++) {
        if (kmers.count(s2.substr(i, k)) > 0) {
            return true;
        }
    }

    return false;
}

// based on https://en.wikipedia.org/wiki/Levenshtein_distance
long long int levenshtein_distance(const std::string& s1, const std::string& s2, long long int maximum = 0) {
    // shortcuts
    if (s1 == s2) return 0;

    long long int s1_length = s1.size();
    long long int s2_length = s2.size();

    if (s1_length == 0) {
        return s2_length;
    }

    if (s2_length == 0) {
        return s1_length;
    }

    if (maximum) {
        long long int minimum_distance = llabs(s1_length - s2_length);
        if (maximum < minimum_distance) {
            return minimum_distance;
        }
    }

    // ok, not getting out of the work that easily
    std::vector<long long int> v1(s2_length + 1);
    std::vector<long long int> v2(s2_length + 1);

    for (long long int i = 0; i < s2_length; i++) {
        v1[i] = i;
    }

    long long int v1_size = v1.size();
    for (long long int i = 0; i < s1_length; i++) {
        v2[0] = i + 1;

        for (long long int j = 0; j < s2_length; j++) {
            int cost = s1[i] == s2[j] ? 0 : 1;
            v2[j + 1] = std::min(v2[j] + 1, std::min(v1[j + 1] + 1, v1[j] + cost));
        }

        long long int minimum = maximum + 1;
        for (long long int j = 0; j < v1_size; j++) {
            if (v2[j] < minimum) {
                minimum = v2[j];
            }
            v1[j] = v2[j];
        }
        if (minimum > maximum) {
            // abort early
            return maximum + 1;
        }
    }

    return v2[s2_length];
}

/// Returns the position of the best alignment of the two sequences
long long int align(const std::string& seq1, const std::string& seq2, long long int max_edit_distance, bool verbose = false) {
    long long int seq1_length = seq1.size();
    long long int seq2_length = seq2.size();
    std::string window;
    long long int distance = -1;
    long long int best_position = -1;
    std::vector<std::pair<long long int, long long int>> alignments;

    int min_shared_kmer_size = seq1_length / (max_edit_distance + 1); // For two strings to have edit distance <= max_edit_distance, they must share a kmer of at least this size

    if (verbose) {
        std::cout << "align: \n" << seq1 << '\n' << seq2 << '\n';
    }

    if (!shares_kmer(seq1, seq2, min_shared_kmer_size)) {
        if (verbose) {
            std::cout << "align: no shared kmer of length >= " << min_shared_kmer_size << " between these sequences, so cannot have edit distance <= " << max_edit_distance << "\n";
        }
        return best_position;
    }

    for (long long int i = 0; i < seq2_length; i++) {
        window = seq2.substr(i, seq1_length);
        distance = levenshtein_distance(window, seq1, max_edit_distance);

        if (verbose) {
            std::cout << "at position " << i << " distance " << distance << " between\n" << seq1 << '\n' << window << '\n';
        }

        if (distance <= max_edit_distance) {
            alignments.push_back(std::pair<long long int, long long int>(distance, i));
        }
    }

    std::sort(alignments.begin(), alignments.end());
    if (verbose) {
        for (auto it : alignments ) {
            std::cout << it.first << ':' << it.second << '\n';
        }
    }

    if (!alignments.empty()) {
        best_position = alignments[0].second;
    }

    return best_position;
}

void trim_fastq_record(FASTQRecord& record, long long int start, long long int end) {
    record.sequence = record.sequence.substr(start, end);
    record.quality = record.quality.substr(start, end);
}

void trim_pair(FASTQRecord& rec1, FASTQRecord& rec2, long long int rc_length, long long int max_edit_distance, long long int fudge = 0, long long int trim_start = 0, bool verbose = false) {
    long long int alignment_start = -1;
    std::string rec2_rc = reverse_complement(rec2.sequence.substr(0, rc_length));
    if (verbose) {
        std::cout << rec1.sequence << '\n' << rec2.sequence << '\n' << rec2_rc << '\n';
    }

    size_t position_of_rec2_rc_in_rec1 = rec1.sequence.rfind(rec2_rc);
    long long int trim_end = rc_length - fudge;

    if (position_of_rec2_rc_in_rec1 != std::string::npos) {
        alignment_start = position_of_rec2_rc_in_rec1;

        if (verbose) {
            std::cout << "trim_pair: rec2 rc " << rec2_rc << " found at " << position_of_rec2_rc_in_rec1 << " in rec1 " << rec1.sequence  << '\n';
        }
    } else if (max_edit_distance > 0) {
        alignment_start = align(rec2_rc, rec1.sequence, max_edit_distance, verbose);

        if (verbose) {
            std::cout << "trim_pair: rec2 rc " << rec2_rc << " not found in " << rec1.sequence << "; alignment result: " << alignment_start << '\n';
        }
    }

    if (alignment_start > 0) {
        if (verbose) {
            std::cout << "trim_pair: best alignment at " << alignment_start << '\n';
        }
        if (alignment_start > 0) {
            if (verbose) {
                std::cout << "trim_pair: trimming pair " << rec1.name << '\n';
            }
            trim_fastq_record(rec1, trim_start, alignment_start + trim_end);
            trim_fastq_record(rec2, trim_start, alignment_start + trim_end);
        }
    } else if (verbose) {
        std::cout << "not touching pair " << rec1.name << '\n';
    }
}

std::string version_string() {
    std::stringstream ss;
    ss << PROGRAM_NAME << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH;
    return ss.str();
}

void print_usage() {
    std::cout << version_string() << ": Trim adapters from paired-end HTS reads.\n\n"

              << "Usage:\n\n" << PROGRAM_NAME << " [options] input_file1 input_file2 [output_file1 output_file2]\n\n"
              << "where:\n"
              << "    input_file1 contains the first reads of all pairs\n"
              << "    input_file2 contains the second reads of all pairs\n"
              << "    output_file1 will contain the trimmed first reads of all pairs\n"
              << "    output_file2 will contain the trimmed second reads of all pairs\n\n"

              << "You must supply both of the output file names, or neither. If not\n"
              << "specified, they will be inferred from the inputs.\n\n"

              << "Options may include:\n\n"

              << "-h|--help: show this usage message.\n\n"
              << "-v|--verbose: show more details and progress updates.\n\n"
              << "--version: print the version of the program." << std::endl << std::endl

              << "-m|--max-edit-distance\n"
              << "    The maximum edit distance permitted when aligning the paired reads. The default is 1.\n"

              << "-f|--fudge\n"
              << "    An arbitrary number of extra bases to trim, to make sure the reads and mates\n"
              << "    don't overlap exactly. Some aligners (e.g. bowtie) apparently don't like perfectly\n"
              << "    overlapping reads.\n"

              << "-r|--rc-length\n"
              << "    Use the reverse complement of this number of bases from the beginning of the\n"
              << "    reverse read to align the reads. The default is 20.\n"

              << "-t|--trim-from-start\n"
              << "    Trim this number of bases from the start of each sequence. The default is zero.\n"

              << "-a|--append-barcode\n"
              << "    Append barcodes from this fastq file to the end of read names (useful for e.g. single-cell data).\n"

              << "-d|--copy-description\n"
              << "    Copy read descriptions (the part of the read name after the first space) from this fastq file\n"
              << "    and add them as descriptions to the paired reads (useful when e.g. using bwa mem with the -C flag\n"
              << "    for alignment).\n"

              << "-s|--strip-description\n"
              << "    Strip read descriptions (the part of the read name after the first space). For example, read name\n"
              << "    '@LH00346:6:22CGGKLT3:1:1101:2387:1032:rCAGACGCGGCTTAAAGAAATGAGG 2:N:0:AGTTAGTT' would be changed to\n"
              << "    '@LH00346:6:22CGGKLT3:1:1101:2387:1032:rCAGACGCGGCTTAAAGAAATGAGG'. Useful in \n"
              << "     conjunction with --copy-description.\n"

              << "-x|--selected-barcodes\n"
              << "    Select and output reads corresponding to these barcodes from this unzipped text file\n"

              << std::endl;
}

bool getFileContent(std::string fileName, std::map<std::string,bool> & vecOfStrs)
{

    // Open the File
    std::ifstream in(fileName.c_str());

    // Check if object is valid
    if(!in)
        {
            std::cerr << "Cannot open the File : "<<fileName<<std::endl;
            return false;
        }

    std::string str;
    // Read the next line from File untill it reaches the end.
    while (std::getline(in, str))
        {
            // Line contains string of length > 0 then save it in vector
            if(str.size() > 0)
                vecOfStrs[str] = true;
        }
    //Close The File
    in.close();
    return true;
}
 
int main(int argc, char **argv)
{
    int c, option_index = 0;
    bool verbose = 0;
    bool strip_description = false;

    long long int max_edit_distance = 1;
    long long int fudge = 0;
    long long int trim_start = 0;
    long long int rc_length = 20;
    std::string input_filename1;
    std::string input_filename2;
    std::string barcode_filename;
    std::string selected_barcode_filename;
    std::string description_filename;
    std::string output_filename1;
    std::string output_filename2;

    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", no_argument, NULL, 'v'},
        {"version", no_argument, NULL, 0},
        {"max-edit-distance", required_argument, NULL, 'd'},
        {"fudge", required_argument, NULL, 'f'},
        {"trim-from-start", required_argument, NULL, 't'},
        {"rc-length", required_argument, NULL, 'r'},
        {"append-barcode", required_argument, NULL, 'a'},
        {"selected-barcodes", required_argument, NULL, 'x'},
        {"copy-description", required_argument, NULL, 'd'},
        {"strip-description", no_argument, NULL, 's'},
        {0, 0, 0, 0}
    };

    // parse the command line arguments
    while ((c = getopt_long(argc, argv, "d:f:t:r:a:x:d:s:vh?", long_options, &option_index)) != -1) {
        switch (c) {
        case 'h':
        case '?':
            print_usage();
            exit(1);
        case 'v':
            verbose = true;
            break;
        case 'm':
            max_edit_distance = std::stoll(optarg);
            break;
        case 'f':
            fudge = std::stoll(optarg);
            break;
        case 't':
            trim_start = std::stoll(optarg);
            break;
        case 'r':
            rc_length = std::stoll(optarg);
            break;
        case 'a':
            barcode_filename = optarg;
            break;
        case 'x':
            selected_barcode_filename = optarg;
            break;
        case 's':
            strip_description = true;
            break;
        case 'd':
            description_filename = optarg;
            break;
        case 0:
            if (long_options[option_index].flag != 0){
                break;
            }

            if (long_options[option_index].name == std::string("version")) {
                std::cout << version_string() << std::endl;
                exit(1);
            }

        default:
            print_usage();
            exit(1);
        }
    }

    if (argc < (optind + 2)) {
        print_usage();
        exit(1);
    }

    input_filename1 = argv[optind];
    input_filename2 = argv[optind + 1];

    if (argc > (optind + 2)) {
        if (argc == optind + 4) {
            output_filename1 = argv[optind + 2];
            output_filename2 = argv[optind + 3];
        } else {
            print_usage();
            exit(1);
        }
    }

    bio::file_source input_file1(input_filename1);
    bio::stream<bio::file_source> input_stream1(input_file1);

    bio::file_source input_file2(input_filename2);
    bio::stream<bio::file_source> input_stream2(input_file2);

    bio::file_source barcode_file(barcode_filename);
    bio::stream<bio::file_source> barcode_stream(barcode_file);

    bio::file_source description_file(description_filename);
    bio::stream<bio::file_source> description_stream(description_file);

    bio::filtering_stream<bio::input> in1;
    if (is_gzipped_file(input_filename1)) {
        in1.push(bio::gzip_decompressor());
    }
    in1.push(input_stream1);

    bio::filtering_stream<bio::input> in2;
    if (is_gzipped_file(input_filename2)) {
        in2.push(bio::gzip_decompressor());
    }
    in2.push(input_stream2);

    bio::filtering_stream<bio::input> inBarcode;
    if (!barcode_filename.empty() && is_gzipped_file(barcode_filename)) {
        inBarcode.push(bio::gzip_decompressor());
    }
    inBarcode.push(barcode_stream);

    bio::filtering_stream<bio::input> inDescription;
    if (!description_filename.empty() && is_gzipped_file(description_filename)) {
        inDescription.push(bio::gzip_decompressor());
    }
    inDescription.push(description_stream);

    if (output_filename1.empty()) {
        output_filename1 = make_trimmed_filename(input_filename1);
    }

    if (output_filename2.empty()) {
        output_filename2 = make_trimmed_filename(input_filename2);
    }

    if (verbose) {
        std::cout << "Writing to " << output_filename1 << " and " << output_filename2 << "...\n";
    }

    std::ofstream output_file1(output_filename1, std::ios_base::out | std::ios_base::binary);
    std::ofstream output_file2(output_filename2, std::ios_base::out | std::ios_base::binary);

    bio::filtering_stream<bio::output> out1;
    if (is_gzipped_filename(output_filename1)) {
        out1.push(bio::gzip_compressor(bio::gzip_params(bio::gzip::best_speed)));
    }
    out1.push(output_file1);

    bio::filtering_stream<bio::output> out2;
    if (is_gzipped_filename(output_filename2)) {
        out2.push(bio::gzip_compressor(bio::gzip_params(bio::gzip::best_speed)));
    }
    out2.push(output_file2);

    FASTQRecord record1;
    FASTQRecord record2;
    FASTQRecord recordBarcode;
    FASTQRecord recordDescription;

    std::string barcode;

    std::map<std::string,bool> selected_barcodes;
    bool result = (!barcode_filename.empty() & !selected_barcode_filename.empty()) ? getFileContent(selected_barcode_filename, selected_barcodes) : true;


    while ((in1 >> record1) && (in2 >> record2)){

        if (!barcode_filename.empty()) {
            inBarcode >> recordBarcode;
            barcode = recordBarcode.sequence;
        }

        if (!description_filename.empty()) {
            inDescription >> recordDescription;
        }

        if (!barcode_filename.empty() && !selected_barcode_filename.empty()) {
            if (selected_barcodes.count(barcode) == 0) {
                continue;
            }
        }

        if (strip_description) {
            record1.strip_description();
            record2.strip_description();
        }

        if (!description_filename.empty()) {
            std::string description = recordDescription.get_description();
            record1.name += description;
            record2.name += description;
        }

        if (!barcode_filename.empty()) {
            record1.append_barcode(barcode);
            record2.append_barcode(barcode);
        }

        trim_pair(record1, record2, rc_length, max_edit_distance, fudge, trim_start, verbose);

        out1 << record1;
        out2 << record2;
    }

    bio::close(out1);
    bio::close(out2);

    bio::close(in1);
    bio::close(in2);
    bio::close(inBarcode);
    bio::close(inDescription);

    return 0;
}

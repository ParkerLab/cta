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

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>


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
bool is_gzipped_file(std::string filename) {
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

std::regex fastq_filename_re("\\.(fq|fastq)?(\\.gz)?$");
std::string make_trimmed_filename(std::string filename) {
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

// based on https://en.wikipedia.org/wiki/Levenshtein_distance
long long levenshtein_distance(const std::string& s1, const std::string& s2, long long maximum = 0) {
    // shortcuts
    if (s1 == s2) return 0;

    long long s1_length = s1.size();
    long long s2_length = s2.size();

    if (s1_length == 0) {
        return s2_length;
    }

    if (s2_length == 0) {
        return s1_length;
    }

    if (maximum) {
        long long minimum_distance = llabs(s1_length - s2_length);
        if (maximum < minimum_distance) {
            return minimum_distance;
        }
    }

    // ok, not getting out of the work that easily
    std::vector<long long> v1(s2_length + 1);
    std::vector<long long> v2(s2_length + 1);

    for (long long i = 0; i < s2_length; i++) {
        v1[i] = i;
    }

    long long v1_size = v1.size();
    for (long long i = 0; i < s1_length; i++) {
        v2[0] = i + 1;

        for (long long j = 0; j < s2_length; j++) {
            int cost = s1[i] == s2[j] ? 0 : 1;
            v2[j + 1] = std::min(v2[j] + 1, std::min(v1[j + 1] + 1, v1[j] + cost));
        }

        for (long long j = 0; j < v1_size; j++) {
            v1[j] = v2[j];
        }
    }

    return v2[s2_length];
}

/// Returns the position of the best alignment of the two sequences
long long align(const std::string& seq1, const std::string& seq2, long max_edit_distance, bool verbose = false) {
    long long seq1_length = seq1.size();
    long long seq2_length = seq2.size();
    std::string window;
    long long distance = -1;
    std::vector<std::pair<long, long long>> alignments;

    if (verbose) {
        std::cout << "align: \n" << seq1 << '\n' << seq2 << '\n';
    }

    for (long long i = 0; i < seq2_length; i++) {
        window = seq2.substr(i, seq1_length);
        distance = levenshtein_distance(window, seq1, max_edit_distance);

        if (verbose) {
            std::cout << "at position " << i << " distance " << distance << " between\n" << seq1 << '\n' << window << '\n';
        }

        if (distance <= max_edit_distance) {
            alignments.push_back(std::pair<long, long long>(distance, i));
        }
    }

    std::sort(alignments.begin(), alignments.end());
    if (verbose) {
        for (auto it : alignments ) {
            std::cout << it.first << ':' << it.second << '\n';
        }
    }

    long long best_position = -1;
    if (!alignments.empty()) {
        best_position = alignments[0].second;
    }

    return best_position;
}

void trim_fastq_record(FASTQRecord& record, long long start, long long end) {
    record.sequence = record.sequence.substr(start, end);
    record.quality = record.quality.substr(start, end);
}

void trim_pair(FASTQRecord& rec1, FASTQRecord& rec2, long rc_length, long max_edit_distance, long fudge=0, long trim_start=0, bool verbose = false) {
    long long alignment_start = -1;
    std::string rec2_rc = reverse_complement(rec2.sequence.substr(0, rc_length));
    if (verbose) {
        std::cout << rec1.sequence << '\n' << rec2.sequence << '\n' << rec2_rc << '\n';
    }

    size_t position_of_rec2_rc_in_rec1 = rec1.sequence.rfind(rec2_rc);
    long long trim_end = rc_length - fudge;

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

void print_usage() {
    std::cout << "trim_adapters: Trim adapters from paired-end HTS reads.\n\n"

              << "Usage:\n\ntrim_adapters [options] file1 file2\n\n"
              << "where:\n"
              << "    file1 contains the first reads of all the pairs\n"
              << "    file2 contains the second reads of all the pairs\n\n"

              << "Options may include:\n\n"

              << "-h|--help: show this usage message.\n\n"
              << "-v|--verbose: show more details and progress updates.\n\n"

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

              << std::endl;
}

int main(int argc, char **argv)
{
    namespace bio = boost::iostreams;

    int c, option_index = 0;
    bool verbose = 0;

    long long max_edit_distance = 1;
    long long fudge = 0;
    long long trim_start = 0;
    long long rc_length = 20;
    std::string input_filename1;
    std::string input_filename2;

    static struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", no_argument, NULL, 'v'},
        {"max-edit-distance", required_argument, NULL, 'd'},
        {"fudge", required_argument, NULL, 'f'},
        {"trim-from-start", required_argument, NULL, 't'},
        {"rc-length", required_argument, NULL, 'r'},
        {0, 0, 0, 0}
    };

    // parse the command line arguments
    while ((c = getopt_long(argc, argv, "d:f:t:r:vh?", long_options, &option_index)) != -1) {
        switch (c) {
        case 'h':
        case '?':
            print_usage();
            exit(1);
        case 'v':
            verbose = true;
            break;
        case 'd':
            max_edit_distance = std::stoll(optarg);
            break;
        case 'f':
            fudge = std::stoll(optarg);
            break;
        case 's':
            trim_start = std::stoll(optarg);
            break;
        case 'r':
            rc_length = std::stoll(optarg);
            break;
        default:
            print_usage();
            exit(1);
        }
    }

    if ((optind + 2) > argc) {
        print_usage();
        exit(1);
    }

    input_filename1 = argv[optind];
    input_filename2 = argv[optind + 1];

    bio::file_source input_file1(input_filename1);
    bio::stream<bio::file_source> input_stream1(input_file1);

    bio::file_source input_file2(input_filename2);
    bio::stream<bio::file_source> input_stream2(input_file2);

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

    std::string output_filename1 = make_trimmed_filename(input_filename1);
    std::string output_filename2 = make_trimmed_filename(input_filename2);

    if (verbose) {
        std::cout << "Writing to " << output_filename1 << " and " << output_filename2 << "...\n";
    }

    std::ofstream output_file1(output_filename1, std::ios_base::out | std::ios_base::binary);
    std::ofstream output_file2(output_filename2, std::ios_base::out | std::ios_base::binary);

    bio::filtering_stream<bio::output> out1;
    out1.push(bio::gzip_compressor(bio::gzip_params(bio::gzip::best_speed)));
    out1.push(output_file1);

    bio::filtering_stream<bio::output> out2;
    out2.push(bio::gzip_compressor(bio::gzip_params(bio::gzip::best_speed)));
    out2.push(output_file2);

    FASTQRecord record1;
    FASTQRecord record2;

    while ((in1 >> record1) && (in2 >> record2)) {
        trim_pair(record1, record2, rc_length, max_edit_distance, fudge, trim_start, verbose);
        out1 << record1;
        out2 << record2;
    }
    bio::close(out1);
    bio::close(out2);
}

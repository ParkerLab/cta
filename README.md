# cta

`cta` is used to trim adapters from paired-end high-throughput sequencing reads.

## Building

Compilation required `boost` libraries. Use `make` to build binaries in `build/`
directory, and `make install` to install in `$PREFIX/bin`. Adjust `PREFIX`
accordingly to install in a custom location.

```
make install  # installs cta in /usr/bin
```

## Usage

```
cta 0.1.3: Trim adapters from paired-end HTS reads.

Usage:

cta [options] input_file1 input_file2 [output_file1 output_file2]

where:
    input_file1 contains the first reads of all pairs
    input_file2 contains the second reads of all pairs
    output_file1 will contain the trimmed first reads of all pairs
    output_file2 will contain the trimmed second reads of all pairs

You must supply both of the output file names, or neither. If not
specified, they will be inferred from the inputs.

If using barcode options, ensure that your files contain the correct
information. For example, if using -a, cta will only process number of
FASTQ reads equal to the number of records in the barcode index file.

Options may include:

-h|--help: show this usage message.

-v|--verbose: show more details and progress updates.

--version: print the version of the program.

-m|--max-edit-distance
    The maximum edit distance permitted when aligning the paired reads. The default is 1.
-f|--fudge
    An arbitrary number of extra bases to trim, to make sure the reads and mates
    don't overlap exactly. Some aligners (e.g. bowtie) apparently don't like perfectly
    overlapping reads.
-r|--rc-length
    Use the reverse complement of this number of bases from the beginning of the
    reverse read to align the reads. The default is 20.
-t|--trim-from-start
    Trim this number of bases from the start of each sequence. The default is zero.
-a|--append-barcode
    Append barcodes from the given FASTQ file to the end of read names (useful for single-cell data).
-x|--selected-barcodes
    Only output reads corresponding to barcodes in the given file.
```

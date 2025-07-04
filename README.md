[![install with conda](
https://anaconda.org/bioconda/fastp/badges/version.svg)](https://anaconda.org/bioconda/fastp)
[![install with conda](
https://anaconda.org/bioconda/fastp/badges/downloads.svg)](https://anaconda.org/bioconda/fastp)
[![DebianBadge](
https://badges.debian.net/badges/debian/unstable/fastp/version.svg)](https://packages.debian.org/unstable/fastp)

# fastp
A tool designed to provide ultrafast all-in-one preprocessing and quality control for FastQ data.     

This tool is designed for processing short reads (i.e. Illumina NovaSeq, MGI), if you are looking for tools to process long reads (i.e. Nanopore, PacBio, Cyclone), please use [fastplong](https://github.com/OpenGene/fastplong).  

fastp supports batch processing of multiple FASTQ files in a folder, see - [batch processing](#batch-processing)


- [features](#features)
- [simple usage](#simple-usage)
- [examples of report](#examples-of-report)
- [get fastp](#get-fastp)
  - [install with Bioconda](#install-with-bioconda)
  - [or download the latest prebuilt binary for Linux users](#or-download-the-latest-prebuilt-binary-for-linux-users)
  - [or compile from source](#or-compile-from-source)
    - [Step 1: install isa-l](#step-1-install-isa-l)
    - [step 2: install libdeflate](#step-2-install-libdeflate)
    - [Step 3: download and build fastp](#step-3-download-and-build-fastp)
- [input and output](#input-and-output)
  - [output to STDOUT](#output-to-stdout)
  - [input from STDIN](#input-from-stdin)
  - [store the unpaired reads for PE data](#store-the-unpaired-reads-for-pe-data)
  - [store the reads that fail the filters](#store-the-reads-that-fail-the-filters)
  - [process only part of the data](#process-only-part-of-the-data)
  - [do not overwrite exiting files](#do-not-overwrite-exiting-files)
  - [split the output to multiple files for parallel processing](#split-the-output-to-multiple-files-for-parallel-processing)
  - [merge PE reads](#merge-pe-reads)
- [filtering](#filtering)
  - [quality filter](#quality-filter)
  - [length filter](#length-filter)
  - [low complexity filter](#low-complexity-filter)
  - [Other filter](#other-filter)
- [adapters](#adapters)
- [per read cutting by quality score](#per-read-cutting-by-quality-score)
- [base correction for PE data](#base-correction-for-pe-data)
- [global trimming](#global-trimming)
- [polyG tail trimming](#polyg-tail-trimming)
- [polyX tail trimming](#polyx-tail-trimming)
- [unique molecular identifier (UMI) processing](#unique-molecular-identifier-umi-processing)
  - [UMI example](#umi-example)
- [output splitting](#output-splitting)
  - [splitting by limiting file number](#splitting-by-limiting-file-number)
  - [splitting by limiting the lines of each file](#splitting-by-limiting-the-lines-of-each-file)
- [overrepresented sequence analysis](#overrepresented-sequence-analysis)
- [merge paired-end reads](#merge-paired-end-reads)
- [duplication rate and deduplication](#duplication-rate-and-deduplication)
  - [duplication rate evaluation](#duplication-rate-evaluation)
  - [deduplication](#deduplication)
- [batch processing](#batch-processing)
- [all options](#all-options)
- [citations](#citations)

# features
0. comprehensive quality profiling for both before and after filtering data (quality curves, base contents, KMER, Q20/Q30, GC Ratio, duplication, adapter contents...)
1. filter out bad reads (too low quality, too short, or too many N...)
2. cut low quality bases for per read in its 5' and 3' by evaluating the mean quality from a sliding window (like Trimmomatic but faster).
3. trim all reads in front and tail
4. cut adapters. Adapter sequences can be automatically detected, which means you don't have to input the adapter sequences to trim them.
5. correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
6. trim polyG in 3' ends, which is commonly seen in NovaSeq/NextSeq data. Trim polyX in 3' ends to remove unwanted polyX tailing (i.e. polyA tailing for mRNA-Seq data)
7. preprocess unique molecular identifier (UMI) enabled data, shift UMI to sequence name.
8. report JSON format result for further interpreting.
9. visualize quality control and filtering results on a single HTML page (like FASTQC but faster and more informative).
10. split the output to multiple files (0001.R1.gz, 0002.R1.gz...) to support parallel processing. Two modes can be used, limiting the total split file number, or limitting the lines of each split file.
11. support long reads (data from PacBio / Nanopore devices).
12. support reading from STDIN and writing to STDOUT
13. support interleaved input
14. support ultra-fast FASTQ-level deduplication
15. ...

If you find a bug or have additional requirement for `fastp`, please file an issue:https://github.com/OpenGene/fastp/issues/new

# simple usage
* for single end data (not compressed)
```
fastp -i in.fq -o out.fq
```
* for paired end data (gzip compressed)
```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```
By default, the HTML report is saved to `fastp.html` (can be specified with `-h` option), and the JSON report is saved to `fastp.json` (can be specified with `-j` option).

# examples of report
`fastp` creates reports in both HTML and JSON format.
* HTML report: http://opengene.org/fastp/fastp.html
* JSON report: http://opengene.org/fastp/fastp.json

# get fastp
## install with Bioconda
[![install with conda](
https://anaconda.org/bioconda/fastp/badges/version.svg)](https://anaconda.org/bioconda/fastp)
```shell
# note: the fastp version in bioconda may be not the latest
conda install -c bioconda fastp
```
## or download the latest prebuilt binary for Linux users
This binary was compiled on CentOS, and tested on CentOS/Ubuntu
```shell
# download the latest build
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# or download specified version, i.e. fastp v0.23.4
wget http://opengene.org/fastp/fastp.0.23.4
mv fastp.0.23.4 fastp
chmod a+x ./fastp
```
## or compile from source
`fastp` depends on `libdeflate` and `libisal`, while `libisal` is not compatible with gcc 4.8. If you use gcc 4.8, your fastp will fail to run. Please upgrade your gcc before you build the libraries and fastp.

### Step 1: install isa-l
It's recommended that to install it using your package manager, for example `apt install isa-l` on ubuntu, or `brew install isa-l` on Mac. Otherwise you can compile it from source. Please be noted that `isa-l` is not compatible with gcc 4.8 or older versions. See https://github.com/intel/isa-l
`autoconf`, `automake`, `libtools`, `nasm (>=2.11.01)` and `yasm (>=1.2.0)` are required to build isa-l.
```shell
git clone https://github.com/intel/isa-l.git
cd isa-l
./autogen.sh
./configure --prefix=/usr --libdir=/usr/lib64
make -j
sudo make install
```

### step 2: install libdeflate
It's recommended that to install it using your package manager, for example `apt install libdeflate` on ubuntu, or `brew install libdeflate` on Mac. Otherwise you can compile it from source. See https://github.com/ebiggers/libdeflate
```shell
git clone https://github.com/ebiggers/libdeflate.git
cd libdeflate
cmake -B build
cmake --build build
cmake --install build
```

### Step 3: download and build fastp
```shell
# get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/fastp.git

# build
cd fastp
make -j

# Install
sudo make install
```

# input and output
`fastp` supports both single-end (SE) and paired-end (PE) input/output.
* for SE data, you only have to specify read1 input by `-i` or `--in1`, and specify read1 output by `-o` or `--out1`.
* for PE data, you should also specify read2 input by `-I` or `--in2`, and specify read2 output by `-O` or `--out2`.
* if you don't specify the output file names, no output files will be written, but the QC will still be done for both data before and after filtering.
* the output will be gzip-compressed if its file name ends with `.gz`
## output to STDOUT
`fastp` supports streaming the passing-filter reads to STDOUT, so that it can be passed to other compressors like `bzip2`, or be passed to aligners like `bwa` and `bowtie2`.
* specify `--stdout` to enable this mode to stream output to STDOUT
* for PE data, the output will be interleaved FASTQ, which means the output will contain records like `record1-R1 -> record1-R2 -> record2-R1 -> record2-R2 -> record3-R1 -> record3-R2 ... `
## input from STDIN
* specify `--stdin` if you want to read the STDIN for processing.
* if the STDIN is an interleaved paired-end stream, specify `--interleaved_in` to indicate that.
## store the unpaired reads for PE data
* you can specify `--unpaired1` to store the reads that read1 passes filters but its paired read2 doesn't, as well as `--unpaired2` for unpaired read2.
* `--unpaired1` and `--unpaired2` can be the same, so the unpaired read1/read2 will be written to the same single file.
## store the reads that fail the filters
* give `--failed_out` to specify the file name to store the failed reads.
* if one read failed and is written to `--failed_out`, its `failure reason` will be appended to its read name. For example, `failed_quality_filter`, `failed_too_short` etc.
* for PE data, if unpaired reads are not stored (by giving --unpaired1 or --unpaired2), the failed pair of reads will be put together. If one read passes the filters but its pair doesn't, the `failure reason` will be `paired_read_is_failing`.
## process only part of the data
If you don't want to process all the data, you can specify `--reads_to_process` to limit the reads to be processed. This is useful if you want to have a fast preview of the data quality, or you want to create a subset of the filtered data.
## do not overwrite exiting files
You can enable the option `--dont_overwrite` to protect the existing files not to be overwritten by `fastp`. In this case, `fastp` will report an error and quit if it finds any of the output files (read1, read2, json report, html report) already exists before.
## split the output to multiple files for parallel processing
See [output splitting](#output-splitting)
## merge PE reads
See [merge paired-end reads](#merge-paired-end-reads)

# filtering
Multiple filters have been implemented.
## quality filter
Quality filtering is enabled by default, but you can disable it by `-Q` or `disable_quality_filtering`. Currently it supports filtering by limiting the N base number (`-n, --n_base_limit`),  and the percentage of unqualified bases.  

To filter reads by its percentage of unqualified bases, two options should be provided:
* `-q, --qualified_quality_phred`       the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
* `-u, --unqualified_percent_limit`    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%

You can also filter reads by its average quality score
* `-e, --average_qual`   if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])

## length filter
Length filtering is enabled by default, but you can disable it by `-L` or `--disable_length_filtering`. The minimum length requirement is specified with `-l` or `--length_required`.

For some applications like small RNA sequencing, you may want to discard the long reads. You can specify `--length_limit` to discard the reads longer than `length_limit`. The default value 0 means no limitation.

## low complexity filter
Low complexity filter is disabled by default, and you can enable it by `-y` or `--low_complexity_filter`. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]). For example:
```
# a 51-bp sequence, with 3 bases that is different from its next base
seq = 'AAAATTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGCCCC'
complexity = 3/(51-1) = 6%
```
The threshold for low complexity filter can be specified by `-Y` or `--complexity_threshold`. It's range should be `0~100`, and its default value is 30, which means 30% complexity is required.

## Other filter
New filters are being implemented. If you have a new idea or new request, please file an issue.

# adapters
Adapter trimming is enabled by default, but you can disable it by `-A` or `--disable_adapter_trimming`. Adapter sequences can be automatically detected for both PE/SE data.
* For SE data, the adapter can be detected for most cases, but if `fastp` failed to find the adapter sequence, you can specify it by `-a` or `--adapter_sequence` option. If adapter sequence is specified, the auto detection is disabled.
* For PE data, the adapters can be trimmed automatically by per-read overlap analysis, which seeks for the overlap of each pair of reads. This method is robust and fast, so normally you don't have to input the adapter sequence. But you can still specify the adapter sequences for read1 by `--adapter_sequence`, and for read2 by `--adapter_sequence_r2`. In case `fastp` fails to find an overlap for some pairs (i.e. due to low quality bases), it will use these sequences to trim adapters for read1 and read2 respectively.
* For PE data, the auto adapter detection is disabled by default. You can enable it by specifing `-2` or `--detect_adapter_for_pe`. If you want to obtain ultra-clean data, this option is strongly suggested.
* For PE data, `fastp` will run a little slower if you specify the sequence adapters or enable the adapter auto-detection. But it may result in a slightly cleaner output (usually finds 0.1% to 0.5% more adapters), since the overlap analysis may fail due to sequencing errors.
* For PE data, you can specify `--allow_gap_overlap_trimming` to allow up to one gap when trim adapters by overlap analysis for PE data. By default no gap is allowed. This may take more time and usually have very limited effect (finds ~0.01% more adapters).
* The most widely used adapters are Illumina TruSeq adapters. If your data is from the TruSeq library, `fastp` should be able to detect it successfully, otherwise you can add `--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` to your command lines.
* `fastp` contains some built-in known adapter sequences for better auto-detection. If you want to make some adapters to be a part of the built-in adapters, please file an issue or make a change in https://github.com/OpenGene/fastp/blob/master/src/knownadapters.h

You can also specify `--adapter_fasta` to give a FASTA file to tell `fastp` to trim multiple adapters in this FASTA file. Here is a sample of such adapter FASTA file:
```
>Illumina TruSeq Adapter Read 1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina TruSeq Adapter Read 2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>polyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

The adapter sequence in this file should be at least 6bp long, otherwise it will be skipped. And you can give whatever you want to trim, rather than regular sequencing adapters (i.e. polyA).

`fastp` first trims the auto-detected adapter or the adapter sequences given by `--adapter_sequence | --adapter_sequence_r2`, then trims the adapters given by `--adapter_fasta` one by one.

The sequence distribution of trimmed adapters can be found at the HTML/JSON reports.

# per read cutting by quality score
`fastp` supports per read sliding window cutting by evaluating the mean quality scores in the sliding window. From `v0.19.6`, `fastp` supports 3 different operations, and you enable one or all of them:
* `-5, --cut_front`             move a sliding window from front (5') to tail, drop the bases in the window if its mean quality is below cut_mean_quality, stop otherwise. Default is disabled. The leading N bases are also trimmed. Use `cut_front_window_size` to set the widnow size, and `cut_front_mean_quality` to set the mean quality threshold. If the window size is 1, this is similar as the Trimmomatic `LEADING` method.
* `-3, --cut_tail`              move a sliding window from tail (3') to front, drop the bases in the window if its mean quality is below cut_mean_quality, stop otherwise. Default is disabled. The trailing N bases are also trimmed. Use `cut_tail_window_size` to set the widnow size, and `cut_tail_mean_quality` to set the mean quality threshold. If the window size is 1, this is similar as the Trimmomatic `TRAILING` method.
* `-r, --cut_right`             move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop. Use `cut_right_window_size` to set the widnow size, and `cut_right_mean_quality` to set the mean quality threshold.  This is similar as the Trimmomatic `SLIDINGWINDOW` method.


***WARNING: all these three operations will interfere deduplication for SE data, and `--cut_front` or `--cut_right` may also interfere deduplication for PE data. The deduplication algorithms rely on the exact matchment of coordination regions of the grouped reads/pairs.***

If `--cut_right` is enabled, then there is no need to enable `--cut_tail`, since the former is more aggressive. If `--cut_right` is enabled together with `--cut_front`, `--cut_front` will be performed first before `--cut_right` to avoid dropping whole reads due to the low quality starting bases.

Please be noted that `--cut_front` will interfere deduplication for both PE/SE data, and `--cut_tail` will interfere deduplication for SE data, since the deduplication algorithms rely on the exact matchment of coordination regions of the grouped reads/pairs.

If you don't set window size and mean quality threshold for these function respectively, `fastp` will use the values from `-W, --cut_window_size` and `-M, --cut_mean_quality `

# base correction for PE data
`fastp` perform `overlap analysis` for PE data, which try to find an overlap of each pair of reads. If an proper overlap is found, it can correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base will be assigned to it so that they will share the same quality.   

This function is not enabled by default, specify `-c` or `--correction` to enable it. This function is based on overlapping detection, which has adjustable parameters `overlap_len_require (default 30)`, `overlap_diff_limit (default 5)` and `overlap_diff_percent_limit (default 20%)`. Please note that the reads should meet these three conditions simultaneously.

# global trimming
`fastp` supports global trimming, which means trim all reads in the front or the tail. This function is useful since sometimes you want to drop some cycles of a sequencing run.

For example, the last cycle of Illumina sequencing is uaually with low quality, and it can be dropped with `-t 1` or `--trim_tail1=1` option.

* For read1 or SE data, the front/tail trimming settings are given with `-f, --trim_front1` and `-t, --trim_tail1`.
* For read2 of PE data, the front/tail trimming settings are given with `-F, --trim_front2` and `-T, --trim_tail2`. But if these options are not specified, they will be as same as read1 options, which means `trim_front2 = trim_front1` and `trim_tail2 = trim_tail1`.
* If you want to trim the reads to maximum length, you can specify `-b, --max_len1` for read1, and `-B, --max_len2` for read2. If `--max_len1` is specified but `--max_len2` is not, `--max_len2` will be same as `--max_len1`. For example, if `--max_len1` is specified and read1 is longer than `--max_len1`, `fastp` will trim read1 at its tail to make it as long as `--max_len1`.

Please note that the trimming for `--max_len` limitation will be applied at the last step. Following are fastp's processing steps that may orderly affect the read lengthes:
```
1, UMI preprocessing (--umi)
2, global trimming at front (--trim_front)
3, global trimming at tail (--trim_tail)
4, quality pruning at 5' (--cut_front)
5, quality pruning by sliding window (--cut_right)
6, quality pruning at 3' (--cut_tail)
7, trim polyG (--trim_poly_g, enabled by default for NovaSeq/NextSeq data)
8, trim adapter by overlap analysis (enabled by default for PE data)
9, trim adapter by adapter sequence (--adapter_sequence, --adapter_sequence_r2. For PE data, this step is skipped if last step succeeded)
10, trim polyX (--trim_poly_x)
11, trim to max length (---max_len)
```

# polyG tail trimming
For Illumina NextSeq/NovaSeq data, `polyG` can happen in read tails since `G` means no signal in the Illumina two-color systems. `fastp` can detect the polyG in read tails and trim them. This feature is enabled for NextSeq/NovaSeq data by default, and you can specify `-g` or `--trim_poly_g` to enable it for any data, or specify `-G` or `--disable_trim_poly_g` to disable it. NextSeq/NovaSeq data is detected by the machine ID in the FASTQ records.  

A minimum length can be set with `<poly_g_min_len>` for `fastp` to detect polyG. This value is 10 by default.

# polyX tail trimming
This feature is similar as polyG tail trimming, but is disabled by default. Use `-x` or `--trim_poly_x` to enable it. A minimum length can be set with `<poly_x_min_len>` for `fastp` to detect polyX. This value is 10 by default.

When `polyG tail trimming` and `polyX tail trimming` are both enabled, fastp will perform `polyG trimming` first, then perform `polyX trimming`. This setting is useful for trimming the tails having `polyX (i.e. polyA) ` before `polyG`. `polyG` is usually caused by sequencing artifacts, while `polyA` can be commonly found from the tails of mRNA-Seq reads.

# unique molecular identifier (UMI) processing
UMI is useful for duplication elimination and error correction based on generating consensus of reads originated from a same DNA fragment. It's usually used in deep sequencing applications like ctDNA sequencing. Commonly for Illumina platforms, UMIs can be integrated in two different places: `index` or head of `read`.  
To enable UMI processing, you have to enable `-U` or `--umi` option in the command line, and specify `--umi_loc`  to specify the UMI location, it can be one of:
* `index1` the first index is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
* `index2` the second index is used as UMI. PE data only, this UMI will be used for both read1/read2.
* `read1` the head of read1 is used as UMI. If the data is PE, this UMI will be used for both read1/read2.
* `read2` the head of read2 is used as UMI. PE data only, this UMI will be used for both read1/read2.
* `per_index` `index1_index2` is used as UMI for both read1/read2.
* `per_read` define `umi1` as the head of read1, and `umi2` as the head of read2. `umi1_umi2` is used as UMI for both read1/read2.

If `--umi_loc` is specified with `read1`, `read2` or `per_read`, the length of UMI should specified with `--umi_len`.

`fastp` will extract the UMIs, and append them to the first part of read names, so the UMIs will also be presented in SAM/BAM records. If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.

A prefix can be specified with `--umi_prefix`. If prefix is specified, an underline will be used to connect it and UMI. For example, UMI=AATTCCGG, prefix=UMI, then the final string presented in the name will be `UMI_AATTCCGG`.

If the UMI location is read1/read2/per_read, fastp can skip some bases after UMI to trim the UMI separator and A/T tailing. Specify `--umi_skip` to enable the number of bases to skip. By default it is not enabled.

## UMI example
The original read:
```
@NS500713:64:HFKJJBGXY:1:11101:1675:1101 1:N:0:TATAGCCT+GACCCCCA
AAAAAAAAGCTACTTGGAGTACCAATAATAAAGTGAGCCCACCTTCCTGGTACCCAGACATTTCAGGAGGTCGGGAAA
+
6AAAAAEEEEE/E/EA/E/AEA6EE//AEE66/AAE//EEE/E//E/AA/EEE/A/AEE/EEA//EEEEEEEE6EEAA
```
After it's processed with command: `fastp -i R1.fq -o out.R1.fq -U --umi_loc=read1 --umi_len=8`:  
```
@NS500713:64:HFKJJBGXY:1:11101:1675:1101:AAAAAAAA 1:N:0:TATAGCCT+GACCCCCA
GCTACTTGGAGTACCAATAATAAAGTGAGCCCACCTTCCTGGTACCCAGACATTTCAGGAGGTCGGGAAA
+
EEE/E/EA/E/AEA6EE//AEE66/AAE//EEE/E//E/AA/EEE/A/AEE/EEA//EEEEEEEE6EEAA
```

# output splitting
For parallel processing of FASTQ files (i.e. alignment in parallel), `fastp` supports splitting the output into multiple files. The splitting can work with two different modes: `by limiting file number` or `by limiting lines of each file`. These two modes cannot be enabled together.   

The file names of these split files will have a sequential number prefix, adding to the original file name specified by `--out1` or `--out2`, and the width of the prefix is controlled by the `-d` or `--split_prefix_digits` option. For example, `--split_prefix_digits=4`, `--out1=out.fq`, `--split=3`, then the output files will be `0001.out.fq`,`0002.out.fq`,`0003.out.fq`

## splitting by limiting file number
Use `-s` or `--split` to specify how many files you want to have. `fastp` evaluates the read number of a FASTQ by reading its first ~1M reads. This evaluation is not accurate so the file sizes of the last several files can be a little differnt (a bit bigger or smaller). For best performance, it is suggested to specify the file number to be a multiple of the thread number.

## splitting by limiting the lines of each file
Use `-S` or `--split_by_lines` to limit the lines of each file. The last files may have smaller sizes since usually the input file cannot be perfectly divided. The actual file lines may be a little greater than the value specified by `--split_by_lines` since `fastp` reads and writes data by blocks (a block = 1000 reads).

# overrepresented sequence analysis
Overrepresented sequence analysis is disabled by default, you can specify `-p` or `--overrepresentation_analysis` to enable it. For consideration of speed and memory, `fastp` only counts sequences with length of 10bp, 20bp, 40bp, 100bp or (cycles - 2 ).

By default, fastp uses 1/20 reads for sequence counting, and you can change this settings by specifying `-P` or `--overrepresentation_sampling` option. For example, if you set `-P 100`, only 1/100 reads will be used for counting, and if you set `-P 1`, all reads will be used but it will be extremely slow. The default value 20 is a balance of speed and accuracy.

`fastp` not only gives the counts of overrepresented sequence, but also gives the information that how they distribute over cycles. A figure is provided for each detected overrepresented sequence, from which you can know where this sequence is mostly found.

# merge paired-end reads
For paired-end (PE) input, fastp supports stiching them by specifying the `-m/--merge` option. In this `merging` mode:

* `--merged_out` shouuld be given to specify the file to store merged reads, otherwise you should enable `--stdout` to stream the merged reads to STDOUT. The merged reads are also filtered.
* `--out1` and `--out2` will be the reads that cannot be merged successfully, but both pass all the filters.
* `--unpaired1` will be the reads that cannot be merged, `read1` passes filters but `read2` doesn't.
* `--unpaired2` will be the reads that cannot be merged, `read2` passes filters but `read1` doesn't.
* `--include_unmerged` can be enabled to make reads of `--out1`, `--out2`, `--unpaired1` and `--unpaired2` redirected to `--merged_out`. So you will get a single output file. This option is disabled by default.

`--failed_out` can still be given to store the reads (either merged or unmerged) failed to passing filters.

In the output file, a tag like `merged_xxx_yyy`will be added to each read name to indicate that how many base pairs are from read1 and from read2, respectively. For example, `
@NB551106:9:H5Y5GBGX2:1:22306:18653:13119 1:N:0:GATCAG merged_150_15`
means that 150bp are from read1, and 15bp are from read2. `fastp` prefers the bases in read1 since they usually have higher quality than read2.

Same as the [base correction feature](#base-correction-for-pe-data), this function is also based on overlapping detection, which has adjustable parameters `overlap_len_require (default 30)`, `overlap_diff_limit (default 5)` and `overlap_diff_percent_limit (default 20%)`. Please note that the reads should meet these three conditions simultaneously.

# duplication rate and deduplication
For both SE and PE data, fastp supports evaluating its duplication rate and removing duplicated reads/pairs. fastp considers one read as duplicated only if its all base pairs are identical as another one. This means if there is a sequencing error or an N base, the read will not be treated as duplicated.

## duplication rate evaluation
By default, fastp evaluates duplication rate, and this module may use 1G memory and take 10% ~ 20% more running time. If you don't need the duplication rate information, you can set `--dont_eval_duplication` to disable the duplication evaluation. But please be noted that, if deduplication (`--dedup`) option is enabled, then `--dont_eval_duplication` option is ignored.

fastp uses a hash algorithm to find the identical sequences. Due to the possible hash collision, about 0.01% of the total reads may be wrongly recognized as deduplicated reads. Normally this may not impact the downstream analysis. The accuracy of calculating duplication can be improved by increasing the hash buffer number or enlarge the buffer size. The option `--dup_calc_accuracy` can be used to specify the level (1 ~ 6). The higher level means more memory usage and more running time. Please refer to following table:

| dup_calc_accuracy level | hash  buffer number | buffer size | memory usage | speed | |
|- | - | - | - | - | - |
| 1 | 1 | 1G | 1G | ultra-fast | default for no-dedup mode |
| 2 | 1 | 2G | 2G | fast | |
| 3 | 2 | 2G | 4G | fast | default for dedup|
| 4 | 2 | 4G | 8G | fast | |
| 5 | 2 | 8G | 12G | fast | |
| 6 | 3 | 8G | 24G | moderate | |

## deduplication
Since `v0.22.0`, fastp supports deduplication for FASTQ data. Specify `-D` or `--dedup` to enable this option. When `--dedup` is enabled, the `dup_calc_accuracy` level is default to `3`, and it can be changed to any value of 1 ~ 6.

# batch processing
[parallel.py](https://github.com/OpenGene/fastp/blob/master/parallel.py) is a script to preprocess all FASTQ files within a folder in parallel. It will automatically couple the paired-end FASTQ files.  

This script will generate an `overall.html` to present an aggregate summary for all processed FASTQ files.  

## example
```shell
python parallel.py -i /path/to/input/folder -o /path/to/output/folder -r /path/to/reports/folder -a '-f 3 -t 2'
```
which means to  
```
. process all the FASTQ data in /path/to/input/folder
. using fastp in PATH
. with arguments -f 3 and -t 2, which means trimming 3bp in head and 2bp in tail
. output all clean data to /path/to/output/folder
. output all HTML and JSON reports to /path/to/reports/folder
```

See `python parallel.py -h` for details.

# all options
```shell
usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
options:
  # I/O options
  -i, --in1                          read1 input file name (string)
  -o, --out1                         read1 output file name (string [=])
  -I, --in2                          read2 input file name (string [=])
  -O, --out2                           read2 output file name (string [=])
      --unpaired1                      for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. (string [=])
      --unpaired2                      for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
      --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
      --overlapped_out                 for each read pair, output the overlapped region if it has no any mismatched base. (string [=])
  -m, --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
      --merged_out                     in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output (string [=])
      --include_unmerged               in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. Disabled by default.
  -6, --phred64                      indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
      --stdin                          input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
      --stdout                         output passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end input. Disabled by default.
      --interleaved_in                 indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
      --reads_to_process             specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
      --dont_overwrite               don't overwrite existing files. Overwritting is allowed by default.
      --fix_mgi_id                     the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.

  # adapter trimming options
  -A, --disable_adapter_trimming     adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
      --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=])
      --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
  -2, --detect_adapter_for_pe          enable adapter detection for PE data to get ultra-clean data. It takes more time to find just a little bit more adapters.
      --allow_gap_overlap_trimming     allow up to one gap when trim adapters by overlap analysis for PE data. By default no gap is allowed.

  # global trimming options
  -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
  -b, --max_len1                       if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
  -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
  -B, --max_len2                       if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])

  # duplication evaluation and deduplication
  -D, --dedup                          enable deduplication to drop the duplicated reads/pairs
      --dup_calc_accuracy              accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
      --dont_eval_duplication          don't evaluate duplication rate to save time and use less memory.

  # polyG tail trimming, useful for NextSeq/NovaSeq data
  -g, --trim_poly_g                  force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
      --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
  -G, --disable_trim_poly_g          disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data

  # polyX tail trimming
  -x, --trim_poly_x                    enable polyX trimming in 3' ends.
      --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])

  # per read cutting by quality options
  -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -r, --cut_right                      move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
  -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
  -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
      --cut_front_window_size          the window size option of cut_front, default to cut_window_size if not specified (int [=4])
      --cut_front_mean_quality         the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
      --cut_tail_window_size           the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
      --cut_tail_mean_quality          the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
      --cut_right_window_size          the window size option of cut_right, default to cut_window_size if not specified (int [=4])
      --cut_right_mean_quality         the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])

  # quality filtering options
  -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is specified, quality filtering is disabled
  -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
  -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
  -e, --average_qual                 if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])


  # length filtering options
  -L, --disable_length_filtering     length filtering is enabled by default. If this option is specified, length filtering is disabled
  -l, --length_required              reads shorter than length_required will be discarded, default is 15. (int [=15])
      --length_limit                 reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])

  # low complexity filtering
  -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
  -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])

  # filter reads with unwanted indexes (to remove possible contamination)
      --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])

  # base correction by overlap analysis options
  -c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled
      --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
      --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
      --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])

  # UMI processing
  -U, --umi                          enable unique molecular identifier (UMI) preprocessing
      --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])
      --umi_prefix                   if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])

  # overrepresented sequence analysis
  -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
  -P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])

  # reporting options
  -j, --json                         the json format report file name (string [=fastp.json])
  -h, --html                         the html format report file name (string [=fastp.html])
  -R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])

  # threading options
  -w, --thread                       worker thread number, default is 3 (int [=3])

  # output splitting options
  -s, --split                        split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
  -S, --split_by_lines               split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
  -d, --split_prefix_digits          the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])

  # help
  -?, --help                         print this message
```

# citations
### Shifu Chen. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. https://doi.org/10.1002/imt2.107
### Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560



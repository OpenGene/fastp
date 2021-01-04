[![install with conda](
https://anaconda.org/bioconda/fastp/badges/version.svg)](https://anaconda.org/bioconda/fastp)
[![install with conda](
https://anaconda.org/bioconda/fastp/badges/downloads.svg)](https://anaconda.org/bioconda/fastp)
[![Build Status](https://travis-ci.org/OpenGene/fastp.svg?branch=master)](https://travis-ci.org/OpenGene/fastp)
# fastp
A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. 
* [features](#features)
* [simple usage](#simple-usage)
* [examples of report](#examples-of-report)
* [download, compile and install](#get-fastp)
* [input and output](#input-and-output)
* [filtering by quality, length, complexity, etc.](#filtering)
* [adapter trimming](#adapters)
* [per read cutting by quality score](#per-read-cutting-by-quality-score)
* [base correction for paired end (PE) data](#base-correction-for-pe-data)
* [globa trimming](#global-trimming)
* [polyG tail trimming](#polyg-tail-trimming) and [polyX tail trimming](#polyx-tail-trimming)
* [unique molecular identifer (UMI) processing](#unique-molecular-identifer-umi-processing)
* [output splitting](#output-splitting)
* [overrepresented sequence analysis](#overrepresented-sequence-analysis)
* [all options](#all-options)
* [citation](#citation)

# features
0. comprehensive quality profiling for both before and after filtering data (quality curves, base contents, KMER, Q20/Q30, GC Ratio, duplication, adapter contents...)
1. filter out bad reads (too low quality, too short, or too many N...)
2. cut low quality bases for per read in its 5' and 3' by evaluating the mean quality from a sliding window (like Trimmomatic but faster).
3. trim all reads in front and tail
4. cut adapters. Adapter sequences can be automatically detected,which means you don't have to input the adapter sequences to trim them.
5. correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
6. trim polyG in 3' ends, which is commonly seen in NovaSeq/NextSeq data. Trim polyX in 3' ends to remove unwanted polyX tailing (i.e. polyA tailing for mRNA-Seq data)
7. preprocess unique molecular identifer (UMI) enabled data, shift UMI to sequence name.
8. report JSON format result for further interpreting. 
9. visualize quality control and filtering results on a single HTML page (like FASTQC but faster and more informative).
10. split the output to multiple files (0001.R1.gz, 0002.R1.gz...) to support parallel processing. Two modes can be used, limiting the total split file number, or limitting the lines of each split file.
11. support long reads (data from PacBio / Nanopore devices).
12. support reading from STDIN and writing to STDOUT
13. support interleaved input
14. ...

This tool is being intensively developed, and new features can be implemented soon if they are considered useful. If you have any additional requirement for `fastp`, please file an issue:https://github.com/OpenGene/fastp/issues/new

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
## or download binary (only for Linux systems, http://opengene.org/fastp/fastp)
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```
## or compile from source
```shell
# get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/fastp.git

# build
cd fastp
make

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
## process only part of the data
If you don't want to process all the data, you can specify `--reads_to_process` to limit the reads to be processed. This is useful if you want to have a fast preview of the data quality, or you want to create a subset of the filtered data.
## do not overwrite exiting files
You can enable the option `--dont_overwrite` to protect the existing files not to be overwritten by `fastp`. In this case, `fastp` will report an error and quit if it finds any of the output files (read1, read2, json report, html report) already exists before.
## split the output to multiple files for parallel processing
See [output splitting](#output-splitting)

# filtering
Multiple filters have been implemented.
## quality filter
Quality filtering is enabled by default, but you can disable it by `-Q` or `disable_quality_filtering`. Currently it supports filtering by limiting the N base number (`-n, --n_base_limit`),  and the percentage of unqualified bases.  

To filter reads by its percentage of unqualified bases, two options should be provided:
* `-q, --qualified_quality_phred`       the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. 
* `-u, --unqualified_percent_limit`    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%

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
Adapter trimming is enabled by default, but you can disable it by `-A` or `--disable_adapter_trimming`. Adapter sequences can be automatically detected for both PE/SE data, which means you don't have to input the adapter sequences to trim them.
* For SE data, the adapters are evaluated by analyzing the tails of first ~1M reads. This evaluation may be inacurrate, and you can specify the adapter sequence by `-a` or `--adapter_sequence` option. If adapter sequence is specified, the auto detection for SE data will be disabled.
* For PE data, the adapters can be detected by per-read overlap analysis, which seeks for the overlap of each pair of reads. This method is robust and fast, so normally you don't have to input the adapter sequence even you know it. But you can still specify the adapter sequences for read1 by `--adapter_sequence`, and for read2 by `--adapter_sequence_r2`. If `fastp` fails to find an overlap (i.e. due to low quality bases), it will use these sequences to trim adapters for read1 and read2 respectively.

The sequence distribution of trimmed adapters can be found at the HTML/JSON reports.

# per read cutting by quality score
`fastp` supports per read sliding window cutting by evaluate the mean quality scores in the sliding window, which is similar with how `Trimmomatic` does, but `fastp` is much faster. This function is disabled by default, to enable it, you can specify either or both of:
* `-5, --cut_by_quality5`              enable per read cutting by quality in front (5'), and trim leading N bases.
* `-3, --cut_by_quality3`              enable per read cutting by quality in tail (3'), and trim trailing N bases.

Please be noted that `--cut_by_quality5` will interfere deduplication for both PE/SE data, and `--cut_by_quality3` will interfere deduplication for SE data, since the deduplication algorithms rely on the exact matchment of coordination regions of the grouped reads/pairs.

The size of sliding window can be specified with `-W, --cut_window_size`, and the mean quality requirement can be specified with `-M, --cut_mean_quality `.

# base correction for PE data
`fastp` perform `overlap analysis` for PE data, which try to find an overlap of each pair of reads. If an proper overlap is found, it can correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base will be assigned to it so that they will share the same quality.   

This function is not enabled by default, specify `-c` or `--correction` to enable it.

# global trimming
`fastp` supports global trimming, which means trim all reads in the front or the tail. This function is useful since sometimes you want to drop some cycles of a sequencing run.

For example, the last cycle of Illumina sequencing is uaually with low quality, and it can be dropped with `-t 1` or `--trim_tail1=1` option.

* For read1 or SE data, the front/tail trimming settings are given with `-f, --trim_front1` and `-t, --trim_tail1`.
* For read2 of PE data, the front/tail trimming settings are given with `-F, --trim_front2` and `-T, --trim_tail2`. But if these options are not specified, they will be as same as read1 options, which means `trim_front2 = trim_front1` and `trim_tail2 = trim_tail1`.

# polyG tail trimming
For Illumina NextSeq/NovaSeq data, `polyG` can happen in read tails since `G` means no signal in the Illumina two-color systems. `fastp` can detect the polyG in read tails and trim them. This feature is enabled for NextSeq/NovaSeq data by default, and you can specify `-g` or `--trim_poly_g` to enable it for any data, or specify `-G` or `--disable_trim_poly_g` to disable it. NextSeq/NovaSeq data is detected by the machine ID in the FASTQ records.   

A minimum length can be set with `<poly_g_min_len>` for `fastp` to detect polyG. This value is 10 by default.

# polyX tail trimming
This feature is similar as polyG tail trimming, but is disabled by default. Use `-x` or `--polyX` to enable it. A minimum length can be set with `<poly_x_min_len>` for `fastp` to detect polyX. This value is 10 by default.   

When `polyG tail trimming` and `polyX tail trimming` are both enabled, fastp will perform `polyG trimming` first, then perform `polyX trimming`. This setting is useful for trimming the tails having `polyX (i.e. polyA) ` before `polyG`. `polyG` is usually caused by sequencing artifacts, while `polyA` can be commonly found from the tails of mRNA-Seq reads.

# unique molecular identifer (UMI) processing
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

# all options
```shell
usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
options:
  # I/O options
  -i, --in1                          read1 input file name (string)
  -o, --out1                         read1 output file name (string [=])
  -I, --in2                          read2 input file name (string [=])
  -O, --out2                         read2 output file name (string [=])
  -6, --phred64                      indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
      --stdin                          input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
      --stdout                         output passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end input. Disabled by defaut.
      --interleaved_in                 indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by defaut.
      --reads_to_process             specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
      --dont_overwrite               don't overwrite existing files. Overwritting is allowed by default.

  # adapter trimming options
  -A, --disable_adapter_trimming     adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence             the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
      --adapter_sequence_r2          the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=])
      --adapter_mm_freq              allowed mismatched within every n bases to detect adapter. 8 by defaults. (int [=8])

  # global trimming options
  -f, --trim_front1                  trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                   trimming how many bases in tail for read1, default is 0 (int [=0])
  -F, --trim_front2                  trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                   trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])

  # polyG tail trimming, useful for NextSeq/NovaSeq data
  -g, --trim_poly_g                  force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
      --poly_g_min_len               the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
      --poly_g_mm_freq               allowed mismatched within every n bases to detect polyG. 8 by defaults. (int [=8])
      --poly_g_mm_max                the maximum number of mismatched allowed in detecting polyG. 5 by defaults. (int [=5])
  -G, --disable_trim_poly_g          disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data

  # polyX tail trimming
  -x, --trim_poly_x                  enable polyX trimming in 3' ends.
      --poly_x_min_len               the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
      --poly_x_mm_freq               allowed mismatched within every n bases to detect polyX. 8 by defaults. (int [=8])
      --poly_x_mm_max                the maximum number of mismatched allowed in detecting polyX. 5 by defaults. (int [=5])

  # per read cutting by quality options
  -5, --cut_by_quality5              enable per read cutting by quality in front (5'), default is disabled (WARNING: this will interfere deduplication for both PE/SE data)
  -3, --cut_by_quality3              enable per read cutting by quality in tail (3'), default is disabled (WARNING: this will interfere deduplication for SE data)
  -W, --cut_window_size              the size of the sliding window for sliding window trimming, default is 4 (int [=4])
  -M, --cut_mean_quality             the bases in the sliding window with mean quality below cutting_quality will be cut, default is Q20 (int [=20])

  # quality filtering options
  -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is specified, quality filtering is disabled
  -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
  -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])

  # length filtering options
  -L, --disable_length_filtering     length filtering is enabled by default. If this option is specified, length filtering is disabled
  -l, --length_required              reads shorter than length_required will be discarded, default is 15. (int [=15])
      --length_limit                 reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])

  # low complexity filtering
  -y, --low_complexity_filter        enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
  -Y, --complexity_threshold         the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])

  # filter reads with unwanted indexes (to remove possible contamination)
      --filter_by_index1             specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      --filter_by_index2             specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      --filter_by_index_threshold    the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])

  # base correction by overlap analysis options
  -c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled
      --overlap_len_require          the minimum length of the overlapped region for overlap analysis based adapter trimming and correction. 30 by default. (int [=30])
      --overlap_diff_limit           the maximum difference of the overlapped region for overlap analysis based adapter trimming and correction. 5 by default. (int [=5])

  # UMI processing
  -U, --umi                          enable unique molecular identifer (UMI) preprocessing
      --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])
      --umi_prefix                   if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      --umi_skip                     if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])

  # overrepresented sequence analysis
  -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
  -P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])

  # reporting options
  -j, --json                         the json format report file name (string [=fastp.json])
  -h, --html                         the html format report file name (string [=fastp.html])
  -R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])

  # threading options
  -w, --thread                       worker thread number, default is 2 (int [=2])

  # output splitting options
  -s, --split                        split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
  -S, --split_by_lines               split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
  -d, --split_prefix_digits          the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])

  # help
  -?, --help                         print this message
```

# citation
`fastp` paper has been accepted by ECCB 2018 conference and will be published in `Bioinformatics` journal. Before it's final publication, you can cite fastp's bioRxiv paper:  

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu. fastp: an ultra-fast all-in-one FASTQ preprocessor. bioRxiv 274100; doi: https://doi.org/10.1101/274100



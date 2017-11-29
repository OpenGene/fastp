# fastp
A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. 
* [features](#features)
* [simple usage](#simple-usage)
* [examples of report](#examples-of-report)
* [download, compile and install](#get-fastp)
* [adapter trimming](#adapters)
* [per read cutting by quality score](#per-read-cutting-by-quality-score)
* [base correction for paired end (PE) data](#base-correction-for-pe-data)
* [globa trimming](#global-trimming)
* [unique molecular identifer (UMI) processing](#unique-molecular-identifer-umi-processing)
* [output splitting](#output-splitting)
* [all options](#all-options)
* [citation](#citation)

# features
1. filter out bad reads (too low quality, too short, or too many N...)
2. cut low quality bases for per read in its 5' and 3' by evaluating the mean quality from a sliding window (like Trimmomatic but faster).
3. trim all reads in front and tail
4. cut adapters. Adapter sequences can be automatically detected,which means you don't have to input the adapter sequences to trim them.
5. correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality
6. preprocess unique molecular identifer (UMI) enabled data, shift UMI to sequence name.
7. report JSON format result for further interpreting. 
8. visualize quality control and filtering results on a single HTML page (like FASTQC but faster and more informative).
9. split the output to multiple files (0001.R1.gz, 0002.R1.gz...) to support parallel processing. Two modes can be used, limiting the total split file number, or limitting the lines of each split file.
10. support long reads (data from PacBio / Nanopore devices).
11. ...

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
## download binary (only for Linux systems, http://opengene.org/fastp/fastp)
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
```
## compile from source
```shell
# get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/fastp.git

# build
cd fastp
make

# Install
sudo make install
```

# filtering
Quality filtering is enabled by default, but you can disable it by `-Q` or `disable_quality_filtering`. Currently it supports filtering by limiting the N base number (`-n, --n_base_limit`),  and the percentage of unqualified bases.  

To filter reads by its percentage of unqualified bases, two options should be provided:
* `-q, --qualified_quality_phred`       the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. 
* `-u, --unqualified_percent_limit`    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%

Length filtering is enabled by default, but you can disable it by `-L` or `--disable_length_filtering`. The minimum length requirement is specified with `-l` or `--length_required`.

New filters are being implemented, such like `polyX` filter and `low complexity` filter. If you have a new idea or new request, please file an issue.

# adapters
Adapter trimming is enabled by default, but you can disable it by `-A` or `--disable_adapter_trimming`. Adapter sequences can be automatically detected for both PE/SE data, which means you don't have to input the adapter sequences to trim them.
* For SE data, the adapters are evaluated by analyzing the tails of first ~1M reads. This evaluation may be inacurrate, and you can specify the adapter sequence by `-a` or `--adapter_sequence` option. If adapter sequence is specified, the auto detection for SE data will be disabled.
* For PE data, the adapters are detected by overlap analysis, which seeks for the overlap of each pair of reads. This method is robust and fast so that you don't have to input the adapter sequence even you know it. For PE data, the `-a` or `--adapter_sequence` option is ignored.  

The sequence distribution of trimmed adapters can be found at the HTML/JSON reports.

# per read cutting by quality score
`fastp` supports per read sliding window cutting by evaluate the mean quality scores in the sliding window, which is similar with how `Trimmomatic` does, but `fastp` is much faster. This function is disabled by default, to enable it, you can specify either or both of:
* `-5, --cut_by_quality5`              enable per read cutting by quality in front (5')
* `-3, --cut_by_quality3`              enable per read cutting by quality in tail (3')

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

# unique molecular identifer (UMI) processing
UMI is useful for duplication elimination and error correction based on generating consesus of reads originated from a same DNA fragment. It's usually used in deep sequencing applications like ctDNA sequencing. Commonly for Illumina platforms, UMIs can be integrated in two different places: `index` or head of `read`.  
To enable UMI processing, you have to enable `-U` or `--umi` option in the command line, and specify `--umi_loc`  to specify the UMI location, it can be one of:
* `index1` the first index is used as UMI.
* `index2` the second index is used as UMI. PE data only.
* `read1` the head of read1 is used as UMI.
* `read2` the head of read2 is used as UMI. PE data only.  

If `--umi_loc` is specifie dwith `read1` or `read2`, the length of UMI should specified with `--umi_len`. 

`fastp` will extract the UMIs, and append them to the first part of read names, so that it will also be presented in SAM/BAM records. A perfix `UMI_` will be added to each UMI. If the UMI is in the reads, then it will be shifted from read so that the read will become shorter. If the UMI is in the index, it will be kept.

## UMI example
original read:
```
@NS500713:64:HFKJJBGXY:1:11101:1675:1101 1:N:0:TATAGCCT+GACCCCCA
AAAAAAAAGCTACTTGGAGTACCAATAATAAAGTGAGCCCACCTTCCTGGTACCCAGACATTTCAGGAGGTCGGGAAA
+
6AAAAAEEEEE/E/EA/E/AEA6EE//AEE66/AAE//EEE/E//E/AA/EEE/A/AEE/EEA//EEEEEEEE6EEAA
```
processed with: `fastp -i testdata/R1.fq -o testdata/out.R1.fq -U --umi_loc=read1 --umi_len=8`:  
```
@NS500713:64:HFKJJBGXY:1:11101:1675:1101:UMI_AAAAAAAA 1:N:0:TATAGCCT+GACCCCCA
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

# all options
```shell
usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
options:
  # I/O options
  -i, --in1                          read1 input file name (string)
  -o, --out1                         read1 output file name (string [=])
  -I, --in2                          read2 input file name (string [=])
  -O, --out2                         read2 output file name (string [=])
  -6, --phred64                      indicates the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 2. (int [=2])
  
  # adapter trimming options
  -A, --disable_adapter_trimming     adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence             the adapter for SE data, default is auto (automatic detection). For PE data adapters can be trimmed without knowing the sequences. (string [=auto])
    
  # global trimming options
  -f, --trim_front1                  trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                   trimming how many bases in tail for read1, default is 0 (int [=0])
  -F, --trim_front2                  trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                   trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
  
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
  -l, --length_required              reads shorter than length_required will be discarded. (int [=30])

  # base correction by overlap analysis options
  -c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled
  
  # UMI processing
  -U, --umi                          enable unique molecular identifer (UMI) preprocessing
      --umi_loc                      specify the location of UMI, can be (index1/index2/read1/read2), default is none (string [=])
      --umi_len                      if the UMI is in read1/read2, its length should be provided (int [=0])

  # reporting options
  -j, --json                         the json format report file name (string [=fastp.json])
  -h, --html                         the html format report file name (string [=fastp.html])
  -R, --report_title                 should be quoted with ' or ", default is "fastp report" (string [=fastp report])
  
  # thread options
  -w, --thread                       worker thread number, default is 3 (int [=3])
  
  # output splitting options
  -s, --split                        split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
  -S, --split_by_lines               split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
  -d, --split_prefix_digits          the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
  
  # help
  -?, --help                         print this message
```

# citation
A paper for this tool is being written, and may be soon available in bioRxiv. If you want to cite this tool before it 
is completed, cite like this:
```
Shifu Chen, A fast FASTQ preprocessor with full features, (2017), GitHub repository, https://github.com/OpenGene/fastp
```


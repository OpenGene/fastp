# fastp
This tool is designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. It has following features:
* low quality reads filtering (too low quality, or too many N)
* adapter cutting automatically (works paired-end data)
* quality control with HTML visualization and JSON reports (HTML/JSON being implemented)
* possible error correction (being implemented )
* ...

This tool is being intensively developed, and new features can be implemented soon if they are considered useful. If you have any additional requirement for `fastp`, please file an issue:https://github.com/OpenGene/fastp/issues/new

# usage
```shell
shifus-MacBook-Pro:fastp shifu$ ./fastp 
usage: ./fastp --in1=string [options] ... 
options:
  -i, --in1                          read1 input file name (string)
  -o, --out1                         read1 output file name (string [=])
  -I, --in2                          read2 input file name (string [=])
  -O, --out2                         read2 output file name (string [=])
  -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 2. (int [=2])
  -f, --trim_front                   trimming how many bases in front, default is 0 (int [=0])
  -t, --trim_tail                    trimming how many bases in tail, default is 0 (int [=0])
  -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is enabled, quality filtering is disabled
  -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
  -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
  -l, --length_required              length filtering will be enabled if this argument is specified, reads shorter than length_required will be discarded. (int [=30])
  -j, --json                         the json format report file name (string [=fastp.json])
  -h, --html                         the html format report file name (string [=fastp.html])
  -T, --thread                       worker thread number, default is 3 (int [=3])
  -?, --help                         print this message
```

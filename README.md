# fastp
This tool is designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. It has following features:
* filter low quality reads (too low quality, or too many N)
* apply possible error correction (being implemented )
* cut adapters automatically (works for paired-end data)
* report JSON format result for further interpreting 
* Visualize quality profiles and filtering results on a single HTML page (being implemented)
* ...

This tool is being intensively developed, and new features can be implemented soon if they are considered useful. If you have any additional requirement for `fastp`, please file an issue:https://github.com/OpenGene/fastp/issues/new

# examples of report
`fastp` creates reports in both HTML and JSON format.
* HTML report: http://opengene.org/fastp/fastp.html
* JSON report: http://opengene.org/fastp/fastp.json

# Get fastp
## Download
Get latest
```shell
# download by http
https://github.com/OpenGene/fastp/archive/master.zip

# or clone by git
git clone https://github.com/OpenGene/fastp.git
```
Get the stable releases  
https://github.com/OpenGene/fastp/releases/latest

## Build
fastp only depends on `libz`, which is always available on Linux or Mac systems. If your system has no `libz`, install it first.
```shell
cd fastp
make
```

## Install
After build is done, run
```
sudo make install
```

# usage
```shell
usage: fastp -i <in1> -o <out1> [-I <in1> -O <out2>] [options...]
options:
  -i, --in1                          read1 input file name (string)
  -o, --out1                         read1 output file name (string [=])
  -I, --in2                          read2 input file name (string [=])
  -O, --out2                         read2 output file name (string [=])
  -z, --compression                  compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 2. (int [=2])
  -f, --trim_front1                  trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                   trimming how many bases in tail for read1, default is 0 (int [=0])
  -F, --trim_front2                  trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                   trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
  -Q, --disable_quality_filtering    quality filtering is enabled by default. If this option is enabled, quality filtering is disabled
  -q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
  -u, --unqualified_percent_limit    how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                 if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
  -l, --length_required              length filtering will be enabled if this argument is specified, reads shorter than length_required will be discarded. (int [=30])
  -j, --json                         the json format report file name (string [=fastp.json])
  -h, --html                         the html format report file name (string [=fastp.html])
  -w, --thread                       worker thread number, default is 3 (int [=3])
  -s, --split                        if this option is specified, the output will be split to multiple (--split) files (i.e. 0001.out.fq, 0002.out.fq...).  (int [=0])
  -d, --split_prefix_digits          the digits for the slice number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
  -?, --help                         print this message
```

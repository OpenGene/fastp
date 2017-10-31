# fastp
This tool will be designed to provide fast all-in-one preprocessing for FastQ files, with features:
* low quality reads filtering (too low quality, or too many N)
* adapter cutting (adapters will be detected automatically)
* quality control with HTML visualization and JSON reports.
* possible error correction
* reads reordering to achieve higher gzip compression ratio
* ...

If you have any additional requirement for `fastp`, please file an issue:https://github.com/OpenGene/fastp/issues/new

# Implementation
`fastp` will be written in C++, with multiple threading supported. It's being developed, and the first release is expected to be published before Nov 10.

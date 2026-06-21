# fastp Benchmark Suite

End-to-end benchmark for fastp: measures throughput, peak memory, and verifies output correctness across I/O format combinations and feature modes.

## Quick Start

```bash
# Run all I/O format modes with generated data (10M pairs)
python3 -m benchmark --opt ./fastp

# Compare optimized build against baseline
python3 -m benchmark --opt ./fastp_opt --orig ./fastp_baseline

# Use your own FASTQ files
python3 -m benchmark --r1 /data/sample_R1.fq.gz --r2 /data/sample_R2.fq.gz --opt ./fastp
```

## Modes

### I/O Format Modes (default: `all`)

Test every input/output compression combination:

| Mode       | Type | Input  | Output |
|------------|------|--------|--------|
| fq-fq      | PE   | .fq    | .fq    |
| fq-gz      | PE   | .fq    | .fq.gz |
| gz-fq      | PE   | .fq.gz | .fq    |
| gz-gz      | PE   | .fq.gz | .fq.gz |
| se-fq-fq   | SE   | .fq    | .fq    |
| se-fq-gz   | SE   | .fq    | .fq.gz |
| se-gz-fq   | SE   | .fq.gz | .fq    |
| se-gz-gz   | SE   | .fq.gz | .fq.gz |
| stdin-stdout| SE  | stdin  | stdout |

### Feature Modes (`all-feat`)

Exercise specific fastp processing paths:

| Mode       | fastp Flags                               | Data        | Tests                              |
|------------|-------------------------------------------|-------------|------------------------------------|
| merge      | `--merge --correction`                    | overlapping | PE merge + base correction (#676)  |
| correction | `--correction`                            | overlapping | overlap analysis + base correction |
| dedup-pe   | `--dedup`                                 | random      | duplicate detection (PE)           |
| dedup-se   | `--dedup`                                 | random      | duplicate detection (SE)           |
| qualcut-pe | `--cut_front --cut_tail -W 4 -M 20`      | random      | sliding window quality cut (PE)    |
| qualcut-se | `--cut_front --cut_tail -W 4 -M 20`      | random      | sliding window quality cut (SE)    |
| ht-pe      | (default) `-w 32`                         | random      | high-thread concurrency stress     |

### Mode Groups

| Alias      | Expands To                                     |
|------------|------------------------------------------------|
| all        | all 9 I/O format modes (default)               |
| all-pe     | PE format modes only                           |
| all-se     | SE format modes + stdin-stdout                 |
| all-feat   | all 7 feature modes                            |
| everything | all I/O + all feature modes (16 total)         |

```bash
python3 -m benchmark --mode all-feat --opt ./fastp
python3 -m benchmark --mode merge,correction --opt ./fastp --orig ./fastp_v110
python3 -m benchmark --mode everything --opt ./fastp
```

## Options

```
--opt PATH       Optimized binary path (default: /tmp/fastp_opt)
--orig PATH      Baseline binary for comparison (optional)
--r1 PATH        Custom R1 input (.fq or .fq.gz), skips data generation
--r2 PATH        Custom R2 input (.fq or .fq.gz), required for PE modes
--mode MODE      Comma-separated mode list (default: all)
--pairs N        Read pairs to generate (default: 10,000,000)
--threads N      Worker threads, 0 = auto (default: 0)
--runs N         Repeat count, reports median (default: 3)
--seed N         RNG seed for data generation (default: 2026)
--json PATH      Load saved JSON results and print summary table
--merge A B      Merge two opt-only result JSONs into comparison table
```

## Test Data

When `--r1`/`--r2` are not provided, the benchmark generates synthetic data:

- **Random PE data** (10M pairs x 150bp): independent R1/R2 with realistic quality profiles. Used for I/O format modes, dedup, qualcut, and high-thread tests.
- **Overlapping PE data** (10M pairs x 150bp): R1 and R2 derived from the same fragment (insert ~220bp, overlap ~80bp). R2 has injected mismatches at low quality (Q10-14) to trigger base correction. Used for merge and correction modes.

Generated data is cached in `/tmp/fastp_benchmark/data/` and reused across runs. Delete the directory to regenerate.

When custom files are provided:
- Format conversion is automatic (decompress .gz to .fq or compress .fq to .gz as needed)
- Merge/correction modes use the custom data directly (real sequencing data has natural overlap)

## Output

- Per-run timing and peak RSS printed live
- Markdown summary table at the end
- JSON results saved to `/tmp/fastp_benchmark/results.json`
- MD5 verification: when baseline is provided, checks that optimized output matches baseline byte-for-byte

## Module Structure

```
benchmark/
├── __main__.py    # python3 -m benchmark entry point
├── run.py         # CLI argument parsing and orchestration
├── modes.py       # mode definitions, constants, path config
├── datagen.py     # synthetic FASTQ data generation
├── sysinfo.py     # hardware/OS info collection
├── runner.py      # benchmark execution (build_cmd, run, measure)
├── verify.py      # output MD5 verification
└── report.py      # summary table formatting
```

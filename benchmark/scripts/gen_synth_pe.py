#!/usr/bin/env python3
"""Deterministic synthetic paired-end Illumina reads with Nextera adapter read-through.

Used as one corpus entry alongside real public accessions: gives a reproducible,
license-free stressor for adapter-detection code paths without depending on any
external download.
"""
import argparse
import gzip
import random

NEXTERA_R1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
NEXTERA_R2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGAT"
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def qual(rng, n):
    return "".join(chr(33 + max(2, min(40, int(rng.gauss(36, 4))))) for _ in range(n))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=500000)
    ap.add_argument("--len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=721)
    ap.add_argument("--out", default="synth")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    genome = "".join(rng.choice("ACGT") for _ in range(2_000_000))
    with gzip.open(f"{a.out}_R1.fastq.gz", "wt", compresslevel=1) as f1, \
         gzip.open(f"{a.out}_R2.fastq.gz", "wt", compresslevel=1) as f2:
        for i in range(a.pairs):
            ins = int(rng.lognormvariate(5.0, 0.6))
            ins = max(20, min(ins, 1000))
            pos = rng.randrange(0, len(genome) - ins)
            frag = genome[pos:pos + ins]
            r1 = (frag + NEXTERA_R1 + "A" * a.len)[:a.len]
            r2 = (revcomp(frag) + NEXTERA_R2 + "A" * a.len)[:a.len]
            name = f"@SYN:1:FC:1:{1 + i // 1000000}:{i % 100000}:{i}"
            f1.write(f"{name} 1:N:0:1\n{r1}\n+\n{qual(rng, a.len)}\n")
            f2.write(f"{name} 2:N:0:1\n{r2}\n+\n{qual(rng, a.len)}\n")


if __name__ == "__main__":
    main()

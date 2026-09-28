#!/usr/bin/env python3
"""Deterministic synthetic paired-end FASTQ for CI.

Models the features fastp's trimming and QC code paths react to: TruSeq adapter
read-through from short inserts, 2-color poly-G tails, N bases, and quality that
degrades along the read. Same arguments always produce byte-identical files.

usage: gen_reads.py OUT_PREFIX [--pairs N] [--len L] [--seed S]
writes OUT_PREFIX_R1.fastq.gz, OUT_PREFIX_R2.fastq.gz, OUT_PREFIX_interleaved.fastq.gz
"""
import argparse, gzip, random

ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGT", "TGCA")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--pairs", type=int, default=40000)
    ap.add_argument("--len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    L = a.len

    def qual():
        # high quality early, degrading and noisier towards the 3' end
        return "".join(chr(33 + max(2, min(41, int(38 - 14 * i / L + rng.gauss(0, 3))))) for i in range(L))

    def read(insert, adapter):
        s = insert[:L]
        if len(s) < L:
            s = (s + adapter + "".join(rng.choices("ACGT", k=L)))[:L]
        s = list(s)
        if rng.random() < 0.05:  # 2-color chemistry: signal loss reads as G
            start = rng.randrange(L // 2, L)
            s[start:] = "G" * (L - start)
        for _ in range(rng.choice([0, 0, 0, 1, 2])):
            s[rng.randrange(L)] = "N"
        return "".join(s)

    with gzip.open(a.out + "_R1.fastq.gz", "wt", compresslevel=1) as f1, \
         gzip.open(a.out + "_R2.fastq.gz", "wt", compresslevel=1) as f2, \
         gzip.open(a.out + "_interleaved.fastq.gz", "wt", compresslevel=1) as fi:
        for i in range(a.pairs):
            size = int(rng.lognormvariate(5.3, 0.45))  # median ~200bp, ~30% shorter than L
            frag = "".join(rng.choices("ACGT", k=max(20, size)))
            r1 = read(frag, ADAPTER_R1)
            r2 = read(frag.translate(COMP)[::-1], ADAPTER_R2)
            rec1 = f"@read{i} 1:N:0:ACGTACGT\n{r1}\n+\n{qual()}\n"
            rec2 = f"@read{i} 2:N:0:ACGTACGT\n{r2}\n+\n{qual()}\n"
            f1.write(rec1); f2.write(rec2); fi.write(rec1 + rec2)


if __name__ == "__main__":
    main()

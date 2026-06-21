"""FASTQ test data generation for benchmarking."""
import gzip
import random
import time
from pathlib import Path


def generate_data(r1_gz: Path, r2_gz: Path, num_pairs: int, seed: int = 2026):
    """Generate random paired-end FASTQ reads (independent R1/R2)."""
    READ_LEN = 150
    BATCH = 100_000
    BASES = b"ACGT"
    QUAL_POOL = 512
    rng = random.Random(seed)

    pool = []
    for _ in range(QUAL_POOL):
        q = bytearray(READ_LEN)
        for j in range(READ_LEN):
            if j < 5 or j > READ_LEN - 10:
                q[j] = rng.randint(20, 35) + 33
            else:
                q[j] = rng.randint(30, 40) + 33
        pool.append(bytes(q))

    t0 = time.time()
    written = 0
    with gzip.open(r1_gz, "wb", compresslevel=1) as f1, \
         gzip.open(r2_gz, "wb", compresslevel=1) as f2:
        while written < num_pairs:
            n = min(BATCH, num_pairs - written)
            b1 = bytearray()
            b2 = bytearray()
            for i in range(n):
                rid = written + i + 1
                name = f"@SIM:BENCH:1:{1101 + rid // 10000000}:{rid % 50000}:{(rid * 7) % 50000}".encode()
                s1 = bytes(rng.choices(BASES, k=READ_LEN))
                s2 = bytes(rng.choices(BASES, k=READ_LEN))
                q1 = pool[rng.randint(0, QUAL_POOL - 1)]
                q2 = pool[rng.randint(0, QUAL_POOL - 1)]
                b1 += name + b" 1:N:0:ATCG\n" + s1 + b"\n+\n" + q1 + b"\n"
                b2 += name + b" 2:N:0:ATCG\n" + s2 + b"\n+\n" + q2 + b"\n"
            f1.write(bytes(b1))
            f2.write(bytes(b2))
            written += n
            if written % 1_000_000 == 0:
                e = time.time() - t0
                eta = (num_pairs - written) / (written / e)
                print(f"  {100 * written / num_pairs:5.1f}%  {written / 1e6:.0f}M pairs"
                      f"  {written / e / 1e6:.2f}M/s  ETA {eta:.0f}s", flush=True)

    print(f"  Generated in {time.time() - t0:.1f}s")


def generate_merge_data(r1_gz: Path, r2_gz: Path, num_pairs: int, seed: int = 2026):
    """Generate overlapping PE reads for merge/correction testing.

    Each pair derives from a random fragment with insert size ~220bp (sd=30).
    R2 is the reverse complement of the fragment's 3' end.  Mismatches are
    injected in the overlap region with low quality on R2 to trigger base
    correction (R1 Q30+ vs R2 Q10-14 at mismatch positions).
    """
    READ_LEN = 150
    BATCH = 100_000
    BASES = b"ACGT"
    COMP = bytes.maketrans(b"ACGT", b"TGCA")
    QUAL_POOL = 512
    INSERT_MEAN = 220
    INSERT_SD = 30
    MIN_INSERT = READ_LEN + 30   # need >= 30bp overlap for fastp
    MAX_INSERT = 2 * READ_LEN - 1
    MISMATCH_RATE = 0.015

    rng = random.Random(seed)

    pool = []
    for _ in range(QUAL_POOL):
        q = bytearray(READ_LEN)
        for j in range(READ_LEN):
            if j < 5 or j > READ_LEN - 10:
                q[j] = rng.randint(25, 35) + 33
            else:
                q[j] = rng.randint(30, 40) + 33
        pool.append(bytes(q))

    def revcomp(seq: bytes) -> bytes:
        return seq.translate(COMP)[::-1]

    t0 = time.time()
    written = 0
    with gzip.open(r1_gz, "wb", compresslevel=1) as f1, \
         gzip.open(r2_gz, "wb", compresslevel=1) as f2:
        while written < num_pairs:
            n = min(BATCH, num_pairs - written)
            b1 = bytearray()
            b2 = bytearray()
            for i in range(n):
                rid = written + i + 1
                name = f"@SIM:MERGE:1:{1101 + rid // 10000000}:{rid % 50000}:{(rid * 7) % 50000}".encode()

                # Sample insert size from truncated normal distribution
                insert = int(rng.gauss(INSERT_MEAN, INSERT_SD))
                insert = max(MIN_INSERT, min(MAX_INSERT, insert))

                # Generate fragment and derive reads
                frag = bytearray(rng.choices(BASES, k=insert))
                s1 = bytes(frag[:READ_LEN])
                s2 = bytearray(revcomp(bytes(frag[insert - READ_LEN:])))

                q1 = bytearray(pool[rng.randint(0, QUAL_POOL - 1)])
                q2 = bytearray(pool[rng.randint(0, QUAL_POOL - 1)])

                # Inject mismatches in R2 overlap with low quality to trigger correction
                overlap_len = 2 * READ_LEN - insert
                n_mm = max(1, int(overlap_len * MISMATCH_RATE))
                n_mm = min(n_mm, min(4, int(overlap_len * 0.19)))
                overlap_start = READ_LEN - overlap_len
                for p in rng.sample(range(overlap_start, READ_LEN), n_mm):
                    alts = [b for b in BASES if b != s2[p]]
                    s2[p] = rng.choice(alts)
                    q2[p] = rng.randint(10 + 33, 14 + 33)  # Q10-14

                b1 += name + b" 1:N:0:ATCG\n" + s1 + b"\n+\n" + bytes(q1) + b"\n"
                b2 += name + b" 2:N:0:ATCG\n" + bytes(s2) + b"\n+\n" + bytes(q2) + b"\n"
            f1.write(bytes(b1))
            f2.write(bytes(b2))
            written += n
            if written % 1_000_000 == 0:
                e = time.time() - t0
                eta = (num_pairs - written) / (written / e)
                print(f"  {100 * written / num_pairs:5.1f}%  {written / 1e6:.0f}M pairs"
                      f"  {written / e / 1e6:.2f}M/s  ETA {eta:.0f}s", flush=True)

    print(f"  Generated in {time.time() - t0:.1f}s")

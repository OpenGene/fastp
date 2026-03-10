#!/usr/bin/env bash
set -euo pipefail

# End-to-end benchmark: optimized vs unoptimized fastp
# Tests both compressed (gz→gz) and uncompressed (fq→fq) modes.
# Usage: bash scripts/bench_e2e.sh [num_pairs] [threads]

NUM_PAIRS=${1:-10000000}   # default 10M pairs (~1x human)
THREADS=${2:-4}
RUNS=3                     # repeat each benchmark N times

BENCH_DIR="/tmp/fastp_benchmark"
DATA_DIR="$BENCH_DIR/data"
OUT_DIR="$BENCH_DIR/output"
FASTP_OPT="/tmp/fastp_opt"
FASTP_ORIG="/tmp/fastp_orig"

mkdir -p "$DATA_DIR" "$OUT_DIR"

echo "============================================"
echo "  fastp End-to-End Benchmark"
echo "============================================"
echo "  Pairs:   ${NUM_PAIRS}"
echo "  Threads: $THREADS"
echo "  Runs:    $RUNS (median reported)"
echo "============================================"
echo

# --- Step 1: Generate test data ---
R1_GZ="$DATA_DIR/bench_R1.fq.gz"
R2_GZ="$DATA_DIR/bench_R2.fq.gz"
R1_FQ="$DATA_DIR/bench_R1.fq"
R2_FQ="$DATA_DIR/bench_R2.fq"

if [[ -f "$R1_GZ" && -f "$R2_GZ" ]]; then
    echo "[data] Reusing existing compressed test data"
    echo "  R1.gz: $(du -h "$R1_GZ" | cut -f1)"
    echo "  R2.gz: $(du -h "$R2_GZ" | cut -f1)"
else
    echo "[data] Generating $NUM_PAIRS pairs..."
    python3 -c "
import gzip, random, os, sys, time

NUM_PAIRS = $NUM_PAIRS
READ_LEN = 150
BATCH = 100_000
BASES = b'ACGT'
QUAL_POOL = 512
rng = random.Random(42)

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
with gzip.open('$R1_GZ', 'wb', compresslevel=1) as f1, \
     gzip.open('$R2_GZ', 'wb', compresslevel=1) as f2:
    while written < NUM_PAIRS:
        n = min(BATCH, NUM_PAIRS - written)
        b1 = bytearray()
        b2 = bytearray()
        for i in range(n):
            rid = written + i + 1
            name = f'@SIM:BENCH:1:{1101 + rid//10000000}:{rid%50000}:{(rid*7)%50000}'.encode()
            s1 = bytes(rng.choices(BASES, k=READ_LEN))
            s2 = bytes(rng.choices(BASES, k=READ_LEN))
            q1 = pool[rng.randint(0, QUAL_POOL-1)]
            q2 = pool[rng.randint(0, QUAL_POOL-1)]
            b1 += name + b' 1:N:0:ATCG\n' + s1 + b'\n+\n' + q1 + b'\n'
            b2 += name + b' 2:N:0:ATCG\n' + s2 + b'\n+\n' + q2 + b'\n'
        f1.write(bytes(b1))
        f2.write(bytes(b2))
        written += n
        if written % 1000000 == 0:
            e = time.time() - t0
            print(f'  {100*written/NUM_PAIRS:5.1f}%  {written/1e6:.0f}M pairs  {written/e/1e6:.2f}M/s  ETA {(NUM_PAIRS-written)/(written/e):.0f}s', flush=True)

print(f'  Generated in {time.time()-t0:.1f}s')
"
    echo "  R1.gz: $(du -h "$R1_GZ" | cut -f1)"
    echo "  R2.gz: $(du -h "$R2_GZ" | cut -f1)"
fi

# Prepare uncompressed data
if [[ -f "$R1_FQ" && -f "$R2_FQ" ]]; then
    echo "[data] Reusing existing uncompressed test data"
else
    echo "[data] Decompressing to plain FASTQ..."
    gunzip -c "$R1_GZ" > "$R1_FQ"
    gunzip -c "$R2_GZ" > "$R2_FQ"
fi
echo "  R1.fq: $(du -h "$R1_FQ" | cut -f1)"
echo "  R2.fq: $(du -h "$R2_FQ" | cut -f1)"
echo

# --- Step 2: Warmup filesystem cache ---
echo "[cache] Warming up filesystem cache..."
cat "$R1_GZ" "$R2_GZ" "$R1_FQ" "$R2_FQ" > /dev/null
echo

# --- Step 3: Benchmark function ---
# Sets BENCH_MEDIAN after each call.
BENCH_MEDIAN=""

run_bench() {
    local label=$1
    local binary=$2
    local in1=$3
    local in2=$4
    local out_ext=$5   # "fq.gz" or "fq"
    local t_file="$BENCH_DIR/times_${label}.txt"
    rm -f "$t_file"

    echo "[bench] $label ($RUNS runs)"

    for run in $(seq 1 $RUNS); do
        rm -f "$OUT_DIR"/${label}_* "$OUT_DIR"/${label}.*

        local t_start=$(python3 -c "import time; print(time.time())")

        "$binary" \
            -i "$in1" -I "$in2" \
            -o "$OUT_DIR/${label}_R1.${out_ext}" \
            -O "$OUT_DIR/${label}_R2.${out_ext}" \
            -j "$OUT_DIR/${label}.json" \
            -h "$OUT_DIR/${label}.html" \
            -w "$THREADS" \
            2>/dev/null

        local t_end=$(python3 -c "import time; print(time.time())")
        local wall_s=$(python3 -c "print(f'{$t_end - $t_start:.2f}')")

        echo "$wall_s" >> "$t_file"
        echo "  Run $run: ${wall_s}s"
    done

    # Sort and pick median
    local median=$(sort -g "$t_file" | head -$(( (RUNS + 1) / 2 )) | tail -1)
    echo "  >> Median: ${median}s"
    echo
    BENCH_MEDIAN="$median"
}

# --- Step 4: Run benchmarks ---

echo "=========================================="
echo "  Mode 1: Compressed (gz -> gz)"
echo "=========================================="
echo
run_bench "orig_gz" "$FASTP_ORIG" "$R1_GZ" "$R2_GZ" "fq.gz"
MEDIAN_ORIG_GZ="$BENCH_MEDIAN"
run_bench "opt_gz" "$FASTP_OPT" "$R1_GZ" "$R2_GZ" "fq.gz"
MEDIAN_OPT_GZ="$BENCH_MEDIAN"

echo "=========================================="
echo "  Mode 2: Uncompressed (fq -> fq)"
echo "=========================================="
echo
run_bench "orig_fq" "$FASTP_ORIG" "$R1_FQ" "$R2_FQ" "fq"
MEDIAN_ORIG_FQ="$BENCH_MEDIAN"
run_bench "opt_fq" "$FASTP_OPT" "$R1_FQ" "$R2_FQ" "fq"
MEDIAN_OPT_FQ="$BENCH_MEDIAN"

# --- Step 5: Verify output correctness ---
echo "[verify] Checking output consistency..."

# md5 helper: works on both macOS (md5) and Linux (md5sum)
md5hash() { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }

# Compressed mode
ORIG_GZ_R1=$(gunzip -c "$OUT_DIR/orig_gz_R1.fq.gz" | md5)
OPT_GZ_R1=$(gunzip -c "$OUT_DIR/opt_gz_R1.fq.gz" | md5)
if [[ "$ORIG_GZ_R1" == "$OPT_GZ_R1" ]]; then
    echo "  gz R1: IDENTICAL"
else
    echo "  gz R1: DIFFERENT (WARNING)"
fi

ORIG_GZ_R2=$(gunzip -c "$OUT_DIR/orig_gz_R2.fq.gz" | md5)
OPT_GZ_R2=$(gunzip -c "$OUT_DIR/opt_gz_R2.fq.gz" | md5)
if [[ "$ORIG_GZ_R2" == "$OPT_GZ_R2" ]]; then
    echo "  gz R2: IDENTICAL"
else
    echo "  gz R2: DIFFERENT (WARNING)"
fi

# Uncompressed mode
ORIG_FQ_R1=$(md5hash "$OUT_DIR/orig_fq_R1.fq")
OPT_FQ_R1=$(md5hash "$OUT_DIR/opt_fq_R1.fq")
if [[ "$ORIG_FQ_R1" == "$OPT_FQ_R1" ]]; then
    echo "  fq R1: IDENTICAL"
else
    echo "  fq R1: DIFFERENT (WARNING)"
fi

ORIG_FQ_R2=$(md5hash "$OUT_DIR/orig_fq_R2.fq")
OPT_FQ_R2=$(md5hash "$OUT_DIR/opt_fq_R2.fq")
if [[ "$ORIG_FQ_R2" == "$OPT_FQ_R2" ]]; then
    echo "  fq R2: IDENTICAL"
else
    echo "  fq R2: DIFFERENT (WARNING)"
fi
echo

# --- Step 6: Summary ---
echo "============================================"
echo "  RESULTS"
echo "============================================"

python3 -c "
num = $NUM_PAIRS

gz_orig = float('${MEDIAN_ORIG_GZ}')
gz_opt  = float('${MEDIAN_OPT_GZ}')
fq_orig = float('${MEDIAN_ORIG_FQ}')
fq_opt  = float('${MEDIAN_OPT_FQ}')

def report(label, orig, opt):
    speedup = orig / opt
    pct = (1 - opt / orig) * 100
    print(f'  --- {label} ---')
    print(f'  Original:   {orig:.2f}s  ({num/orig/1e6:.2f}M pairs/s)')
    print(f'  Optimized:  {opt:.2f}s  ({num/opt/1e6:.2f}M pairs/s)')
    print(f'  Speedup:    {speedup:.2f}x ({pct:+.1f}%)')
    print(f'  Time saved: {orig - opt:.2f}s per run')
    print()

report('Compressed (gz -> gz)', gz_orig, gz_opt)
report('Uncompressed (fq -> fq)', fq_orig, fq_opt)
"
echo "============================================"

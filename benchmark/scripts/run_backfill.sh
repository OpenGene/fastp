#!/bin/bash
# Sweeps a set of fastp binaries x thread counts x corpus datasets, records
# wall-clock time and an output-record digest (for correctness comparison)
# to a CSV. Intended to run on a dedicated multi-core Linux box; a run that
# doesn't finish within --timeout is recorded as HUNG rather than blocking
# the sweep forever (relevant since some historical releases deadlock above
# 32 threads -- see https://github.com/OpenGene/fastp/pull/723).
#
# Usage: run_backfill.sh <bin_dir> <data_dir> <out_csv> [threads_csv] [timeout_s]
#   bin_dir/  must contain one executable per row of versions.tsv (see below)
#   data_dir/ must contain <dataset>_R1.fastq.gz / _R2.fastq.gz per benchmark/corpus.md
set -uo pipefail
BIN=${1:?usage: run_backfill.sh <bin_dir> <data_dir> <out_csv> [threads] [timeout_s]}
DATA=${2:?}
OUT=${3:?}
THREADS=${4:-"1,16,32,48"}
TIMEOUT=${5:-75}
REPS=${6:-1}
HERE=$(cd "$(dirname "$0")" && pwd)
VERSIONS="$HERE/versions.tsv"
DATASETS="atac wgs synth"

echo -e "version\tdataset\tthreads\trep\tresult\tsecs\tdigest" > "$OUT"

run_one() {
  local ver=$1 bin=$2 ds=$3 t=$4 rep=$5
  local r1="$DATA/${ds}_R1.fastq.gz" r2="$DATA/${ds}_R2.fastq.gz"
  [ -f "$r1" ] || { echo "skip $ds: missing $r1" >&2; return; }
  local work; work=$(mktemp -d)
  local start; start=$(date +%s.%N)
  timeout -s KILL "$TIMEOUT" "$BIN/$bin" -w "$t" --detect_adapter_for_pe \
    -i "$r1" -I "$r2" -o "$work/o1.fq.gz" -O "$work/o2.fq.gz" \
    -j "$work/r.json" -h "$work/r.html" > "$work/out.log" 2>"$work/err.log"
  local ec=$? secs; secs=$(python3 -c "print(round($(date +%s.%N)-$start,2))")
  local result=OK digest=-
  if [ $ec -eq 137 ]; then result=HUNG
  elif [ $ec -ne 0 ]; then result="EXIT$ec"
  else
    digest=$(for f in "$work/o1.fq.gz" "$work/o2.fq.gz"; do zcat "$f" | paste - - - - | sort | md5sum | cut -c1-12; done | tr '\n' '_')
  fi
  echo -e "$ver\t$ds\t$t\t$rep\t$result\t$secs\t$digest" | tee -a "$OUT"
  rm -rf "$work"
}

IFS=',' read -ra TLIST <<< "$THREADS"
while IFS=$'\t' read -r ver bin; do
  [ -z "$ver" ] && continue
  [[ "$ver" == \#* ]] && continue
  if [ ! -x "$BIN/$bin" ]; then echo "skip $ver: $BIN/$bin not found/executable" >&2; continue; fi
  for ds in $DATASETS; do
    for t in "${TLIST[@]}"; do
      for rep in $(seq 1 $REPS); do run_one "$ver" "$bin" "$ds" "$t" "$rep"; done
    done
  done
done < "$VERSIONS"

echo "BACKFILL DONE"

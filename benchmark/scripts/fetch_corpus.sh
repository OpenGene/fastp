#!/bin/bash
# Downloads and truncates the public corpus entries described in benchmark/corpus.md.
# Usage: fetch_corpus.sh <out_dir> [read_pairs_per_dataset]
set -euo pipefail
OUT=${1:?usage: fetch_corpus.sh <out_dir> [read_pairs]}
PAIRS=${2:-4000000}
LINES=$((PAIRS * 4))
mkdir -p "$OUT"

fetch() {
  local acc=$1 name=$2
  for mate in 1 2; do
    local dest="$OUT/${name}_R${mate}.fastq.gz"
    [ -s "$dest" ] && { echo "skip ${name}_R${mate}: already present"; continue; }
    local url="https://ftp.sra.ebi.ac.uk/vol1/fastq/${acc:0:6}/${acc}/${acc}_${mate}.fastq.gz"
    echo "fetching ${name}_R${mate} <- $url ($PAIRS pairs)"
    # `head -n` closes its read end once satisfied, sending SIGPIPE upstream;
    # with pipefail that looks like a failure even though the output is
    # complete and correct, so this pipeline's own exit status is ignored --
    # the file-size check right after is the real success signal.
    (curl -sf "$url" | gzip -dc | head -n "$LINES" | gzip -1 > "$dest") || true
    [ -s "$dest" ] || { echo "FAILED: $dest is empty" >&2; exit 1; }
  done
}

fetch SRR891268 atac
fetch SRR952827 wgs
ls -la "$OUT"

#!/usr/bin/env bash
set -euo pipefail

# Repro for issue #426: with -6/--phred64, every read must be converted to
# phred33, including reads whose buffers are recycled from the read pool
# (which only happens once enough reads have been processed).

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp/in.fq" <<'PY'
import sys
seq = 'ACGT' * 25
qual = 'h' * len(seq)   # Q40 in phred64, 'I' in phred33
with open(sys.argv[1], 'w') as f:
    for i in range(200000):
        f.write('@r%d\n%s\n+\n%s\n' % (i, seq, qual))
PY

./fastp -6 -w 4 \
  --disable_adapter_trimming --disable_trim_poly_g \
  --disable_quality_filtering --disable_length_filtering \
  -i "$tmp/in.fq" -o "$tmp/out.fq" \
  -j "$tmp/out.json" -h "$tmp/out.html" 2> "$tmp/log"

total=$(awk 'NR % 4 == 0' "$tmp/out.fq" | wc -l)
bad=$(awk 'NR % 4 == 0 && $0 !~ /^I+$/' "$tmp/out.fq" | wc -l)
if [[ "$total" -ne 200000 || "$bad" -ne 0 ]]; then
  echo "FAIL: $bad of $total output reads were not converted from phred64 to phred33"
  exit 1
fi

echo "PASS: issue #426 repro converted all $total reads to phred33"

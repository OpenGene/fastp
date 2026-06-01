#!/usr/bin/env bash
set -euo pipefail

# Repro for issue #697: --stdout should emit merged reads in merge mode.

python - <<'PY' > /tmp/fp_repro_697.interleaved.fq
seq='AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGC'
comp=str.maketrans('ACGTN','TGCAN')
rc=seq.translate(comp)[::-1]
qual='I'*len(seq)
print('@r1/1'); print(seq); print('+'); print(qual)
print('@r1/2'); print(rc);  print('+'); print(qual)
PY

./fastp \
  --disable_adapter_trimming --disable_trim_poly_g \
  --disable_quality_filtering --disable_length_filtering \
  --stdin --interleaved_in --merge --stdout \
  < /tmp/fp_repro_697.interleaved.fq \
  > /tmp/fp_697_stdout.fq 2> /tmp/fp_697_stdout.log

lines=$(wc -l < /tmp/fp_697_stdout.fq)
if [[ "$lines" -ne 4 ]]; then
  echo "FAIL: expected 4 FASTQ lines from --stdout merge mode, got $lines"
  exit 1
fi

echo "PASS: issue #697 repro produced $lines lines as expected"

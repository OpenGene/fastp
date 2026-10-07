#!/usr/bin/env bash
set -euo pipefail

# Repro for issue #638 (also #528): with --dedup, "reads passed filter" must count
# only the reads that are actually written, the removed duplicates must be reported
# as "duplicated_reads", and --merge mode must drop duplicates too.

python - <<'PY' > /tmp/fp_repro_638.fq
import random
random.seed(638)
reads=[''.join(random.choice('ACGT') for _ in range(100)) for _ in range(500)]
for tag in ('', '_dup'):                      # every read twice -> 1000 reads, 500 distinct
    for i, s in enumerate(reads):
        print(f'@r{i}{tag}'); print(s); print('+'); print('I' * 100)
PY

python - <<'PY' > /tmp/fp_repro_638.interleaved.fq
import random
comp = str.maketrans('ACGTN', 'TGCAN')
for tag in ('a', 'b'):                        # every overlapping pair twice -> 400 pairs, 200 distinct
    for i in range(200):
        random.seed(i); ins = ''.join(random.choice('ACGT') for _ in range(120))
        r1 = ins[:100]; r2 = ins.translate(comp)[::-1][:100]
        print(f'@p{i}{tag}/1'); print(r1); print('+'); print('I' * 100)
        print(f'@p{i}{tag}/2'); print(r2); print('+'); print('I' * 100)
PY

COMMON="--disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --dedup"

./fastp $COMMON -i /tmp/fp_repro_638.fq -o /tmp/fp_638_se.fq \
  -j /tmp/fp_638_se.json -h /tmp/fp_638_se.html 2> /tmp/fp_638_se.log
written=$(( $(wc -l < /tmp/fp_638_se.fq) / 4 ))
read -r passed duplicated after < <(python - <<'PY'
import json
j = json.load(open('/tmp/fp_638_se.json')); fr = j['filtering_result']
print(fr['passed_filter_reads'], fr.get('duplicated_reads', -1), j['summary']['after_filtering']['total_reads'])
PY
)
if [[ "$passed" -ne "$written" || "$after" -ne "$written" || "$duplicated" -ne $(( 1000 - written )) ]]; then
  echo "FAIL: SE --dedup wrote $written reads but reported passed_filter_reads=$passed duplicated_reads=$duplicated after_filtering.total_reads=$after"
  exit 1
fi
if ! grep -q "reads failed due to duplication: $duplicated" /tmp/fp_638_se.log; then
  echo "FAIL: SE --dedup did not report 'reads failed due to duplication: $duplicated' on stderr"
  exit 1
fi

./fastp $COMMON --stdin --interleaved_in --merge --merged_out /tmp/fp_638_merged.fq \
  -j /tmp/fp_638_pe.json -h /tmp/fp_638_pe.html < /tmp/fp_repro_638.interleaved.fq 2> /tmp/fp_638_pe.log
merged=$(( $(wc -l < /tmp/fp_638_merged.fq) / 4 ))
read -r passed duplicated < <(python - <<'PY'
import json
fr = json.load(open('/tmp/fp_638_pe.json'))['filtering_result']
print(fr['passed_filter_reads'], fr.get('duplicated_reads', -1))
PY
)
if [[ "$merged" -gt 200 ]]; then
  echo "FAIL: --merge --dedup wrote $merged merged reads from 200 distinct pairs"
  exit 1
fi
if [[ "$passed" -ne $(( 2 * merged )) || "$duplicated" -ne $(( 800 - 2 * merged )) ]]; then
  echo "FAIL: --merge --dedup wrote $merged merged reads but reported passed_filter_reads=$passed duplicated_reads=$duplicated"
  exit 1
fi

echo "PASS: issue #638 repro: SE wrote $written reads, merge wrote $merged reads, counts agree with the reports"

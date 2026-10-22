#!/usr/bin/env bash
set -euo pipefail

# Repro for issue #349: the "command" field of the JSON report must stay valid
# JSON when arguments contain quotes, backslashes or other special characters.

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

title='566136_S11 "fastp" report \ with	tab'

./fastp -i testdata/R1.fq -o "$tmp/out.fq" \
  -j "$tmp/out.json" -h "$tmp/out.html" -R "$title" 2> "$tmp/log"

if ! python3 - "$tmp/out.json" "$title" <<'PY'
import json, sys
with open(sys.argv[1]) as f:
    report = json.load(f)
assert sys.argv[2] in report["command"], report["command"]
PY
then
  echo "FAIL: JSON report is not valid or its command field is wrong"
  exit 1
fi

echo "PASS: issue #349 repro produced a valid JSON report"

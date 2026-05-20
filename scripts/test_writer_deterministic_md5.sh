#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

FASTP="${FASTP:-./fastp}"
if [ ! -x "$FASTP" ]; then
    echo "ERROR: fastp binary not found at '$FASTP' (run 'make' first)" >&2
    exit 2
fi

R1=testdata/R1.fq
R2=testdata/R2.fq
for f in "$R1" "$R2"; do
    [ -r "$f" ] || { echo "ERROR: missing fixture $f" >&2; exit 2; }
done

TMP="$(mktemp -d 2>/dev/null || mktemp -d -t fastp-md5)"
trap 'rm -rf "$TMP"' EXIT

run_fastp() {
    local run_id="$1"
    perl -e 'alarm shift @ARGV; exec @ARGV or die "exec: $!"' 30 \
        "$FASTP" \
        -i "$R1" -I "$R2" \
        -o "$TMP/o1.$run_id.fq" -O "$TMP/o2.$run_id.fq" \
        -w 16 \
        --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
        --json "$TMP/r.$run_id.json" --html "$TMP/r.$run_id.html" \
        >"$TMP/run.$run_id.log" 2>&1
}

md5_file() {
    if command -v md5sum >/dev/null 2>&1; then
        set -- $(md5sum "$1")
        printf '%s\n' "$1"
    else
        md5 -q "$1"
    fi
}

printf '[case] %-50s ' "PE w=16 repeat decompressed MD5"
run_fastp 1
base1="$(md5_file "$TMP/o1.1.fq")"
base2="$(md5_file "$TMP/o2.1.fq")"
for i in 2 3 4 5; do
    run_fastp "$i"
    md51="$(md5_file "$TMP/o1.$i.fq")"
    md52="$(md5_file "$TMP/o2.$i.fq")"
    if [ "$md51" != "$base1" ] || [ "$md52" != "$base2" ]; then
        echo "FAIL"
        echo "R1 baseline=$base1 run$i=$md51" >&2
        echo "R2 baseline=$base2 run$i=$md52" >&2
        exit 1
    fi
done
echo "OK"

echo
echo "ALL PASSED"

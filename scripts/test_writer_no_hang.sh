#!/usr/bin/env bash
# Regression test for WriterThread round-robin SPSC deadlock.
# Background: src/writerthread.cpp:output() rotates per-tid SPSC lists. When a
# slot is permanently empty (producerFinished && isEmpty) the rotation must
# skip the slot; otherwise the writer loops on usleep forever while other slots
# still hold items. See: project_master_writer_skip_drained.md.
#
# Each case must finish within the per-case timeout. SIGALRM (perl alarm) kills
# the binary on hang, producing a non-zero exit that this script reports.

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

TMP="$(mktemp -d 2>/dev/null || mktemp -d -t fastp-spsc)"
trap 'rm -rf "$TMP"' EXIT

run_case() {
    local label="$1"; shift
    local timeout_s="$1"; shift
    printf '[case] %-50s ' "$label"
    if perl -e 'alarm shift @ARGV; exec @ARGV or die "exec: $!"' "$timeout_s" \
            "$FASTP" "$@" >"$TMP/out.log" 2>&1; then
        echo "OK"
    else
        local rc=$?
        echo "FAIL (exit $rc, possibly SIGALRM=hang)"
        echo "----- tail of $TMP/out.log -----" >&2
        tail -30 "$TMP/out.log" >&2 || true
        exit 1
    fi
}

# Trailing-empty slot: 8 reads -> 1 pack, w=16 workers => 15 lists never produce.
run_case "PE w=16 trailing-empty"  15 \
  -i "$R1" -I "$R2" -o "$TMP/o1.fq" -O "$TMP/o2.fq" -w 16 \
  --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
  --json "$TMP/r.json" --html "$TMP/r.html"

# Optional writers (--failed_out, --unpaired*) take the conditional-input path
# that historically had the highest deadlock blast radius.
run_case "PE w=16 + failed_out + unpaired" 15 \
  -i "$R1" -I "$R2" -o "$TMP/o1.fq" -O "$TMP/o2.fq" -w 16 \
  --failed_out "$TMP/failed.fq" --unpaired1 "$TMP/u1.fq" --unpaired2 "$TMP/u2.fq" \
  --json "$TMP/r.json" --html "$TMP/r.html"

run_case "SE w=16 trailing-empty"  15 \
  -i "$R1" -o "$TMP/o.fq" -w 16 \
  --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering \
  --json "$TMP/r.json" --html "$TMP/r.html"

# Race-prone repeat to catch timing-dependent hangs.
for i in 1 2 3 4 5; do
    run_case "PE w=64 trailing-empty repeat $i" 20 \
      -i "$R1" -I "$R2" -o "$TMP/o1.fq" -O "$TMP/o2.fq" -w 64 \
      --json "$TMP/r.json" --html "$TMP/r.html"
done

echo
echo "ALL PASSED"

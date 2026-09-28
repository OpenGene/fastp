#!/usr/bin/env python3
"""Correctness checks that don't need a golden output file.

1. Thread-count invariance: every output mode must produce byte-identical
   output (after decompression) at every thread count, and must finish.
   Thread counts above 32 would have caught #721 (reader/worker deadlock).
   fastp caps -w to the machine's core count, which on a CI runner would
   silently turn -w 48 into -w 4; on Linux the fake_nprocs.c LD_PRELOAD shim
   makes fastp see 64 CPUs so high worker counts are really exercised, and any
   run that still reports capping fails the check.
2. Adapter auto-detection finds the adapter planted by gen_reads.py.

usage: check_correctness.py FASTP [--threads 1,2,4,9,17,33,48] [--pairs 40000]
Writes a markdown summary to stdout and to $GITHUB_STEP_SUMMARY if set.
Exits non-zero on any mismatch, hang, crash, or missed adapter.
"""
import argparse, glob, gzip, hashlib, os, platform, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
TIMEOUT = 180

# name -> (args, files to digest, digest per file vs. over sorted records of all files)
MODES = {
    "pe_gz":           ("-i {r1} -I {r2} -o o1.fq.gz -O o2.fq.gz", ["o1.fq.gz", "o2.fq.gz"], False),
    "pe_plain":        ("-i {r1} -I {r2} -o o1.fq -O o2.fq", ["o1.fq", "o2.fq"], False),
    "se_gz":           ("-i {r1} -o o1.fq.gz", ["o1.fq.gz"], False),
    "interleaved_in":  ("--interleaved_in -i {il} -o o1.fq.gz -O o2.fq.gz", ["o1.fq.gz", "o2.fq.gz"], False),
    "merge":           ("-i {r1} -I {r2} -m --merged_out m.fq.gz -o o1.fq.gz -O o2.fq.gz", ["m.fq.gz", "o1.fq.gz", "o2.fq.gz"], False),
    "failed_unpaired": ("-i {r1} -I {r2} -o o1.fq.gz -O o2.fq.gz -q 30 -u 30 --failed_out failed.fq.gz "
                        "--unpaired1 u1.fq.gz --unpaired2 u2.fq.gz", ["o1.fq.gz", "o2.fq.gz", "failed.fq.gz", "u1.fq.gz", "u2.fq.gz"], False),
    "stdout":          ("-i {r1} -I {r2} --stdout", ["stdout"], False),
    # split output assigns files per worker, so compare record content, not file layout
    "split":           ("-i {r1} -I {r2} -o o1.fq -O o2.fq -s 3", ["*.o1.fq", "*.o2.fq"], True),
}


def read_bytes(path) -> bytes:
    if path.endswith(".gz"):
        with gzip.GzipFile(path, "rb") as f:
            return f.read()
    with open(path, "rb") as f:
        return f.read()


def digest(workdir, patterns, sort_records):
    out = []
    for pat in patterns:
        paths = sorted(glob.glob(os.path.join(workdir, pat)))
        if not paths:
            return "missing:" + pat
        if sort_records:
            recs = []
            for p in paths:
                lines = read_bytes(p).split(b"\n")
                recs += [b"\n".join(lines[i:i + 4]) for i in range(0, len(lines) - 3, 4)]
            data = b"\n".join(sorted(recs))
        else:
            data = b"".join(read_bytes(p) for p in paths)
        out.append(hashlib.md5(data).hexdigest()[:10])
    return "_".join(out)


ENV = dict(os.environ)


def prepare_high_thread_env(tmp):
    """On Linux, build the nprocs shim so fastp doesn't cap -w to the runner's cores.
    Returns the max worker count that will really run."""
    if platform.system() != "Linux":
        return os.cpu_count() or 1
    so = os.path.join(tmp, "fake_nprocs.so")
    subprocess.run(["cc", "-shared", "-fPIC", "-O2", "-o", so, os.path.join(HERE, "fake_nprocs.c"), "-ldl"], check=True)
    ENV["LD_PRELOAD"] = so
    ENV["FAKE_NPROCS"] = "64"
    return 64


def run(fastp, args, workdir):
    stdout = open(os.path.join(workdir, "stdout"), "wb")
    try:
        p = subprocess.run([fastp] + args.split() + ["-j", "r.json", "-h", "r.html"], cwd=workdir,
                           stdout=stdout, stderr=subprocess.PIPE, timeout=TIMEOUT, env=ENV)
        err = p.stderr.decode(errors="replace")
        if "Reduce worker threads" in err:
            return "CAPPED", err
        return ("OK" if p.returncode == 0 else f"EXIT{p.returncode}"), err
    except subprocess.TimeoutExpired:
        return "HUNG", ""
    finally:
        stdout.close()


def detected(stderr, mate):
    """Adapter lines printed between 'Detecting adapter sequence for readN...' and the next blank line."""
    lines = stderr.split("\n")
    marker = f"Detecting adapter sequence for read{mate}..."
    if marker not in lines:
        return None
    block = lines[lines.index(marker) + 1:]
    block = block[:block.index("")] if "" in block else block
    seqs = [l for l in block if l and set(l) <= set("ACGTN")]
    return seqs[-1] if seqs else ""


def compatible(found, planted):
    # detection may legitimately report a shorter/longer known adapter from the same family
    return bool(found) and len(found) >= 12 and (planted.startswith(found) or found.startswith(planted))


def main():
    global TIMEOUT
    ap = argparse.ArgumentParser()
    ap.add_argument("fastp")
    ap.add_argument("--threads", default="1,2,4,9,17,33,48")
    ap.add_argument("--pairs", type=int, default=40000)
    ap.add_argument("--timeout", type=int, default=TIMEOUT)
    a = ap.parse_args()
    TIMEOUT = a.timeout
    fastp = os.path.abspath(a.fastp)
    tmp = tempfile.mkdtemp(prefix="fastp-ci-")
    max_workers = prepare_high_thread_env(tmp)
    requested = [int(t) for t in a.threads.split(",")]
    threads = [t for t in requested if t <= max_workers or t <= 4]
    skipped = [t for t in requested if t not in threads]
    subprocess.run([sys.executable, os.path.join(HERE, "gen_reads.py"), os.path.join(tmp, "in"), "--pairs", str(a.pairs)], check=True)
    files = {"r1": os.path.join(tmp, "in_R1.fastq.gz"), "r2": os.path.join(tmp, "in_R2.fastq.gz"),
             "il": os.path.join(tmp, "in_interleaved.fastq.gz")}

    failures = []
    rows = ["| mode | " + " | ".join(f"-w {t}" for t in threads) + " |", "|---" * (len(threads) + 1) + "|"]
    for mode, (args, pats, sort_records) in MODES.items():
        ref, cells = None, []
        for t in threads:
            wd = tempfile.mkdtemp(dir=tmp)
            res, _ = run(fastp, args.format(**files) + f" -w {t}", wd)
            if res != "OK":
                cells.append(f"**{res}**"); failures.append(f"{mode} -w {t}: {res}"); continue
            d = digest(wd, pats, sort_records)
            ref = ref or d
            if d != ref:
                cells.append(f"**differs** `{d}`"); failures.append(f"{mode} -w {t}: output differs from -w {threads[0]}")
            else:
                cells.append(f"`{d}`" if t == threads[0] else "same")
        rows.append(f"| {mode} | " + " | ".join(cells) + " |")

    det_rows = ["| input | read1 adapter | read2 adapter | result |", "|---|---|---|---|"]
    for name, args, expect in [("SE (default detection)", "-i {r1} -o o.fq.gz", [ADAPTER_R1]),
                               ("PE --detect_adapter_for_pe", "-i {r1} -I {r2} -o o1.fq.gz -O o2.fq.gz --detect_adapter_for_pe",
                                [ADAPTER_R1, ADAPTER_R2])]:
        res, err = run(fastp, args.format(**files) + " -w 4", tempfile.mkdtemp(dir=tmp))
        found = [detected(err, m + 1) for m in range(len(expect))]
        ok = res == "OK" and all(compatible(f, e) for f, e in zip(found, expect))
        if not ok:
            failures.append(f"adapter detection ({name}): got {found}, expected {expect} [{res}]")
        det_rows.append(f"| {name} | `{found[0]}` | " + (f"`{found[1]}`" if len(found) > 1 else "n/a") + f" | {'OK' if ok else '**FAIL**'} |")

    note = ([f"Skipped -w {','.join(map(str, skipped))}: this runner can only run {max_workers} real workers "
             "(no nprocs shim on this OS).", ""] if skipped else [])
    report = ["## fastp correctness checks", "",
              f"Thread-count invariance on {a.pairs} synthetic read pairs (output digests, decompressed):", ""] + note + rows + \
             ["", "Adapter auto-detection on the same reads (TruSeq adapters planted):", ""] + det_rows + [""]
    report += ["**All checks passed.**"] if not failures else ["**Failures:**"] + [f"- {f}" for f in failures]
    text = "\n".join(report) + "\n"
    print(text)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            f.write(text)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()

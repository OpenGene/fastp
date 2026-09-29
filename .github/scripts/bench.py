#!/usr/bin/env python3
"""Benchmark fastp builds against each other, interleaved on the same machine.

Shared CI runners are too noisy to compare against numbers stored from another
run, so a PR is measured as base vs head in one job, alternating builds each
repetition. Each run is split into stages from fastp's stderr timestamps:
adapter auto-detection (serial, before processing) and processing; plus CPU
time and peak RSS.

Subsets make fixed-cost stages look bigger than they are: pre-processing,
adapter detection (at most 256K reads / 39.6M bases per mate) and the report
cost the same on a 300K-read subset as on a 50M-read run. So each run is also
projected to a full-size run of --project-reads reads (pairs for PE):
wall = fixed stages + processing x (reads / subset reads). Fixed stages run on
one thread, so projected CPU = fixed + (CPU - fixed) x the same factor.
Validated against 30 full-size public runs (6 datasets, 2 builds, -w 8/16/48):
median error 8% on wall time, 3 percentage points on base-vs-head deltas.

A gated metric (projected wall, projected CPU, peak RSS) regresses when its median is worse than
base by more than --threshold percent AND every head run is worse than every
base run, so a single noisy run can't trip it.

  bench.py prepare DATA_DIR            # synthetic sets + a cached public subset
  bench.py run --builds base=PATH,head=PATH --data DATA_DIR [--threshold 10] [--fail-on-regression] ...
"""
import argparse, gzip, json, shutil, os, platform, statistics, subprocess, sys, tempfile, threading, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
# input sets: name -> public ENA run (first N pairs streamed) or None for gen_reads.py output
SOURCES = {
    "synthetic": (None, 300000),  # 45M bases: past the detection cap, as on a full file
    "atac_hiseq": ("SRR891268", 500000),  # Buenrostro 2013 GM12878 ATAC-seq, Nextera adapters
}
# (benchmark name, layout, source)
DATASETS = [
    ("synthetic_pe", "PE", "synthetic"),
    ("synthetic_se", "SE", "synthetic"),
    ("atac_hiseq_pe", "PE", "atac_hiseq"),
]
METRICS = [("proj_wall", "s", "projected wall"), ("proj_cpu", "s", "projected CPU"), ("rss_mb", "MB", "peak RSS"),
           ("wall", "s", "wall time"), ("cpu", "s", "CPU (user+sys)"),
           ("detect", "s", "adapter detection"), ("process", "s", "processing")]
GATED = ["proj_wall", "proj_cpu", "rss_mb"]
MEASURED = ["wall", "cpu", "detect", "process"]
MIN_BASE = 0.2  # seconds; shorter timings are too noisy to judge
NOTABLE_PCT = 5.0  # improvements are reported (never gated) from here
SIGNOFF_LABEL = "perf-regression-ok"
MARKER = "<!-- fastp-benchmark -->"  # identifies the sticky PR comment


def prepare(data):
    os.makedirs(data, exist_ok=True)
    for src, (acc, n) in SOURCES.items():
        if os.path.exists(os.path.join(data, f"{src}_R2.fastq.gz")):
            continue
        if acc is None:
            subprocess.run([sys.executable, os.path.join(HERE, "gen_reads.py"), os.path.join(data, src),
                            "--pairs", str(n), "--seed", "7"], check=True)
            continue
        for mate in (1, 2):
            url = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{acc[:6]}/{acc}/{acc}_{mate}.fastq.gz"
            dest = os.path.join(data, f"{src}_R{mate}.fastq.gz")
            with urllib.request.urlopen(url, timeout=120) as resp, gzip.GzipFile(fileobj=resp) as inp, \
                 gzip.open(dest + ".tmp", "wb", compresslevel=1) as out:
                for _ in range(n * 4):
                    line = inp.readline()
                    if not line:
                        break
                    out.write(line)
            os.rename(dest + ".tmp", dest)


def run_once(fastp, args, timeout):
    """Returns (exit code, wall s, [(t, stderr line)], cpu s, peak RSS MB) for one run."""
    t0 = time.time()
    p = subprocess.Popen([fastp] + args, stderr=subprocess.PIPE, stdout=subprocess.DEVNULL, text=True)
    timer = threading.Timer(timeout, p.kill)
    timer.start()
    assert p.stderr is not None
    lines = [(time.time() - t0, l.rstrip("\n")) for l in p.stderr]  # fastp's stderr is unbuffered
    _, status, ru = os.wait4(p.pid, 0)  # this child's own rusage, not the cumulative RUSAGE_CHILDREN
    wall = time.time() - t0
    timer.cancel()
    p.returncode = os.waitstatus_to_exitcode(status)
    rss_mb = ru.ru_maxrss / (1 << 20 if platform.system() == "Darwin" else 1 << 10)
    return p.returncode, wall, lines, ru.ru_utime + ru.ru_stime, rss_mb


def stages(lines, wall):
    """(detect, process) seconds from stderr timestamps; detect is 0 when fastp skipped it."""
    stat = next((t for t, l in lines if l.startswith("Read1 before filtering")), wall)
    det = [i for i, (_, l) in enumerate(lines) if l.startswith("Detecting adapter")]
    if not det:
        return 0.0, stat
    end = next((lines[i][0] for i in range(det[-1] + 1, len(lines)) if lines[i][1] == ""), stat)
    return end - lines[det[0]][0], stat - end


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare"); p.add_argument("data")
    r = sub.add_parser("run")
    r.add_argument("--builds", required=True, help="name=PATH,...; the first is the baseline, the last is compared to it")
    r.add_argument("--data", required=True)
    r.add_argument("--threads", default="1,4")
    r.add_argument("--reps", type=int, default=3)
    r.add_argument("--timeout", type=int, default=600)
    r.add_argument("--threshold", type=float, default=10.0, help="regression threshold, percent")
    r.add_argument("--project-reads", type=float, default=50e6, help="full-size run to project to (reads, or pairs for PE)")
    r.add_argument("--fail-on-regression", action="store_true")
    r.add_argument("--signed-off", action="store_true", help=f"regressions accepted (the {SIGNOFF_LABEL} label)")
    r.add_argument("--summary", help="append the report here (e.g. $GITHUB_STEP_SUMMARY)")
    r.add_argument("--comment", help="write the report here, for posting as a PR comment")
    r.add_argument("--compare-json", help="base vs head per cell, with verdicts")
    r.add_argument("--json", help="head-only values for github-action-benchmark")
    r.add_argument("--tsv", help="every run")
    a = ap.parse_args()
    if a.cmd == "prepare":
        return prepare(a.data)

    builds = [tuple(b.split("=", 1)) for b in a.builds.split(",")]
    threads = [int(t) for t in a.threads.split(",")]
    rows, failures = [], []
    for name, layout, src in DATASETS:
        d = lambda m: os.path.join(a.data, f"{src}_R{m}.fastq.gz")
        for t in threads:
            for rep in range(a.reps):
                for bname, path in (builds if rep % 2 == 0 else builds[::-1]):  # alternate order: no first-run bias
                    w = tempfile.mkdtemp()
                    if layout == "PE":
                        args = ["-i", d(1), "-I", d(2), "-o", f"{w}/o1.fq.gz", "-O", f"{w}/o2.fq.gz", "--detect_adapter_for_pe"]
                    else:
                        args = ["-i", d(1), "-o", f"{w}/o1.fq.gz"]
                    args += ["-w", str(t), "-j", f"{w}/r.json", "-h", f"{w}/r.html"]
                    rc, wall, lines, cpu, rss = run_once(os.path.abspath(path), args, a.timeout)
                    reads = 0
                    if rc == 0:
                        with open(f"{w}/r.json") as f:
                            reads = json.load(f)["summary"]["before_filtering"]["total_reads"] // (2 if layout == "PE" else 1)
                    shutil.rmtree(w, ignore_errors=True)
                    if rc != 0 or not reads:
                        failures.append(f"{bname} {name} -w {t}: exit {rc}")
                        continue
                    det, proc = stages(lines, wall)
                    fixed, scale = wall - proc, a.project_reads / reads
                    rows.append(dict(build=bname, dataset=name, threads=t, rep=rep, wall=wall,
                                     detect=det, process=proc, cpu=cpu, rss_mb=rss,
                                     proj_wall=fixed + proc * scale,
                                     proj_cpu=fixed + max(cpu - fixed, 0) * scale))

    names = [b for b, _ in builds]
    base, head = names[0], names[-1]
    compare = len(names) > 1

    def vals(b, ds, t, k):
        return [r[k] for r in rows if r["build"] == b and r["dataset"] == ds and r["threads"] == t]

    def med(b, ds, t, k):
        v = vals(b, ds, t, k)
        return statistics.median(v) if v else float("nan")

    def verdict(ds, t, k):
        vb, vh = vals(base, ds, t, k), vals(head, ds, t, k)
        if not vb or not vh or (k != "rss_mb" and statistics.median(vb) < MIN_BASE):
            return None, "n/a"
        delta = 100 * (statistics.median(vh) - statistics.median(vb)) / statistics.median(vb)
        if delta > a.threshold and min(vh) > max(vb):
            return delta, "regression"
        if delta < -NOTABLE_PCT and max(vh) < min(vb):
            return delta, "improvement"
        return delta, "noise"

    cells = [(n, t, k) for n, _, _ in DATASETS for t in threads for k, _, _ in METRICS]
    results = {c: verdict(*c) for c in cells} if compare else {}
    regressions = [c for c in cells if c[2] in GATED and results.get(c, (0, ""))[1] == "regression"]
    improvements = [c for c in cells if c[2] in GATED and results.get(c, (0, ""))[1] == "improvement"]
    icon = {"regression": " 🔴", "improvement": " 🟢"}
    size = f"{a.project_reads / 1e6:g}M reads (pairs for PE)"

    def table(keys):
        out = ["| dataset | -w | " + " | ".join(label for k, _, label in METRICS if k in keys) + " |",
               "|---|---|" + "---|" * len(keys)]
        for n, _, _ in DATASETS:
            for t in threads:
                row = []
                for k, unit, _ in METRICS:
                    if k not in keys:
                        continue
                    if not compare:
                        row.append(f"{med(head, n, t, k):.2f} {unit}")
                        continue
                    d, v = results[(n, t, k)]
                    row.append(f"{med(base, n, t, k):.2f} → {med(head, n, t, k):.2f} {unit}"
                               + (f" ({d:+.1f}%){icon.get(v, '')}" if d is not None else ""))
                out.append(f"| {n} | {t} | " + " | ".join(row) + " |")
        return out

    out = [MARKER, "## fastp benchmark", ""]
    if compare:
        if failures:
            out.append("❌ **Some benchmark runs failed** (listed below).")
        elif regressions and a.signed_off:
            out.append(f"🟡 **{len(regressions)} regression(s) over {a.threshold:g}%, accepted** via the `{SIGNOFF_LABEL}` label.")
        elif regressions:
            out.append(f"🔴 **{len(regressions)} regression(s) over {a.threshold:g}%.** "
                       f"If intended, a maintainer can accept them by adding the `{SIGNOFF_LABEL}` label.")
        elif improvements:
            out.append(f"🟢 **No regressions; {len(improvements)} improvement(s).**")
        else:
            out.append(f"✅ **No regressions** over {a.threshold:g}%.")
        out += ["", f"`{base}` → `{head}`, median of {a.reps} interleaved runs on one {os.cpu_count()}-CPU runner. "
                f"🔴 = worse by more than {a.threshold:g}%, 🟢 = better by more than {NOTABLE_PCT:g}%, in both cases "
                "with no overlap between the base and head runs; unmarked changes are within runner noise.", "",
                f"**Projected to a full-size run of {size}** from the measured subset: fixed stages (pre-processing, "
                "adapter detection, report) as measured, processing scaled by read count. On 30 full-size public "
                "runs this was within 8% (median) of measured wall time and within 3 points on base-vs-head deltas. "
                "Projected wall, projected CPU and peak RSS are gated.", ""]
    out += table(GATED) + ["", "<details><summary>measured on the subsets</summary>", ""] + \
        table(MEASURED) + ["", "</details>", ""]
    if failures:
        out += ["**Failed runs:**"] + [f"- {f}" for f in failures]
    text = "\n".join(out) + "\n"
    print(text)
    if a.summary:
        with open(a.summary, "a") as f:
            f.write(text)
    if a.comment:
        with open(a.comment, "w") as f:
            f.write(text)
    if a.tsv:
        cols = ("build", "dataset", "threads", "rep", "wall", "detect", "process", "cpu", "rss_mb", "proj_wall", "proj_cpu")
        with open(a.tsv, "w") as f:
            f.write("\t".join(cols) + "\n")
            for row in rows:
                f.write("\t".join(str(row[k]) for k in cols) + "\n")
    if a.json:
        entries = [{"name": f"{n} -w{t} {label}", "unit": unit, "value": round(med(head, n, t, k), 3)}
                   for n, _, _ in DATASETS for t in threads for k, unit, label in METRICS]
        with open(a.json, "w") as f:
            json.dump(entries, f, indent=1)
    if a.compare_json and compare:
        with open(a.compare_json, "w") as f:
            json.dump({"base": base, "head": head, "threshold_pct": a.threshold, "signed_off": a.signed_off,
                       "project_reads": a.project_reads,
                       "failures": failures, "regressions": len(regressions), "improvements": len(improvements),
                       "cells": [{"dataset": n, "threads": t, "metric": k, "gated": k in GATED,
                                  "base": med(base, n, t, k), "head": med(head, n, t, k),
                                  "delta_pct": results[(n, t, k)][0], "verdict": results[(n, t, k)][1]}
                                 for n, t, k in cells]}, f, indent=1)
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as f:
            f.write(f"regressions={len(regressions)}\nimprovements={len(improvements)}\n")
    sys.exit(1 if failures or (regressions and a.fail_on_regression and not a.signed_off) else 0)


if __name__ == "__main__":
    main()

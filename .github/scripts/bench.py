#!/usr/bin/env python3
"""Benchmark fastp builds against each other, interleaved on the same machine.

Shared CI runners are too noisy to compare against numbers stored from another
run, so a PR is measured as base vs head in one job, alternating builds each
repetition. Each run is split into stages from fastp's stderr timestamps:
adapter auto-detection (serial, before processing) and processing; plus CPU
time and peak RSS.

  bench.py prepare DATA_DIR            # synthetic sets + a cached public subset
  bench.py run --builds base=PATH,head=PATH --data DATA_DIR [--threads 1,4] [--reps 3]
               [--summary FILE] [--json FILE] [--tsv FILE]
"""
import argparse, gzip, json, shutil, os, platform, statistics, subprocess, sys, tempfile, threading, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
# input sets: name -> public ENA run (first N pairs streamed) or None for gen_reads.py output
SOURCES = {
    "synthetic": (None, 200000),
    "atac_hiseq": ("SRR891268", 500000),  # Buenrostro 2013 GM12878 ATAC-seq, Nextera adapters
}
# (benchmark name, layout, source)
DATASETS = [
    ("synthetic_pe", "PE", "synthetic"),
    ("synthetic_se", "SE", "synthetic"),
    ("atac_hiseq_pe", "PE", "atac_hiseq"),
]


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
    r.add_argument("--builds", required=True)
    r.add_argument("--data", required=True)
    r.add_argument("--threads", default="1,4")
    r.add_argument("--reps", type=int, default=3)
    r.add_argument("--timeout", type=int, default=600)
    r.add_argument("--summary"); r.add_argument("--json"); r.add_argument("--tsv")
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
                for bname, path in builds:
                    w = tempfile.mkdtemp()
                    if layout == "PE":
                        args = ["-i", d(1), "-I", d(2), "-o", f"{w}/o1.fq.gz", "-O", f"{w}/o2.fq.gz", "--detect_adapter_for_pe"]
                    else:
                        args = ["-i", d(1), "-o", f"{w}/o1.fq.gz"]
                    args += ["-w", str(t), "-j", f"{w}/r.json", "-h", f"{w}/r.html"]
                    rc, wall, lines, cpu, rss = run_once(os.path.abspath(path), args, a.timeout)
                    if rc != 0:
                        failures.append(f"{bname} {name} -w {t}: exit {rc}")
                        continue
                    det, proc = stages(lines, wall)
                    rows.append(dict(build=bname, dataset=name, threads=t, rep=rep, wall=wall,
                                     detect=det, process=proc, cpu=cpu, rss_mb=rss))
                    shutil.rmtree(w, ignore_errors=True)

    def med(b, ds, t, k):
        v = [r[k] for r in rows if r["build"] == b and r["dataset"] == ds and r["threads"] == t]
        return statistics.median(v) if v else float("nan")

    names = [b for b, _ in builds]
    metrics = [("wall", "s", "wall time"), ("detect", "s", "adapter detection"), ("process", "s", "processing"),
               ("cpu", "s", "CPU time (user+sys)"), ("rss_mb", "MB", "peak RSS")]
    out = ["## fastp benchmark", "",
           f"Median of {a.reps} interleaved runs per cell on this runner ({os.cpu_count()} CPUs). "
           + (f"Delta is `{names[-1]}` vs `{names[0]}`; runner noise is typically a few percent, "
              "so treat small deltas as noise." if len(names) > 1 else ""), ""]
    for key, unit, label in metrics:
        out += [f"**{label}** ({unit})", "",
                "| dataset | -w | " + " | ".join(names) + (" | delta |" if len(names) > 1 else " |"),
                "|---|---|" + "---|" * len(names) + ("---|" if len(names) > 1 else "")]
        for name, _, _ in DATASETS:
            for t in threads:
                vals = [med(b, name, t, key) for b in names]
                cells = " | ".join(f"{v:.2f}" for v in vals)
                delta = ""
                if len(names) > 1:
                    base, head = vals[0], vals[-1]
                    delta = " | " + (f"{100 * (head - base) / base:+.1f}%" if base > 0.05 else "n/a")
                out.append(f"| {name} | {t} | {cells}{delta} |")
        out.append("")
    if failures:
        out += ["**Failed runs:**"] + [f"- {f}" for f in failures]
    text = "\n".join(out) + "\n"
    print(text)
    if a.summary:
        with open(a.summary, "a") as f:
            f.write(text)
    if a.tsv:
        with open(a.tsv, "w") as f:
            f.write("build\tdataset\tthreads\trep\twall\tdetect\tprocess\tcpu\trss_mb\n")
            for r in rows:
                f.write("\t".join(str(r[k]) for k in ("build", "dataset", "threads", "rep", "wall", "detect", "process", "cpu", "rss_mb")) + "\n")
    if a.json:
        head = names[-1]
        entries = [{"name": f"{name} -w{t} {label}", "unit": unit, "value": round(med(head, name, t, key), 3)}
                   for name, _, _ in DATASETS for t in threads for key, unit, label in metrics]
        with open(a.json, "w") as f:
            json.dump(entries, f, indent=1)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()

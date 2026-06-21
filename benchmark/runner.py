"""Benchmark execution engine."""
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from statistics import median

from .modes import parse_mode, FEAT_DEFS, OUT_DIR


def output_files(label: str, mode: str) -> list[Path]:
    """Return list of output file paths for a given label+mode."""
    mtype, _, out_fmt = parse_mode(mode)
    if mtype == "stdin":
        return []
    ext = "fq" if out_fmt == "fq" else "fq.gz"
    if mode == "merge":
        return [OUT_DIR / f"{label}_merged.{ext}",
                OUT_DIR / f"{label}_R1.{ext}",
                OUT_DIR / f"{label}_R2.{ext}"]
    if mtype == "se":
        return [OUT_DIR / f"{label}_R1.{ext}"]
    return [OUT_DIR / f"{label}_R1.{ext}", OUT_DIR / f"{label}_R2.{ext}"]


def build_cmd(binary: str, mode: str, label: str, threads: int,
              r1_fq: Path, r2_fq: Path, r1_gz: Path, r2_gz: Path) -> list[str]:
    """Return argv list for the given mode."""
    inp = {"fq": r1_fq, "gz": r1_gz}
    inp2 = {"fq": r2_fq, "gz": r2_gz}
    out_ext = {"fq": "fq", "gz": "fq.gz"}

    mtype, in_fmt, out_fmt = parse_mode(mode)

    if mtype == "stdin":
        return [
            binary, "--stdin", "--stdout",
            "-j", "/dev/null", "-h", "/dev/null",
            "-w", str(threads),
        ]

    ext = out_ext[out_fmt]
    if mtype == "se":
        cmd = [
            binary,
            "-i", str(inp[in_fmt]),
            "-o", str(OUT_DIR / f"{label}_R1.{ext}"),
            "-j", "/dev/null", "-h", "/dev/null",
            "-w", str(threads),
        ]
        if mode in FEAT_DEFS:
            cmd += FEAT_DEFS[mode]["extra"]
        return cmd
    # pe
    cmd = [
        binary,
        "-i", str(inp[in_fmt]), "-I", str(inp2[in_fmt]),
    ]
    if mode == "merge":
        cmd += ["--merged_out", str(OUT_DIR / f"{label}_merged.{ext}")]
    cmd += [
        "-o", str(OUT_DIR / f"{label}_R1.{ext}"),
        "-O", str(OUT_DIR / f"{label}_R2.{ext}"),
        "-j", "/dev/null", "-h", "/dev/null",
        "-w", str(threads),
    ]
    if mode in FEAT_DEFS:
        cmd += FEAT_DEFS[mode]["extra"]
    return cmd


def kill_stale_fastp(binaries: list[str]) -> None:
    """Kill any leftover fastp processes from previous runs."""
    names = {os.path.basename(b) for b in binaries if b}
    names.add("fastp")
    for name in names:
        try:
            r = subprocess.run(["pgrep", "-f", name], capture_output=True, text=True)
            for pid_str in r.stdout.strip().splitlines():
                pid = int(pid_str)
                if pid == os.getpid():
                    continue
                try:
                    os.kill(pid, signal.SIGKILL)
                    print(f"[cleanup] Killed stale process {name} (PID {pid})")
                except ProcessLookupError:
                    pass
        except Exception:
            pass


def run_once(cmd: list[str], mode: str, r1_fq: Path) -> tuple[float, int]:
    """Run a single benchmark iteration, return (wall-clock seconds, peak RSS KB)."""
    stdin_f = None
    kwargs: dict = dict(stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if mode == "stdin-stdout":
        stdin_f = open(r1_fq, "rb")
        kwargs["stdin"] = stdin_f

    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, **kwargs)
    _, status, rusage = os.wait4(proc.pid, 0)
    elapsed = time.perf_counter() - t0
    proc.returncode = os.waitstatus_to_exitcode(status)

    if stdin_f:
        stdin_f.close()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)

    peak_kb = rusage.ru_maxrss
    if sys.platform == "darwin":
        peak_kb //= 1024  # macOS reports bytes
    return elapsed, peak_kb


def run_bench(label: str, binary: str, mode: str, runs: int, threads: int,
              r1_fq: Path, r2_fq: Path, r1_gz: Path, r2_gz: Path) -> tuple[float, int]:
    """Run benchmark N times, print per-run times, return (median time, max peak RSS KB)."""
    cmd = build_cmd(binary, mode, label, threads, r1_fq, r2_fq, r1_gz, r2_gz)
    print(f"[bench] {label} ({mode}, {runs} runs)")

    # warmup
    run_once(cmd, mode, r1_fq)

    times = []
    peak_rss = 0
    for i in range(1, runs + 1):
        for f in OUT_DIR.glob(f"{label}_*"):
            f.unlink()
        t, rss = run_once(cmd, mode, r1_fq)
        times.append(t)
        peak_rss = max(peak_rss, rss)
        print(f"  Run {i}: {t:.2f}s  peak={rss // 1024}MB")

    med = median(times)
    print(f"  >> Median: {med:.4f}s  peak={peak_rss // 1024}MB")
    print()
    return med, peak_rss

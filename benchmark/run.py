"""CLI entry point for fastp benchmark suite."""
import argparse
import gzip
import json
import os
import shutil
import sys
from pathlib import Path

from .modes import (
    ALL_MODES, CPU_COUNT, FEAT_DEFS, MODE_ALIASES,
    DATA_DIR, OUT_DIR, RESULTS_JSON,
    parse_mode, auto_threads,
)
from .datagen import generate_data, generate_merge_data
from .sysinfo import system_info, format_cores
from .runner import run_bench, kill_stale_fastp
from .verify import verify_outputs
from .report import print_summary


def human_size(path: Path) -> str:
    sz = path.stat().st_size
    for unit in ("B", "KB", "MB", "GB"):
        if sz < 1024:
            return f"{sz:.1f} {unit}"
        sz /= 1024
    return f"{sz:.1f} TB"


def banner(text: str, char: str = "=", width: int = 44):
    print(char * width)
    print(f"  {text}")
    print(char * width)


def main():
    parser = argparse.ArgumentParser(description="fastp end-to-end benchmark")
    parser.add_argument("--json", metavar="PATH",
                        help="Load results from JSON and print summary table (skip benchmark)")
    parser.add_argument("--merge", nargs=2, metavar=("ORIG_JSON", "OPT_JSON"),
                        help="Merge two opt-only JSONs into a comparison table")
    parser.add_argument("--mode", default="all",
                        help="Comma-separated modes: fq-fq,fq-gz,gz-fq,gz-gz,"
                             "se-fq-fq,se-fq-gz,se-gz-fq,se-gz-gz,stdin-stdout,"
                             "merge,correction,dedup-pe,dedup-se,qualcut-pe,qualcut-se,ht-pe,"
                             "all,all-pe,all-se,all-feat,everything (default: all)")
    parser.add_argument("--r1", metavar="PATH",
                        help="Custom R1 input (.fq or .fq.gz). Skips data generation.")
    parser.add_argument("--r2", metavar="PATH",
                        help="Custom R2 input (.fq or .fq.gz). Required for PE modes.")
    parser.add_argument("--pairs", type=int, default=10_000_000,
                        help="Number of read pairs (default: 10000000)")
    parser.add_argument("--threads", "-w", type=int, default=0,
                        help="Worker threads (default: auto, cpu_count minus reader/writer per mode)")
    parser.add_argument("--runs", type=int, default=3,
                        help="Repeat count, median reported (default: 3)")
    parser.add_argument("--seed", type=int, default=2026,
                        help="Random seed for data generation (default: 2026)")
    parser.add_argument("--orig", default=None,
                        help="Baseline binary path (default: none, opt-only mode)")
    parser.add_argument("--opt", default="/tmp/fastp_opt",
                        help="Optimized binary path (default: /tmp/fastp_opt)")
    args = parser.parse_args()

    # --- Merge mode: combine two opt-only JSONs into comparison ---
    if args.merge:
        orig_p, opt_p = Path(args.merge[0]), Path(args.merge[1])
        for p in (orig_p, opt_p):
            if not p.exists():
                print(f"File not found: {p}", file=sys.stderr)
                sys.exit(1)
        orig_report = json.loads(orig_p.read_text())
        opt_report = json.loads(opt_p.read_text())
        merged = {
            "system": opt_report.get("system", orig_report.get("system", {})),
            "config": {
                **opt_report["config"],
                "orig": orig_report["config"]["opt"],
            },
            "modes": {},
        }
        for mode in opt_report["modes"]:
            opt_entry = opt_report["modes"][mode]
            orig_entry = orig_report["modes"].get(mode, {})
            merged["modes"][mode] = {
                "threads": opt_entry.get("threads", orig_entry.get("threads")),
                "opt_median": opt_entry["opt_median"],
                "opt_peak_mb": opt_entry.get("opt_peak_mb", 0),
                "orig_median": orig_entry.get("opt_median", 0),
                "orig_peak_mb": orig_entry.get("opt_peak_mb", 0),
                "verify": opt_entry.get("verify", {}),
            }
        print_summary(merged)
        return

    # --- JSON-only mode: load and display ---
    if args.json:
        p = Path(args.json)
        if not p.exists():
            print(f"File not found: {p}", file=sys.stderr)
            sys.exit(1)
        report = json.loads(p.read_text())
        print_summary(report)
        return

    # Expand mode aliases
    raw_modes = args.mode.split(",")
    modes = []
    for m in raw_modes:
        m = m.strip()
        if m == "all":
            modes.extend(ALL_MODES)
        elif m in MODE_ALIASES:
            modes.extend(MODE_ALIASES[m])
        else:
            modes.append(m)
    seen = set()
    modes = [m for m in modes if not (m in seen or seen.add(m))]

    fastp_orig = args.orig
    fastp_opt = args.opt

    if not os.access(fastp_opt, os.X_OK):
        print(f"Missing binary: {fastp_opt}", file=sys.stderr)
        sys.exit(1)
    if fastp_orig and not os.access(fastp_orig, os.X_OK):
        print(f"Missing binary: {fastp_orig}", file=sys.stderr)
        sys.exit(1)
    has_baseline = fastp_orig is not None

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Resolve per-mode thread counts ---
    threads_map: dict[str, int] = {}
    for mode in modes:
        threads_map[mode] = args.threads if args.threads > 0 else auto_threads(mode)

    # --- System info ---
    sysinfo = system_info()

    # --- Header ---
    banner("fastp End-to-End Benchmark")
    if "cpu" in sysinfo:
        print(f"  CPU:     {sysinfo['cpu']}")
    print(f"  Cores:   {format_cores(sysinfo)}")
    mem_str = f"{sysinfo.get('mem_total_gb', '?')}GB"
    if "mem_type" in sysinfo:
        mem_str += f" {sysinfo['mem_type']}"
    if "mem_avail_gb" in sysinfo:
        mem_str += f", {sysinfo['mem_avail_gb']}GB avail"
    print(f"  Memory:  {mem_str}")
    disk = sysinfo.get("disk", {})
    if disk:
        d_parts = [disk.get("model", "?"), disk.get("size", "")]
        proto = disk.get("protocol", "")
        iface = disk.get("interface", "")
        link = disk.get("link_speed", "")
        if proto:
            conn = proto
            if iface and iface != proto:
                conn += f"/{iface}"
            if link:
                conn += f" {link}"
            d_parts.append(conn)
        print(f"  Disk:    {', '.join(p for p in d_parts if p)}")
    if "load_avg" in sysinfo:
        la = sysinfo["load_avg"]
        print(f"  Load:    {la[0]} (1m)  {la[1]} (5m)  {la[2]} (15m)")
    print(f"  OS:      {sysinfo['os']} ({sysinfo['arch']})")
    if args.r1:
        print(f"  Data:    custom ({args.r1})")
    else:
        print(f"  Pairs:   {args.pairs}")
    unique_w = sorted(set(threads_map.values()))
    if len(unique_w) == 1:
        print(f"  Threads: {unique_w[0]}")
    else:
        print(f"  Threads: auto (PE={threads_map.get(modes[0], unique_w[0])}"
              f", SE={threads_map.get('stdin-stdout', unique_w[-1])})")
    print(f"  Runs:    {args.runs} (median reported)")
    print(f"  Modes:   {' '.join(modes)}")
    print(f"  Seed:    {args.seed}")
    print("=" * 44)
    print()

    # --- Data files ---
    custom_data = args.r1 is not None
    need_fq = any(parse_mode(m)[1] == "fq" or m == "stdin-stdout" for m in modes)
    need_gz = any(parse_mode(m)[1] == "gz" for m in modes)
    needs_merge = any(
        m in FEAT_DEFS and FEAT_DEFS[m].get("data") == "merge" for m in modes
    )

    if custom_data:
        r1_src = Path(args.r1)
        r2_src = Path(args.r2) if args.r2 else None
        if not r1_src.exists():
            print(f"File not found: {r1_src}", file=sys.stderr)
            sys.exit(1)
        if r2_src and not r2_src.exists():
            print(f"File not found: {r2_src}", file=sys.stderr)
            sys.exit(1)
        has_pe = any(parse_mode(m)[0] == "pe" for m in modes)
        if has_pe and not r2_src:
            print("ERROR: PE modes require --r2", file=sys.stderr)
            sys.exit(1)

        is_gz = r1_src.name.endswith(".gz")
        if is_gz:
            r1_gz, r2_gz = r1_src, (r2_src or r1_src)
            r1_fq = DATA_DIR / "custom_R1.fq"
            r2_fq = DATA_DIR / "custom_R2.fq"
        else:
            r1_fq, r2_fq = r1_src, (r2_src or r1_src)
            r1_gz = DATA_DIR / "custom_R1.fq.gz"
            r2_gz = DATA_DIR / "custom_R2.fq.gz"

        print(f"[data] Custom R1: {r1_src} ({human_size(r1_src)})")
        if r2_src:
            print(f"       Custom R2: {r2_src} ({human_size(r2_src)})")

        if need_fq and not r1_fq.exists():
            print("[data] Decompressing custom data to plain FASTQ...")
            for gz_p, fq_p in ((r1_gz, r1_fq), (r2_gz, r2_fq)):
                with gzip.open(gz_p, "rb") as fi, open(fq_p, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
        if need_gz and not r1_gz.exists():
            print("[data] Compressing custom data to .gz...")
            for fq_p, gz_p in ((r1_fq, r1_gz), (r2_fq, r2_gz)):
                with open(fq_p, "rb") as fi, gzip.open(gz_p, "wb", compresslevel=1) as fo:
                    shutil.copyfileobj(fi, fo)

        merge_r1_gz = r1_gz
        merge_r2_gz = r2_gz
        if needs_merge:
            print("[data] Merge/correction will use custom data (assumes natural overlap)")
    else:
        r1_gz = DATA_DIR / "bench_R1.fq.gz"
        r2_gz = DATA_DIR / "bench_R2.fq.gz"
        r1_fq = DATA_DIR / "bench_R1.fq"
        r2_fq = DATA_DIR / "bench_R2.fq"

        if r1_gz.exists() and r2_gz.exists():
            print("[data] Reusing existing compressed test data")
            print(f"  R1.gz: {human_size(r1_gz)}")
            print(f"  R2.gz: {human_size(r2_gz)}")
        else:
            print(f"[data] Generating {args.pairs} pairs...")
            generate_data(r1_gz, r2_gz, args.pairs, args.seed)
            print(f"  R1.gz: {human_size(r1_gz)}")
            print(f"  R2.gz: {human_size(r2_gz)}")

        if need_fq:
            if r1_fq.exists() and r2_fq.exists():
                print("[data] Reusing existing uncompressed test data")
            else:
                print("[data] Decompressing to plain FASTQ...")
                for gz_path, fq_path in ((r1_gz, r1_fq), (r2_gz, r2_fq)):
                    with gzip.open(gz_path, "rb") as fi, open(fq_path, "wb") as fo:
                        shutil.copyfileobj(fi, fo)
            print(f"  R1.fq: {human_size(r1_fq)}")
            print(f"  R2.fq: {human_size(r2_fq)}")

        merge_r1_gz = DATA_DIR / "merge_R1.fq.gz"
        merge_r2_gz = DATA_DIR / "merge_R2.fq.gz"
        if needs_merge:
            if merge_r1_gz.exists() and merge_r2_gz.exists():
                print("[data] Reusing existing merge test data (overlapping PE)")
                print(f"  merge_R1.gz: {human_size(merge_r1_gz)}")
                print(f"  merge_R2.gz: {human_size(merge_r2_gz)}")
            else:
                print(f"[data] Generating {args.pairs} overlapping pairs for merge/correction...")
                generate_merge_data(merge_r1_gz, merge_r2_gz, args.pairs, args.seed)
                print(f"  merge_R1.gz: {human_size(merge_r1_gz)}")
                print(f"  merge_R2.gz: {human_size(merge_r2_gz)}")
    print()

    # --- Cache warmup ---
    print("[cache] Warming up filesystem cache...")
    cache_files = [r1_gz, r2_gz]
    if need_fq:
        cache_files += [r1_fq, r2_fq]
    if needs_merge and merge_r1_gz != r1_gz:
        cache_files += [merge_r1_gz, merge_r2_gz]
    for p in cache_files:
        if p.exists():
            with open(p, "rb") as f:
                while f.read(1 << 20):
                    pass
    print()

    # --- Run benchmarks ---
    results: dict[str, dict] = {}

    for mode in modes:
        w = threads_map[mode]
        banner(f"Mode: {mode}  (W={w})", char="-")
        kill_stale_fastp([fastp_orig or "", fastp_opt])
        print()

        # Select data files: merge/correction modes use overlapping PE data
        if mode in FEAT_DEFS and FEAT_DEFS[mode].get("data") == "merge":
            d_r1_fq, d_r2_fq = merge_r1_gz, merge_r2_gz  # merge modes use gz only
            d_r1_gz, d_r2_gz = merge_r1_gz, merge_r2_gz
        else:
            d_r1_fq, d_r2_fq = r1_fq, r2_fq
            d_r1_gz, d_r2_gz = r1_gz, r2_gz

        if has_baseline:
            t_orig, mem_orig = run_bench(
                f"orig_{mode}", fastp_orig, mode, args.runs, w,
                d_r1_fq, d_r2_fq, d_r1_gz, d_r2_gz)
        t_opt, mem_opt = run_bench(
            f"opt_{mode}", fastp_opt, mode, args.runs, w,
            d_r1_fq, d_r2_fq, d_r1_gz, d_r2_gz)
        entry = {
            "threads": w,
            "opt_median": round(t_opt, 4),
            "opt_peak_mb": mem_opt // 1024,
        }
        if has_baseline:
            entry["orig_median"] = round(t_orig, 4)
            entry["orig_peak_mb"] = mem_orig // 1024
        results[mode] = entry

    # --- Verify ---
    verification = verify_outputs(modes, has_baseline)
    for mode in modes:
        results[mode]["verify"] = verification[mode]

    # --- Save JSON ---
    report = {
        "system": sysinfo,
        "config": {
            "pairs": args.pairs,
            "cpus": CPU_COUNT,
            "runs": args.runs,
            "seed": args.seed,
            "orig": fastp_orig,
            "opt": fastp_opt,
        },
        "modes": results,
    }
    RESULTS_JSON.write_text(json.dumps(report, indent=2) + "\n")
    print(f"[json] Results saved to {RESULTS_JSON}")
    print()

    print_summary(report)

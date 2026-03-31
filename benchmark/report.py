"""Summary table generation."""
import subprocess

from .modes import parse_mode, auto_threads


def get_version(binary: str) -> str:
    """Extract version string from fastp binary."""
    try:
        r = subprocess.run([binary, "--version"], capture_output=True, text=True, timeout=5)
        out = (r.stderr + r.stdout).strip()
        for line in out.splitlines():
            if "fastp" in line.lower():
                return line.strip().replace("fastp ", "")
        return out.splitlines()[0] if out else "?"
    except Exception:
        return "?"


def print_summary(report: dict):
    """Print flat markdown table from full report dict."""
    cfg = report["config"]
    results = report["modes"]
    num_pairs = cfg["pairs"]

    si = report.get("system", {})
    has_baseline = cfg.get("orig") is not None
    orig_ver = get_version(cfg["orig"]) if has_baseline else None
    opt_ver = get_version(cfg["opt"])

    from .sysinfo import format_cores
    print()
    print("## fastp Benchmark")
    print()
    cores_str = format_cores(si)
    if "cpu" in si:
        print(f"- **CPU:** {si['cpu']}  **Cores:** {cores_str}")
    else:
        print(f"- **Cores:** {cores_str}")
    mem_parts = []
    if "mem_total_gb" in si:
        mem_parts.append(f"{si['mem_total_gb']}GB")
    if "mem_type" in si:
        mem_parts.append(si["mem_type"])
    if "mem_avail_gb" in si:
        mem_parts.append(f"{si['mem_avail_gb']}GB avail")
    if mem_parts:
        print(f"- **Mem:** {', '.join(mem_parts)}")
    if "load_avg" in si:
        la = si["load_avg"]
        print(f"- **Load:** {la[0]} (1m)  {la[1]} (5m)  {la[2]} (15m)")
    disk = si.get("disk", {})
    if disk:
        disk_parts = [disk.get("model", "?"), disk.get("size", "")]
        proto = disk.get("protocol", "")
        iface = disk.get("interface", "")
        link = disk.get("link_speed", "")
        if proto:
            conn = proto
            if iface and iface != proto:
                conn += f"/{iface}"
            if link:
                conn += f" {link}"
            disk_parts.append(conn)
        print(f"- **Disk:** {', '.join(p for p in disk_parts if p)}")
    if "os" in si:
        print(f"- **OS:** {si['os']} ({si.get('arch', '?')})")
    print(f"- **Pairs:** {num_pairs:,}  **Runs:** {cfg['runs']}  **Seed:** {cfg['seed']}")
    if has_baseline:
        print(f"- **Orig:** `{cfg['orig']}` ({orig_ver})")
    print(f"- **Opt:**  `{cfg['opt']}` ({opt_ver})")
    print()

    # build rows
    global_w = cfg.get("threads")
    rows = []
    for mode, entry in results.items():
        w = entry.get("threads", global_w or auto_threads(mode))
        t_opt = entry["opt_median"]
        tp_opt = num_pairs / t_opt / 1e6
        mem_opt = entry.get("opt_peak_mb", 0)
        mtype = parse_mode(mode)[0]
        io_type = "PE" if mtype == "pe" else "SE"

        verify = entry.get("verify", {})
        v_status = verify.get("status", "n/a")
        v_tag = {"pass": "PASS", "fail": "FAIL", "ok": "OK"}.get(v_status, "n/a")

        if has_baseline:
            t_orig = entry["orig_median"]
            speedup = t_orig / t_opt
            tp_orig = num_pairs / t_orig / 1e6
            saved = t_orig - t_opt
            mem_orig = entry.get("orig_peak_mb", 0)
            mem_ratio = f"{mem_opt / mem_orig:.2f}x" if mem_orig else "-"
            rows.append((mode, io_type, w, t_orig, t_opt, speedup, tp_orig, tp_opt, saved,
                          mem_orig, mem_opt, mem_ratio, v_tag))
        else:
            rows.append((mode, io_type, w, t_opt, tp_opt, mem_opt, v_tag))

    if has_baseline:
        hdr = ("Mode", "Type", "W", "Orig (s)", "Opt (s)", "Speedup",
               "Orig (M/s)", "Opt (M/s)", "Saved (s)",
               "Orig MB", "Opt MB", "Mem", "Verify")
        fmt = [
            lambda r: r[0],
            lambda r: r[1],
            lambda r: str(r[2]),
            lambda r: f"{r[3]:.2f}",
            lambda r: f"{r[4]:.2f}",
            lambda r: f"{r[5]:.2f}x",
            lambda r: f"{r[6]:.2f}",
            lambda r: f"{r[7]:.2f}",
            lambda r: f"{r[8]:.2f}",
            lambda r: str(r[9]) if r[9] else "-",
            lambda r: str(r[10]) if r[10] else "-",
            lambda r: r[11],
            lambda r: r[12],
        ]
    else:
        hdr = ("Mode", "Type", "W", "Time (s)", "M/s", "MB", "Verify")
        fmt = [
            lambda r: r[0],
            lambda r: r[1],
            lambda r: str(r[2]),
            lambda r: f"{r[3]:.2f}",
            lambda r: f"{r[4]:.2f}",
            lambda r: str(r[5]) if r[5] else "-",
            lambda r: r[6],
        ]

    cols = []
    for i, h in enumerate(hdr):
        w = len(h)
        for row in rows:
            w = max(w, len(fmt[i](row)))
        cols.append(w)

    def fmt_row(cells):
        parts = []
        for i, c in enumerate(cells):
            if i == 0:
                parts.append(f" {c:<{cols[i]}} ")
            else:
                parts.append(f" {c:>{cols[i]}} ")
        return "|" + "|".join(parts) + "|"

    print(fmt_row(hdr))
    seps = ["-" * (cols[i] + 2) for i in range(len(hdr))]
    print("|" + "|".join(seps) + "|")
    for row in rows:
        cells = [fmt[i](row) for i in range(len(hdr))]
        print(fmt_row(cells))

    any_fail = False
    for mode, entry in results.items():
        verify = entry.get("verify", {})
        if verify.get("status") == "fail":
            if not any_fail:
                print()
                print("**Verification failures:**")
                any_fail = True
            for fname, fv in verify.get("files", {}).items():
                if not fv["match"]:
                    print(f"- {mode} {fname}: orig=`{fv['orig_md5'][:16]}..` opt=`{fv['opt_md5'][:16]}..`")

    print()

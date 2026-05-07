"""Mode definitions, constants, and path configuration."""
import os
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────

BENCH_DIR = Path("/tmp/fastp_benchmark")
DATA_DIR = BENCH_DIR / "data"
OUT_DIR = BENCH_DIR / "output"
RESULTS_JSON = BENCH_DIR / "results.json"

CPU_COUNT = os.cpu_count() or 4

# ── Mode definitions ─────────────────────────────────────────────────────────

PE_MODES = ["fq-fq", "fq-gz", "gz-fq", "gz-gz"]
SE_MODES = ["se-fq-fq", "se-fq-gz", "se-gz-fq", "se-gz-gz"]
ALL_MODES = PE_MODES + SE_MODES + ["stdin-stdout"]

# Feature test modes: exercise specific fastp features beyond default filtering
FEAT_MODES = [
    "merge", "correction", "dedup-pe", "dedup-se",
    "qualcut-pe", "qualcut-se", "ht-pe",
]

FEAT_DEFS = {
    "merge":      {"base": "pe", "in_fmt": "gz", "out_fmt": "gz", "data": "merge",
                   "extra": ["--merge", "--correction"]},
    "correction": {"base": "pe", "in_fmt": "gz", "out_fmt": "gz", "data": "merge",
                   "extra": ["--correction"]},
    "dedup-pe":   {"base": "pe", "in_fmt": "gz", "out_fmt": "gz",
                   "extra": ["--dedup"]},
    "dedup-se":   {"base": "se", "in_fmt": "gz", "out_fmt": "gz",
                   "extra": ["--dedup"]},
    "qualcut-pe": {"base": "pe", "in_fmt": "gz", "out_fmt": "gz",
                   "extra": ["--cut_front", "--cut_tail", "-W", "4", "-M", "20"]},
    "qualcut-se": {"base": "se", "in_fmt": "gz", "out_fmt": "gz",
                   "extra": ["--cut_front", "--cut_tail", "-W", "4", "-M", "20"]},
    "ht-pe":      {"base": "pe", "in_fmt": "gz", "out_fmt": "gz",
                   "extra": [], "threads": 32},
}

MODE_ALIASES = {
    "all-pe": PE_MODES,
    "all-se": SE_MODES + ["stdin-stdout"],
    "all-feat": FEAT_MODES,
    "everything": ALL_MODES + FEAT_MODES,
}


def parse_mode(mode: str) -> tuple[str, str, str]:
    """Parse mode string -> (type, in_fmt, out_fmt).

    Returns ("stdin", "", "") for stdin-stdout,
            ("se", "fq"|"gz", "fq"|"gz") for se-* modes,
            ("pe", "fq"|"gz", "fq"|"gz") for bare fq-fq etc.
    Feature modes (merge, correction, etc.) resolve via FEAT_DEFS.
    """
    if mode in FEAT_DEFS:
        d = FEAT_DEFS[mode]
        return (d["base"], d["in_fmt"], d["out_fmt"])
    if mode == "stdin-stdout":
        return ("stdin", "", "")
    if mode.startswith("se-"):
        parts = mode[3:].split("-")
        return ("se", parts[0], parts[1])
    parts = mode.split("-")
    return ("pe", parts[0], parts[1])


def auto_threads(mode: str) -> int:
    """Calculate worker threads, reserving reader/writer for each mode.

    PE:  2 readers (L+R) + 2 writers (L+R) = reserve 4
    SE:  1 reader + 1 writer               = reserve 2
    Feature modes may override thread count (e.g. ht-pe forces 32).
    """
    if mode in FEAT_DEFS and "threads" in FEAT_DEFS[mode]:
        return FEAT_DEFS[mode]["threads"]
    mtype = parse_mode(mode)[0]
    reserved = 4 if mtype == "pe" else 2
    return max(1, CPU_COUNT - reserved)

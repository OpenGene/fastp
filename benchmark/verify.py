"""Output verification via MD5 comparison."""
import gzip
import hashlib
from pathlib import Path

from .runner import output_files


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_gz_content(path: Path) -> str:
    """MD5 of decompressed content (compares FASTQ content, not gz bytes)."""
    h = hashlib.md5()
    with gzip.open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fq_md5(path: Path) -> str:
    """MD5 of FASTQ content -- decompress if .gz, otherwise hash directly."""
    if path.suffix == ".gz":
        return md5_gz_content(path)
    return md5_file(path)


def verify_outputs(modes: list[str], has_baseline: bool = True) -> dict[str, dict]:
    """Verify output FASTQ content (decompress gz first). Returns per-mode verification dict."""
    if has_baseline:
        print("[verify] Computing FASTQ content md5 (decompressing gz if needed)...")
    else:
        print("[verify] Computing opt output md5 (no baseline to compare)...")
    verification: dict[str, dict] = {}

    for mode in modes:
        if mode == "stdin-stdout":
            verification[mode] = {"status": "skipped", "reason": "no output files"}
            print(f"  {mode}: skipped (no output files)")
            continue

        opt_files = output_files(f"opt_{mode}", mode)
        mode_result: dict = {"files": {}}

        if has_baseline:
            orig_files = output_files(f"orig_{mode}", mode)
            all_match = True
            for f_orig, f_opt in zip(orig_files, opt_files):
                name = f_orig.name.split("_", 2)[-1]
                h_orig = fq_md5(f_orig)
                h_opt = fq_md5(f_opt)
                match = h_orig == h_opt
                if not match:
                    all_match = False
                mode_result["files"][name] = {
                    "orig_md5": h_orig,
                    "opt_md5": h_opt,
                    "match": match,
                }
                status = "MATCH" if match else "DIFFER"
                print(f"  {mode} {name}: {status}  orig={h_orig[:12]}  opt={h_opt[:12]}")
            mode_result["status"] = "pass" if all_match else "fail"
        else:
            for f_opt in opt_files:
                name = f_opt.name.split("_", 2)[-1]
                h_opt = fq_md5(f_opt)
                mode_result["files"][name] = {"opt_md5": h_opt}
                print(f"  {mode} {name}: md5={h_opt[:12]}")
            mode_result["status"] = "ok"

        verification[mode] = mode_result

    n_pass = sum(1 for v in verification.values() if v.get("status") == "pass")
    n_check = sum(1 for v in verification.values() if v.get("status") in ("pass", "fail"))
    if n_check > 0:
        print(f"  Passed: {n_pass} / {n_check} modes")
    print()
    return verification

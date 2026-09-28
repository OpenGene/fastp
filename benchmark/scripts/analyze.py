#!/usr/bin/env python3
"""Summarize benchmark/results/*.csv: median wall-clock per (version, dataset,
threads), flags HUNG results, and flags output-digest mismatches within a
(dataset, threads) group (a real correctness regression, not just a perf one).
Usage: analyze.py <results.tsv>
"""
import csv
import statistics
import sys
from collections import defaultdict

path = sys.argv[1] if len(sys.argv) > 1 else "results/results.tsv"
rows = list(csv.DictReader(open(path), delimiter="\t"))

times = defaultdict(list)
hangs = defaultdict(int)
digests = defaultdict(set)
for r in rows:
    key = (r["version"], r["dataset"], int(r["threads"]))
    if r["result"] == "OK":
        times[key].append(float(r["secs"]))
        digests[(r["dataset"], int(r["threads"]))].add(r["digest"])
    elif r["result"] == "HUNG":
        hangs[key] += 1

versions = sorted({k[0] for k in times} | {k[0] for k in hangs}, key=lambda v: (v.startswith("fix"), v))
datasets = sorted({k[1] for k in times} | {k[1] for k in hangs})
threads = sorted({k[2] for k in times} | {k[2] for k in hangs})

print("## Median wall-clock seconds (HUNG = timed out)\n")
for ds in datasets:
    print(f"### {ds}\n")
    print("| version | " + " | ".join(f"-w {t}" for t in threads) + " |")
    print("|---|" + "---|" * len(threads))
    for v in versions:
        cells = []
        for t in threads:
            k = (v, ds, t)
            if hangs.get(k):
                cells.append(f"**HUNG** ({hangs[k]}/{hangs[k]+len(times.get(k, []))})")
            elif k in times:
                cells.append(f"{statistics.median(times[k]):.1f}s")
            else:
                cells.append("-")
        print(f"| {v} | " + " | ".join(cells) + " |")
    print()

print("## Output-digest agreement per (dataset, threads)\n")
print("Every version that completed for a given (dataset, threads) should produce")
print("the same digest -- fastp's thread count must not change trimming output.\n")
mismatches = [(k, v) for k, v in digests.items() if len(v) > 1]
if mismatches:
    print("**MISMATCHES FOUND:**\n")
    for k, v in mismatches:
        print(f"- {k}: {len(v)} distinct digests -- {v}")
else:
    print("No digest mismatches across any version/thread/dataset combination that completed.")

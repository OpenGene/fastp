# fastp benchmark corpus & release backfill

A small, reproducible benchmark corpus (see [corpus.md](corpus.md)) plus scripts
to sweep it across fastp releases and thread counts, so performance and
correctness regressions (like the thread-count deadlock fixed in #723) show up
in data instead of being rediscovered by users in production.

## Quick start

```bash
# 1. Fetch/generate the corpus (only needs to be done once)
bash scripts/fetch_corpus.sh data/
python3 scripts/gen_synth_pe.py --out data/synth

# 2. Populate bin/ with the fastp binaries you want to compare, named to match
#    scripts/versions.tsv (prebuilt releases: http://opengene.org/fastp/fastp.<version>)
mkdir -p bin
for v in 1.1.0 1.2.0 1.3.0 1.3.2 1.3.7; do
  curl -sfL -o bin/fastp-$v "http://opengene.org/fastp/fastp.$v" && chmod +x bin/fastp-$v
done

# 3. Run the sweep (thread counts and timeout are overridable; see script header)
bash scripts/run_backfill.sh bin data results/results.tsv

# 4. Summarize
python3 scripts/analyze.py results/results.tsv
```

## What's checked out for each run

For every (version, dataset, thread count): wall-clock time, and an
order-independent digest of the decompressed output (sorted FASTQ records,
md5). A version that hangs is recorded as `HUNG` rather than blocking the
sweep (relevant for the several fastp releases that deadlock above 32
threads — see [#721](https://github.com/OpenGene/fastp/issues/721) /
[#723](https://github.com/OpenGene/fastp/pull/723)). A digest that disagrees
with other versions/thread-counts on the same dataset is a correctness
regression, not just a speed one — `analyze.py` flags this explicitly.

## Backfilled results (this PR)

[`results/backfill-2026-09.tsv`](results/backfill-2026-09.tsv) is a first
real run of this sweep: every non-prerelease v1.x tag with a build available,
plus the [#723](https://github.com/OpenGene/fastp/pull/723) fix branch, across
all three corpus datasets at `-w 1/16/32/48`, on a 48-vCPU AMD EPYC (n2d-highmem-48).
Version selection deliberately keeps every v1.3.x patch rather than collapsing
to "latest patch per minor" — see the comment in `scripts/versions.tsv` for why.
[`results/backfill-2026-09-summary.md`](results/backfill-2026-09-summary.md)
is `analyze.py`'s output on that run.

## Extending this

- Add a dataset: pick a verified ENA/SRA accession, add a `fetch` line to
  `fetch_corpus.sh`, document it in `corpus.md`.
- Add a version to the backfill: add a row to `versions.tsv` and drop the
  matching binary in `bin/`.
- CI: this isn't wired into GitHub Actions yet. The full corpus (~1GB
  compressed, tens of minutes per full sweep at a 12-version x 4-thread x
  3-dataset matrix) is too heavy for every PR, but a 1-dataset x 2-thread
  smoke sweep against the current build would be cheap enough to run on
  every push and would have caught the #721 regression the day v1.3.4 shipped.

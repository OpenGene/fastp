## Median wall-clock seconds (HUNG = timed out)

### atac

| version | -w 1 | -w 16 | -w 32 | -w 48 |
|---|---|---|---|---|
| v1.0.1 | 35.3s | 8.5s | 9.1s | 9.9s |
| v1.1.0 | 35.9s | 8.3s | 8.8s | 10.0s |
| v1.2.0 | 29.0s | 9.2s | 9.5s | 9.8s |
| v1.3.0 | 29.2s | 5.9s | 48.8s | **HUNG** (1/1) |
| v1.3.1 | 28.9s | 6.0s | 46.0s | **HUNG** (1/1) |
| v1.3.2 | 28.5s | 6.0s | 5.9s | 6.5s |
| v1.3.3 | 28.8s | 5.5s | 5.9s | 6.3s |
| v1.3.4 | 28.8s | 5.9s | 47.0s | **HUNG** (1/1) |
| v1.3.5 | 28.6s | 6.2s | 46.7s | **HUNG** (1/1) |
| v1.3.6 | 28.6s | 6.0s | 45.6s | **HUNG** (1/1) |
| v1.3.7 | 28.2s | 5.8s | 45.4s | **HUNG** (1/1) |
| fix-#723 | 28.7s | 6.0s | 6.0s | 6.9s |

### synth

| version | -w 1 | -w 16 | -w 32 | -w 48 |
|---|---|---|---|---|
| v1.0.1 | 11.1s | 3.9s | 4.1s | 4.2s |
| v1.1.0 | 11.2s | 4.0s | 4.1s | 4.2s |
| v1.2.0 | 8.1s | 4.0s | 4.1s | 4.1s |
| v1.3.0 | 8.1s | 2.6s | 11.3s | **HUNG** (1/1) |
| v1.3.1 | 8.3s | 2.6s | 11.4s | **HUNG** (1/1) |
| v1.3.2 | 8.1s | 2.7s | 2.7s | 2.8s |
| v1.3.3 | 8.1s | 2.6s | 2.6s | 2.8s |
| v1.3.4 | 8.1s | 2.6s | 11.3s | **HUNG** (1/1) |
| v1.3.5 | 8.1s | 2.8s | 11.3s | **HUNG** (1/1) |
| v1.3.6 | 8.1s | 2.8s | 11.2s | **HUNG** (1/1) |
| v1.3.7 | 8.1s | 2.7s | 11.2s | **HUNG** (1/1) |
| fix-#723 | 8.2s | 2.7s | 2.8s | 2.9s |

### wgs

| version | -w 1 | -w 16 | -w 32 | -w 48 |
|---|---|---|---|---|
| v1.0.1 | 52.1s | 31.1s | 30.8s | 31.3s |
| v1.1.0 | 52.1s | 28.7s | 28.5s | 30.1s |
| v1.2.0 | 37.5s | 28.9s | 29.9s | 30.1s |
| v1.3.0 | 39.3s | 26.0s | 59.1s | **HUNG** (1/1) |
| v1.3.1 | 38.7s | 26.0s | 60.8s | **HUNG** (1/1) |
| v1.3.2 | 37.4s | 25.2s | 26.9s | 26.0s |
| v1.3.3 | 37.2s | 24.8s | 25.3s | 25.8s |
| v1.3.4 | 36.2s | 25.1s | 58.4s | **HUNG** (1/1) |
| v1.3.5 | 36.9s | 24.7s | 58.8s | **HUNG** (1/1) |
| v1.3.6 | 37.3s | 24.6s | 57.1s | **HUNG** (1/1) |
| v1.3.7 | 37.0s | 24.9s | 59.0s | **HUNG** (1/1) |
| fix-#723 | 36.6s | 25.2s | 25.9s | 26.0s |

## Output-digest agreement per (dataset, threads)

Every version that completed for a given (dataset, threads) should produce
the same digest -- fastp's thread count must not change trimming output.

No digest mismatches across any version/thread/dataset combination that completed.

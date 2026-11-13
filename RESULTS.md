# Technical smoke-test result

Date: 2026-11-13. Status: **PASS** for the bounded technical pipeline.

This is a public PBMC3k blood-cell example, not a human lung/IPF cohort. It does not establish a cell type, a disease mechanism, a spatial niche, cross-donor reproducibility, or large-dataset capacity.

## Verified execution

| Item | Measured result |
|---|---|
| Environment | Independent virtual environment; Python 3.10.20, Scanpy 1.11.5 |
| Input | 2,700 cells x 32,738 genes, raw counts |
| After demonstration QC | 2,698 cells x 13,714 genes |
| Variable genes / PCA | 2,000 variable genes; 30 components |
| Clustering | Eight unannotated Leiden clusters; no biological labels assigned |
| Analysis wall time | 57.143 seconds, including scientific imports and first-run compilation |
| Whole successful job | 250.561 seconds, including acquisition and input checks; package installation excluded |
| Peak analysis-process memory | 819.105 MiB; this is peak RSS, not reserved memory |
| Controls | Two computational threads; 16 GiB virtual-address-space ceiling; low process priority |
| Entire test-root disk use | 1,027,760 KiB, approximately 0.98 GiB, including environment, caches, input and results at measurement |
| Execution status | Exit code 0 and validated `SUCCESS.json` |

The example retained 2,698 cells under tutorial-only cutoffs. This is not a validated QC rule for IPF and is not a biological claim about the eight clusters.

## Checks that passed

1. Environment dependency resolution and core imports.
2. Input H5AD hash matching the acquisition record; expected dimensions and unique cell identifiers.
3. Finite, nonnegative, integer-like input counts.
4. Finite PCA/UMAP embeddings.
5. Reloaded output retains cell/gene identifiers, shape, filtered raw counts, cluster assignments, PCA and UMAP exactly.
6. All saved artifact hashes match on the server; retrieved artifacts match the completion manifest.
7. Executed acquisition, pipeline, runner and requirements hashes match the reviewed local source files.
8. QC and UMAP PNGs were visually inspected for readable labels and intact layout. No cell-type or IPF interpretation was applied.

## Acquisition correction and retained failure

Attempt `pilot_01` exited with code 1 at data acquisition: the original 10x archive endpoint returned HTTP 403. No QC or clustering ran in that attempt. The failure status and original acquisition script remain on the server.

Attempt `pilot_02` used the raw H5AD explicitly referenced by the official Scanpy 1.11.5 dataset implementation, downloaded unchanged from `https://falexwolf.de/data/pbmc3k_raw.h5ad`. Its 5,855,727 bytes have SHA-256 `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1`. The exact URL, pinned source-definition URL, acquisition time, dimensions and hash are in [acquisition.json](results/pilot_02/acquisition.json).

The corrected source was accepted for this technical example only. No dataset selection or biological threshold for Stage 6 was inferred from it.

## Reviewable outputs

- [QC figure](results/pilot_02/qc.png) and [UMAP figure](results/pilot_02/umap.png).
- [Run summary and saved-object checks](results/pilot_02/run_summary.json).
- [Execution status and source hashes](results/pilot_02/execution.json).
- [Completion manifest](results/pilot_02/SUCCESS.json) and [server verification](results/pilot_02/server_verification.json).
- [Cell-level technical results](results/pilot_02/cell_qc_clusters.csv), [run log](results/pilot_02/run.log), and [first-attempt failure](results/pilot_02/failed_attempt_01.json).
- [Resolved package versions](environment.freeze.txt) and [environment validation](ENVIRONMENT.md).

The raw and processed H5AD checkpoints remain under `/home/data/t120631/project/ipf_stage6_smoke_test_20261113/runs/pilot_02/` on the server. They are intentionally absent from the small Git repository; their hashes are retained in the manifest.

## What this establishes for the project

The Mac Codex -> SSH -> isolated Linux environment -> persistent background job -> saved checkpoints -> verified result return loop works for a small real single-cell example. The job completed in tmux without requiring a continuously attached SSH command.

Next, select and inspect suitable IPF dataset metadata before defining biological analysis settings. Spatial processing, independent-cohort validation, state annotation, large-file throughput, GPU methods and large-dataset memory requirements remain untested. No formal IPF discovery analysis or experiment has started.

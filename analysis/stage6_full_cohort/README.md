# Full native cohort preparation and human review

This execution implements the [approved preparation specification](../stage6_specification/README.md). It prepares evidence for human QC, identity and anatomical decisions. It does not lock biological eligibility or perform niche/target discovery.

The run is `/home/data/t120631/project/ipf_stage6_pilot/runs/preparation_001`. Original public inputs remain under `/home/data/t120631/project/ipf_stage6_pilot/data/raw`; the earlier pilot is preserved. Large raw counts, per-cell review tables and original full-resolution H&E stay on the server.

## Native intake

| Resource | Native inputs | Original source profiles | Role restriction |
|---|---:|---:|---|
| GSE135893 | 43 libraries | 114,396 cells | IPF/control discovery; other ILD separate; capture provenance remains a proxy where unresolved |
| GSE159354 | 48 libraries | 184,019 cells | Preserve donor and FACS fraction; matched ordinary EPCAM+/Triple− comparisons start at at most 3 IPF versus 3 control donors each |
| GSE159585 | 18 samples | 148,154 nuclei | Separate robustness component; disease and procurement fully confounded |
| Original Franzén | 25 sections | 124,800 raw capture positions | All-section anatomy review; prioritize within-IPF interfaces; HC comparison retains procurement confounding |

These are measurement units and source-retained profiles, not independent donor counts or post-QC inclusion totals. GSE136831/GSE122960 candidate expression outcomes remain reserved.

## Execution and evidence

1. `prepare_execution.py` derives the full acquisition manifest from locked intake. `bootstrap_full.py` creates a new isolated run once. Do not run the bootstrap against an existing run.
2. The existing pilot acquisition tool downloads/verifies complete files. `adopt_download_relay.py` can adopt exact-manifest, size- and SHA-verified client transfers while preserving interrupted server partials and prior records.
3. `prepare_rna.py` writes one original-count checkpoint per native RNA input and verifies count/axis reload equality. It does not filter, normalize, merge symbols or integrate resources.
4. `rna_review.py` generates per-capture MAD sensitivity flags, doublet diagnostics, fixed marker-seed evidence, original-label comparisons and grouped coverage. All assignments and flags require adjudication; unknown/conflicting profiles remain visible. Ambient correction is not forced without adequate background input.
5. `spatial_review.py` checks all 25 raw/filtered matrices, coordinates and original H&E mappings, and creates whole-field, fixed-grid and registration review images. Technical checks do not certify a target anatomical interface.
6. `genes/` fixes official reference assets and produces row-preserving gene-identity/coverage reports. Symbol-only correspondences retain provisional status. Mapping tables do not assert that count matrices have been merged or harmonized.
7. `validate_full_run.py` independently checks the completed native preparation. `build_review_pack.py` requires exact intake coverage and creates an offline review package with local-link/image checks and file hashes.

Scientific jobs use the existing Python environment at `/home/data/t120631/project/ipf_stage6_smoke_test_20261113/.venv/bin/python`, four threads per process, and isolated additional packages under the run's `dependencies/`. The installation report is retained; the existing environment is not modified.

## Human review

Open `results/review_pack/index.html` after verified export. Review joint QC and marker/source-label evidence by capture, then H&E and geometry for every section before examining candidate expression. Save reviewed copies of the CSV templates separately, with reviewer, date, decision and reason. Record uncertain regions as well as eligible interfaces.

The generated originals deliberately retain `not_adjudicated` / `unreviewed`. Per-cell masks and the corresponding source-qualified IDs remain on the server. Low cell coverage affects a specified endpoint; it is not automatic whole-donor removal. Expert QC/identity decisions, anatomical ROI masks, cross-source independence and endpoint-specific estimability are required before final cohort allocation and Checkpoint A.

See `RESULTS.md` for the verified completion state, exceptions and export links.

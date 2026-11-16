# Full native preparation: verified review material

Completed 2026-11-16 at `/home/data/t120631/project/ipf_stage6_pilot/runs/preparation_001`.

**Status: preparation evidence complete; human QC, identity, anatomy and final endpoint eligibility pending.** No profiles were excluded, no count matrices were integrated, and no disease/state/niche/target discovery was performed. External reserved candidate expression outcomes were not inspected.

## Complete source and native-input coverage

All 58 declared source files passed complete-file acquisition/SHA verification: **3,993,968,410 bytes**. Interrupted downloads were preserved; eight complete public-source client downloads were relayed over SSH with equal sizes and hashes. Original acquisition records and previous attempts remain traceable.

| RNA resource | Native libraries/samples | Original profiles retained | 3-MAD review flags | 4-MAD | 5-MAD |
|---|---:|---:|---:|---:|---:|
| GSE135893 | 43 | 114,396 | 6,727 | 3,695 | 2,184 |
| GSE159354 | 48 | 184,019 | 15,732 | 9,935 | 7,148 |
| GSE159585 | 18 | 148,154 | 21,473 | 15,789 | 12,348 |
| Total measurement units | 109 | 446,569 | 43,932 | 29,419 | 21,680 |

All 109 original-count checkpoints passed parser/reload count and axis comparisons; each was linked by SHA to its review input. No zero-total profiles were present. Flags are per-capture diagnostic queues, not final exclusions. Cells, nuclei, sorted fractions and repeated libraries must not be interpreted as independent donors.

All 109 inputs obtained Scrublet scores and automatic diagnostic thresholds, with zero unscored profiles. An initial Scanpy/SciPy neighbor compatibility error affected larger captures. A minimal change to the official `use_approx_neighbors=False` option passed a real 7,114-nucleus test; all three resources were uniformly rerun with exact Euclidean neighbors. The initial output remains at `rna_review_initial/`. The existing environment was not modified. Thresholds and 0.02/0.05/0.10 prior sensitivity require review; these priors are not observed doublet rates or removal quotas.

Provisional broad marker assignment left 5,689 / 15,857 / 19,567 profiles unknown in GSE135893 / GSE159354 / GSE159585 respectively. Unknown and competing programs remain represented. Source labels, multi-marker evidence and uncertainty are available for human adjudication; these assignments are not expert fine-state annotation or lineage inference.

## Spatial and gene evidence

All **25 sections from eight source donors** passed technical parsing/linkage checks, preserving **124,800 raw capture positions**, including **101,907 author in-tissue positions**. Every section has exact source H&E mapping, whole-field review image, fixed 3×3 grid, coordinate overlay and raw-count/feature maps. The whole SpaceRanger ZIP passed CRC. All section/ROI/mask expert states remain unreviewed. Technical registration and image readability do not certify an airway interface.

All **134 native RNA/spatial units** have a verified feature-axis crosswalk. There are four distinct ordered feature axes and fixed official Ensembl116/HGNC reference assets. The shared unique candidate set including provisional symbol correspondences contains **18,725 genes**. Strict shared identity is zero because two exports contain symbols without verified original stable-ID crosswalks; this does not mean there are no shared biological genes. No symbol-based count collapse or derived-matrix projection occurred. CLIC5 ambiguity and the GSE159354 RGS5 duplicate/source-ID conflict are explicitly exposed in the review package. See [gene results](genes/RESULTS.md).

## Verified deliverables

- [Offline human review page](results/review_pack/index.html): 109 RNA and 25 tissue review cards, **539 readable images**, **1,834 local links**, and 1,410 portable files.
- [Review instructions](REVIEW_GUIDE.md): QC/identity, all-section anatomy, ROI inventory and endpoint-specific allocation order.
- [Independent full-run audit](results/execution_evidence/full_run_validation.json): **13,898 passed, zero failed**. It checks exact intake, original H5AD structure and creation/review hash links, profile-level review rows/count totals, spatial/H&E correspondence, gene assets and all portable-file hashes/links. Large count payloads were checked during creation/review, rather than rehashed again in this final audit.
- [RNA diagnostic totals](results/execution_evidence/RNA_diagnostic_totals.json) and [run status](results/execution_evidence/RUN_STATUS.json).
- Source export `RUN/review_export.tar.gz`: **365,932,680 bytes**, SHA-256 `86a6fd92c04a67fe925366faef7b2089f8b27d27bfc1e63fcdbfb4646bb34283`. [Local transfer verification](results/transfer_verification.json) passed the archive hash, all 1,409 non-self-referential portable-file hashes and 1,834 local references. The package manifest is the 1,410th file. Five canonical execution-script hashes also match the local source; representative QC, marker and histology images were visually inspected for readability without assigning biological eligibility.

Full counts, per-cell review tables, complete row-level gene mapping and original full-resolution H&E remain on the server. The portable package has resized review derivatives; zoom into the source image for detailed anatomical adjudication.

## Remaining decisions

Human QC/identity dispositions, reviewed per-cell masks and qualified anatomical ROI masks must be supplied before final endpoint-specific inclusion. Capture provenance proxies, source-reference uncertainty, cross-source donor independence and endpoint estimability remain open where documented. FACS routes remain separate; GSE159585 and Franzén HC/IPF procurement confounding remains part of interpretation. The final biological cohort lock and Checkpoint A are **not passed**.

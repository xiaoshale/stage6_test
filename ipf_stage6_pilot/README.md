# Native IPF complete-file technical pilot

Authorized and started: 2026-11-15. Server root, as named by the user:

`/home/data/t120631/project/ipf_stage6_pilot`

Purpose: acquire the complete original files for the selected readiness package, verify native parsing and donor/library linkage, calculate descriptive technical QC, and save reloadable raw-count checkpoints. This is not fine-state discovery, integration, differential expression, deconvolution, neighborhood testing or target inference.

The readiness assessment in the research project and this [pilot manifest](config/pilot_manifest.json) define the donor/section selection and evidence boundaries. The complete-file [acquisition manifest](config/acquisition_manifest.json) has 24 required files: four Habermann files, 16 native DePianto core/background H5 libraries, two Franzén files and two snRNA files. Four optional unsorted DePianto libraries are outside this first required pass. A complete global count file may contain more donors/cells than the intended QC subset.

## Server storage and execution

```text
ipf_stage6_pilot/
  code/                       deployed acquisition/QC scripts
  config/                     exact pilot and acquisition manifests
  data/raw/<resource>/        complete original files and preserved partials
  data/acquisition_records/   per-file URL, transfer status, bytes and SHA-256
  data/acquisition.json       aggregate acquisition result
  runs/pilot_001/<resource>/  raw-count H5AD, cell/spot QC and source summaries
  runs/pilot_001/source/      immutable script snapshots used by each component
  runs/fixture_*/             explicit small positive/failure examples
  logs/                      independent acquisition and component logs
  environment.json           actual runtime path and package versions
  environment.freeze.txt     resolved dependencies
```

The existing isolated smoke-test environment is reused without changing its packages. This pilot's code, data, temporary files, logs and results stay in the new directory. The earlier PBMC inputs/results are separate.

Before execution, the account quota was checked: approximately 49 GB used against a 500 GB soft quota. The shared filesystem had approximately 27 TB available, but that is not an individual allocation. The account rule remains at most 48 CPU/200 GB across all user tasks. This pilot uses two download workers and four computational threads, a 48 GiB per-component virtual-memory ceiling, and sequential substantive QC components. Measured peak memory/runtime are reported with results; these settings are not a promise of capacity for an eventual full atlas.

Example monitoring from Termius:

```bash
tail -f /home/data/t120631/project/ipf_stage6_pilot/logs/acquisition_20261115.log
tail -f /home/data/t120631/project/ipf_stage6_pilot/logs/GSE135893_pilot_001.log
```

`acquire_inputs.py` preserves partial transfers, verifies byte counts and local hashes, and checks small metadata against the readiness snapshot where available. A generated SHA-256 records the acquired bytes; it is not an independently published repository checksum. Existing completed inputs are rehashed before reuse. Failed transfers and failed run directories are retained.

`run_component.py` checks the required acquisition records, sizes and hashes, snapshots the executed scripts, enforces resource settings, and records the process exit status. A detached tmux session allows a job to survive an SSH/client disconnect. Submission alone is not success.

## Component-specific interpretation

| Component | Selected measurement units | QC/output boundary |
|---|---|---|
| GSE135893 | Four manifest donors, 8,577 exact source-annotated cell barcodes | Read the full original MTX and axes. Preserve the selected raw counts and every source gene. The 105,817 additional global columns are not automatically called cells and do not enter the selected H5AD. |
| GSE159354 | Four donors, eight EPCAM+/Triple-negative plus eight CD31+/CD45+ libraries | Retain each source-filtered H5 in a separate raw checkpoint, namespace barcodes by GSM, and preserve fractions. A FACS fraction is not a verified cell state or an absolute tissue proportion. |
| Franzén | Six original sections from three IPF and three control donors | Retain all 4,992 raw capture positions per section with author in_tissue and filtered-H5 membership. Verify raw/filtered submatrix equality, coordinate orientation, images and scales. Positions are not single cells. |
| GSE159585 | Four selected donors, 28,970 source nuclei | Read the entire compressed gene-by-cell CSV through EOF/CRC using streaming rows; keep selected counts sparse. All IPF are explants and all controls autopsies, so this is a separate technical/measurement robustness component, not equal disease-effect evidence. |

For scRNA/snRNA, detected genes, total counts and mitochondrial proportions use original integer counts. Flags for fewer than 200 detected genes, more than 25% MT-prefix counts, and zero counts are **descriptive inspection flags only**. They do not remove cells and are not adopted as final, cross-platform biological inclusion criteria. Mitochondrial percentage has no defined denominator for zero-count cells and is retained as missing. FFPE spatial measurements have different detection properties; an identical cell-level mitochondrial rule is not imposed on them.

Source annotations are retained with explicit provenance. Broad Mesenchymal/stroma labels and Triple-negative sorting do not establish usable fibroblast states. Formal state recovery, donor-level power, high-resolution airway-interface adjudication, procurement covariates and cross-study donor independence remain follow-up gates. The external GSE136831/GSE122960 measurements and integrated references containing them are not acquired or used for training in this run.

## Verification

Three small hand-specified fixtures cover material format/mapping risks before full runs:

- Habermann: shuffled annotation order, an unannotated same-library column that must be excluded, four donor/library mappings, duplicate gene symbols and a retained zero-count cell.
- snRNA: full cell-ID mapping including LUNG010, an unselected COVID cell, duplicate gene rows, a retained zero-count cell, and rejection of missing global metadata or fractional selected counts.
- Spatial: asymmetric x/y coordinates, raw versus author-filtered counts, all capture positions/flags, two image resolutions, roundtrip preservation, and rejection of a deliberately inconsistent filtered submatrix.

The actual runs validate original values, sparse structure, source axes/metadata and saved H5AD reloads. Final artifact/transfer checks and component outcomes are recorded in the results report after completion. No raw matrices or large H5AD objects are copied to the small Git repository; code, configurations, logs, summaries and small figures can be versioned at verified milestones.

## Primary format documentation

- [10x feature-barcode HDF5 matrix structure](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices): feature-by-barcode CSC, explicit gene/barcode axes.
- [AnnData H5AD writing](https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.write_h5ad.html): file export interface. Actual execution uses the recorded AnnData 0.11.4 installation, not an assumption that the current online documentation version is installed.

The exact source-file URLs remain in the acquisition manifest, with inherited source provenance in the readiness package. Format documentation does not validate IPF biology or universal QC thresholds.

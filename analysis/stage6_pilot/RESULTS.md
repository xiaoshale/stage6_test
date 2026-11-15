# Native IPF complete-file technical pilot results

Completed: **2026-11-15**. Status: **all four declared technical components passed**. The server root is:

`/home/data/t120631/project/ipf_stage6_pilot`

The 24 required source files were completely acquired: **2,836,607,664 bytes (2.84 GB)**. All four components parsed successfully, retained the selected original counts, and passed their saved-object reload checks. The final [completion manifest](results/pilot_001/SUCCESS.json) records **89 passing checks and 24 sparse H5AD checkpoints**. This is a successful technical pilot on the selected subsets, not a final biological cohort lock or a discovery result.

## Actual scope and resource measurements

| Component | Selected people/libraries | Retained measurement units | Source feature rows | Process time | Peak process RSS |
|---|---|---:|---:|---:|---:|
| GSE135893 | Four donors; one source library each | 8,577 source-annotated cells | 33,694 | 36.4 s | 6.94 GB |
| GSE159354 | Four donors; 16 sorted libraries | 55,189 source-filtered library cells | 58,302 per library | 40.6 s | 0.63 GB |
| GSE159585 | Four donors | 28,970 nuclei | 27,714 | 275.3 s | 2.07 GB |
| E-MTAB-14121 / S-BSST1410 | Six sections from six donors, three IPF/three controls | 29,952 raw capture positions; 25,466 author in-tissue positions | 36,601 per section | 29.2 s | 0.57 GB |

GB is decimal. Times are component-process elapsed times, excluding acquisition and final aggregate review. Measurements apply to these actual runs and do not predict the capacity of an eventual integrated atlas. Cells, nuclei, capture positions, repeated libraries and independent people remain distinct units; do not pool their counts as independent biological replicates.

The original global Habermann MTX and snRNA CSV include more material than the selected QC subsets. The full global payloads were read, while QC/checkpoints retain the explicitly selected barcodes. All selected cells/features were retained. The remaining eligible donors, libraries and spatial sections still need their own full-cohort QC and inclusion assessment before formal discovery.

## What was verified

- **Acquisition:** all 24 source records are complete; actual bytes and SHA-256 match the retained files. Small metadata with readiness-snapshot hashes agrees with that snapshot. Initial successful transfers were full HTTP 200 responses; the spatial continuation used the same source entity/If-Range. The complete spatial ZIP passed CRC for every member.
- **Habermann barcode/annotation join:** all author metadata barcodes map to the original axes. The four selected donors contribute 8,577 exact annotated barcodes. Computed nCount/nFeature agree with the corresponding source metadata with zero disagreements. The 105,817 additional global columns remain in the unchanged raw input; they were not automatically called cells or assigned to donors by a guessed prefix.
- **DePianto native fractions:** all 16 complete H5 files parsed, with library-qualified cell IDs and preserved donor/fraction labels. Original integer counts, source features and metadata survive H5AD roundtrip. Libraries have separate checkpoints rather than being silently integrated or zero-filled across potentially different gene references. FACS labels describe collection conditions, not verified fine states or absolute tissue proportions.
- **snRNA:** all 176,540 source cell IDs in the CSV header and metadata agree. Streaming reached EOF and passed gzip CRC, retaining exactly the four selected patients' 28,970 nuclei and all 27,714 source rows. Sparse values/indices, source axes, cell totals and source labels survived reload. The complete file was temporarily downloaded on the Mac after direct server transfer proved slow; equal source/SSH/server hashes and both transfer histories are retained.
- **Spatial:** every selected raw H5 barcode joins the coordinate table; each author-filtered H5 equals its exact raw-count submatrix. All original positions, in-tissue/filtered flags, source features, coordinates, scale factors and both image resolutions survive H5AD roundtrip. In-tissue centers fall within the scaled image bounds, and low/high-resolution dimensions agree with their scales.
- **Transfer to the Mac:** the 35-file small review bundle matches its server size/SHA-256. All 31 exported run artifacts listed in the completion manifest match their recorded hashes. QC-input hashes in the retained summaries agree with acquisition records. Raw files and large H5ADs were not copied into the local review package or Git repository. See [transfer verification](results/transfer_verification.json).

The generated SHA-256 values identify acquired bytes and verify subsequent preservation. They are not represented as independently published repository checksums.

## Spatial section results

| Section | Donor | Diagnosis | Raw positions | Author in-tissue | Author filtered H5 | Roundtrip/registration checks |
|---|---|---|---:|---:|---:|---|
| V10T03-279-C1 | IPF_1 | IPF | 4,992 | 4,882 | 4,882 | Pass |
| V10T31-015-C1 | IPF_2 | IPF | 4,992 | 4,740 | 4,740 | Pass |
| V10T03-281-C1 | IPF_3 | IPF | 4,992 | 3,773 | 3,773 | Pass |
| V10T31-015-A1 | HC_2 | Control | 4,992 | 3,643 | 3,643 | Pass |
| V10T31-019-A1 | HC_3 | Control | 4,992 | 4,393 | 4,393 | Pass |
| V10T03-281-A1 | HC_4 | Control | 4,992 | 4,035 | 4,035 | Pass |

All six raw/filtered count comparisons passed, with no per-section technical warnings. The six low/high-resolution coordinate overlays and four cohort-level QC distribution plots were visually reviewed. Tissue footprints and coordinate orientation agree between resolutions, with no obvious registration reversal or broken plot. This is a gross technical image review, not specialist pathology or proof of an airway epithelial–fibroblast niche.

The inspected images are the native SpaceRanger low/high-resolution images. The separate full-resolution H&E archive is outside this 24-file pilot. Expert anatomical adjudication and appropriate original image resolution remain necessary before a formal airway-interface inclusion rule.

Example review artifacts:

- [IPF coordinate overlay](results/pilot_001/Franzen/V10T03-279-C1/coordinate_overlay.png)
- [Control coordinate overlay](results/pilot_001/Franzen/V10T31-015-A1/coordinate_overlay.png)
- [Spatial technical distributions](results/pilot_001/figures/Franzen_technical_qc.png)
- [Single-cell technical distributions](results/pilot_001/figures/GSE159354_technical_qc.png)
- [Section/donor QC table](results/pilot_001/Franzen/section_donor_qc.tsv)

Histology and coordinates derive from the original Franzén acquisition, [BioStudies S-BSST1410](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410); overlays were generated for technical review. Source sample identities and attribution are preserved.

## Descriptive QC and capture observations

Technical flags used original counts: fewer than 200 detected source features, more than 25% MT-prefix counts, or zero counts. These are inspection references, not validated biological inclusion criteria. In the 16 DePianto libraries, 167 cells carried the low-feature flag and 507 carried the high-mitochondrial flag; the flags can overlap. They were retained. No selected zero-count cells were observed in the real pilot, but the positive/negative fixtures verify their preservation and undefined mitochondrial proportions as missing values.

The fixed 17-marker technical capture audit completed on 18 scRNA/snRNA checkpoints and produced 952 group/marker records. All fixed marker symbols were represented in the relevant source feature axes; nonzero detection and proportions are reported separately. Published broad labels remain published labels; the script assigns no new cell identity and performs no disease comparison. Duplicate symbol rows use a detection union for this audit without merging or rewriting X. See [marker capture](results/pilot_001/capture_markers/marker_capture.tsv) and [its scope](results/pilot_001/capture_markers/summary.json).

The pilot resolves file accessibility, native format parsing and explicit technical linkage. It does **not** establish adequate fine epithelial/fibroblast-state recovery in every donor. For example, source-labelled Mesenchymal cells remain sparse in some Habermann controls; a broad label or Triple-negative fraction is not equivalent to a usable fibroblast state.

## Storage and reproducibility

- Original files: `data/raw/<resource>/` under the server root.
- Raw QC checkpoints and per-cell/position tables: `runs/pilot_001/<resource>/`.
- Source/code snapshots: `runs/pilot_001/source/` and component execution JSON records.
- Figures: `runs/pilot_001/figures/` plus six per-section coordinate overlays.
- Completion record: `runs/pilot_001/SUCCESS.json`.
- Local small reviews: `analysis/stage6_pilot/results/` in the research project.

The existing isolated Python environment was reused without package changes; [environment.json](environment.json) and [resolved dependencies](environment.freeze.txt) record its actual path/versions. Source-specific scripts and fixtures are versioned in the user-designated `xiaoshale/stage6_test` repository. Failed/slow direct-transfer partials and interruption logs are retained as provenance and are not analysis inputs. Nothing in this task changed the frozen Stage 5 evidence or the baseline resource catalog.

An earlier automatic permission-review timeout interrupted monitoring and postponed the spatial run. Following the user's continuation reply, remote access succeeded and all remaining checks completed. The old interruption snapshot in [RESUME.md](RESUME.md) is retained as history and is superseded by this report.

## Biological entry decision

The four resources pass this **selected technical pilot**. Final biological allocation and Checkpoint A remain conditional on:

1. Whole-cohort donor/specimen/library provenance, clinical/sampling covariates and target-cell recovery, without treating repeated regions or sorted fractions as independent people.
2. A documented source-gene-ID harmonization and ambiguity policy; the different feature-row counts are not already harmonized expression matrices.
3. Donor/library-aware final QC, doublet and ambient-RNA strategy appropriate to the available source data, preserving biologically meaningful rare states and sampling differences.
4. Expert airway/airway-like anatomical coverage review across the intended spatial donors/sections. Visium spots remain mixtures and do not directly prove single-cell adjacency.
5. Continued separation of the procurement-confounded snRNA component and preservation of GSE136831/GSE122960 for conditional external tests. No reserved expression outcomes or atlas-derived disease states were used to train this pilot.

The next verifiable work is to specify those full-cohort cleaning/annotation/measurement rules using these observed inputs, then fix Checkpoint A before formal discovery. No normalization, integration, clustering, disease-effect testing, niche reconstruction, ligand–receptor prediction or molecular ranking was performed in this technical pilot.

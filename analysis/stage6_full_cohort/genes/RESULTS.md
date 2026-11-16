# Full native gene identity preparation results

Executed and verified 2026-11-16 for `/home/data/t120631/project/ipf_stage6_pilot/runs/preparation_001`.
Status: complete identity/coverage metadata for human review. Original counting-reference reconstruction, matrix projection and final biological eligibility remain separate gates.

## Coverage of the locked intake

All **134 native units** have a verified feature mapping crosswalk: 43 GSE135893 libraries, 48 GSE159354 native libraries, 18 GSE159585 samples and 25 Franzén sections.
The 28 physically distinct feature-axis files represent **four distinct axes by ordered content**. Each declared raw feature checksum, axis length, unit identity and donor crosswalk was checked. Neither source row order nor feature symbols were rewritten for content comparison.

| Resource | Source rows | Strict source-ID eligible rows | Symbol-only provisional rows | Source rows sharing another row's canonical ID |
|---|---:|---:|---:|---:|
| GSE135893 | 33,694 | 0 | 21,730 | 155 |
| GSE159354 | 58,302 | 56,079 | 0 | 0 |
| GSE159585 | 27,714 | 0 | 19,786 | 109 |
| Franzén | 36,601 | 35,406 | 0 | 0 |

The zero strict sets for the two symbol-only exports reflect unresolved original ID/reference provenance, not a lack of shared biological genes. Ambiguous, retired, unmapped and colliding rows remain in the ledgers and original objects; they were not aggregated or discarded from count matrices.

The source GSE159354 axis has **1,654 repeated nonblank symbol labels**, with unique source feature IDs. Franzén has 10 repeated symbol labels. GSE135893 and GSE159585 have zero repeated source labels. Synthetic duplicate tests are separate from these real source facts.

The cross-resource candidate union is 59,101 canonical candidates, including held candidates. The intersection of unique source-ID or provisional symbol correspondences contains **18,725 genes**. The strict intersection is zero under the declared source-provenance rule. These are feature-identity coverage sets, not measured nonzero expression, sufficient state/program coverage or approved final inference features.

## Fixed identity-marker audit

All 62 fixed RNA review seed labels occur in each native source axis. Mapping status adds two specific cautions:

- **CLIC5**: its source label is present in GSE135893 and GSE159585, but frozen HGNC approved/alias/history mappings are not unique under the conservative rule. Its identity mapping remains held; a source-label match does not establish an unambiguous cross-resource gene test.
- **RGS5**: GSE159354 has two rows with this label (`ENSG00000143248` and `ENSG00000232995`); one source-ID/symbol mapping conflicts. The marker is recorded as multiple-source-row review. No symbol-based gene collapse was performed.

Franzén has strict source-ID mapping for all 62 markers; the two symbol-only exports each have 61 provisional markers plus CLIC5 held; GSE159354 has 61 single-source-ID markers plus RGS5 multiple-row review. These findings do not remove cells, prove marker specificity or constitute expert cell annotation.

## Artifacts and transfer verification

Server results reside under `RUN/genes/mapping_001/`:

- `batch_mapping_summary.json`: actual intake counts, four axis summaries, source/reference checksums and coverage accounting.
- `unit_feature_mapping_crosswalk.tsv`: all 134 native units linked to raw feature axes and mapping files.
- `axes/*/gene_mapping.tsv` and `mapping_summary.json`: complete row-level identity evidence and holds.
- `coverage/`: common strict and unique/provisional sets and the resource-specific coverage mask.
- `identity_marker_coverage/`: fixed marker availability and mapping cautions.

The local `results_20261116/` contains 12 small result tables/summaries plus transfer verification. The server ZIP was **288,924 bytes**, SHA-256 `529d64f59fa621dac7197b81cd835df8a051cb1810c4258a07978a50f32d2328`; the received ZIP matched and every member passed ZIP CRC on extraction. Large row-level mapping ledgers remain on the server.

No gene-metadata operation opened, filtered, projected, normalized or changed a source count matrix. The raw preparation's separate count roundtrips must remain linked in the overall review pack. Shared-feature restriction and count-conservation checks for any future derived matrix have not been performed here.

# Stage 6 gene identity and measurement standard

Date: 2026-11-15. Scope: full-cohort preparation before biological discovery.
Status: executable project specification; harmonization has not been implemented or passed.
All thresholds, reference choices and acceptance rules below are project choices, not universal official requirements.

## Verified starting point

| Resource | Verified pilot feature axis | Source-reference information | Remaining identity gate |
|---|---|---|---|
| GSE135893 | 33,694 one-column symbols; 33,694 unique; 13 MT-prefix labels | GEO: hg38, Cell Ranger 3.0.2 | Original gene-ID crosswalk and exact annotation/reference package |
| GSE159354 | 58,302 feature IDs per library; all 16 pilot gene-ID axes match; source symbols retained | GEO: GRCh38, GENCODE basic, Cell Ranger 3.0.2 | Exact GENCODE/reference release and full-cohort feature-axis audit |
| GSE159585 | 27,714 source symbol rows; no duplicate source labels | Source counts and metadata roundtrip verified | Original gene-ID crosswalk, snRNA reference release and exon/intron counting definition |
| Original Franzén | 36,601 unique Ensembl IDs per pilot section; 10 duplicate symbol labels; all features Gene Expression | BioStudies: GRCh38, Space Ranger 1.2.2 | Exact reference package and all-section feature-axis audit |

Local evidence: [complete GSE135893 axis](../stage6_readiness/sources/GSE135893_genes.tsv.gz), [GSE135893 summary](../stage6_pilot/results/pilot_001/GSE135893/summary.json), [GSE159354 summary](../stage6_pilot/results/pilot_001/GSE159354/summary.json), [GSE159585 summary](../stage6_pilot/results/pilot_001/GSE159585/summary.json), [Franzén section example](../stage6_pilot/results/pilot_001/Franzen/V10T03-279-C1/section_summary.json), corresponding GEO family SOFT files under `../stage6_readiness/sources/`, and [BioStudies source snapshot](../stage6_readiness/sources/spatial/S-BSST1410.json).
The zero duplicate counts above describe real source axes; duplicate-label fixtures test handling and do not establish real-data duplicates.
The existing 17-marker capture audit is not stable-ID harmonization or proof of adequate cell-state/program coverage.

## Reference freeze and provenance

1. Preserve the source counting reference separately from the common identity reference. Record genome assembly, annotation provider/release, reference package, counting software/version, exon/intron policy, source URL and evidence status. Unknown fields remain unknown; feature count or software version alone cannot identify a reference release.
2. **Project common identity target: human GRCh38, Ensembl release 116 (June 2026).** This is an identity crosswalk target, not a claim that source reads were aligned to release 116 or that quantification differences have disappeared. The official [release archive](https://jun2026.archive.ensembl.org/info/website/archives/index.html) was checked on 2026-11-15.
3. **Blocker:** the release-116 annotation/crosswalk assets have not been acquired, archived or hashed. Before implementation, freeze exact files, URLs, release, acquisition date and SHA-256 in the reference manifest; verify required identity fields and source-ID history support.
4. Acquire and freeze an official HGNC complete-set snapshot plus available withdrawal/history records and Ensembl split/merge history, recording data date, acquisition date, URL and SHA-256. Those assets have not been acquired or pinned. A mutable `latest` URL or live API response without a saved snapshot does not meet this gate. See [HGNC downloads](https://www.genenames.org/download/statistics-and-files/) and [archive policy](https://www.genenames.org/download/archive/), checked 2026-11-15.

## Canonical identity and mapping rules

1. Retain immutable `resource_id`, `source_file_sha256`, `source_row_1based`, full source feature ID, source symbol, feature type, genome and source annotation evidence. Software-unique indices are storage keys, not biological identifiers.
2. Use a human Ensembl **gene** stable ID as the canonical cross-resource key, with approved HGNC symbol and HGNC ID as annotations. Preserve the full source ID, any version suffix, source symbol and mapping method. Do not substitute transcript/protein IDs for gene IDs.
3. Split a terminal `.number` version only after validating that the value is an Ensembl gene ID. Do not strip punctuation from symbols: `RP11-34P13.3` is an actual GSE135893 source label. Treat special reference-specific identifiers/suffixes explicitly rather than applying a general string cleanup.
4. Prefer source Ensembl identity or a documented source-reference gene-ID/symbol crosswalk. A symbol-only match must be unique and consistent with source-reference and frozen nomenclature evidence; current approved-symbol, previous-symbol and alias matches remain distinguishable.
5. Classify every row as `resolved_unique`, `ambiguous`, `retired_unresolved`, `unmapped` or `non_gene_feature`; record candidates, evidence and reason. Historical aliases do not default to the first result. HGNC explicitly documents that one input can match multiple genes: [symbol-checker help](https://www.genenames.org/help/multi-symbol-checker/), checked 2026-11-15.
6. A stable ID can be retired after a gene-model change. Accept historical-to-common conversion only with documented unique correspondence; hold splits, conflicting histories and unresolved merges out of primary cross-resource gene inference. See [Ensembl ID-history notes](https://www.ensembl.org/info/website/archives/index.html), checked 2026-11-15.
7. Keep unresolved rows and their counts in the source object. Do not discard whole donors because some genes are unresolved, and do not relax mapping after inspecting candidate disease effects. Assess usable gene/program coverage before assigning an analysis role.

## Duplicates and count preservation

- Different Ensembl IDs with the same symbol remain separate; the 10 duplicate Franzén symbol labels are not permission to sum counts. The [10x H5 specification](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices) distinguishes feature `id` from `name` (checked 2026-11-15).
- Unexpected repeated source feature IDs trigger review before assembly. Multiple source rows mapping to one canonical ID are withheld from the primary harmonized matrix until a documented review establishes compatible, nonduplicated count definitions. No automatic symbol-based aggregation, count duplication across one-to-many matches, or averaging of integer counts is allowed.
- Keep source integer matrices unchanged. Write any approved one-to-one projection or explicitly justified aggregation as a derived object with a mapping manifest and rule version. Normalized, integrated, imputed or ambient-adjusted values never replace the preserved source counts.
- For one-to-one projection, require exact source-to-derived values and per-cell count totals over the included source rows. Report excluded/held rows and their count mass; do not mislabel their removal as whole-matrix count conservation. Any approved aggregation must also conserve the included integer count mass and retain its row membership.
- Record reference/feature-space differences as measurement limitations. Common identity cannot correct different gene models or exon/intron quantification; [10x reference release notes](https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-reference-release-notes) document changed genes/names despite unchanged genome sequence (checked 2026-11-15).

## Mitochondrial and measurement coverage rules

- Compute technical QC on original source counts and its declared expression-feature universe, before common-feature restriction. Do not change the denominator to the cross-cohort gene intersection.
- Where source annotation is reliable, identify mitochondrial genes by the source mitochondrial chromosome annotation (`MT`/`chrM`) and reconcile source labels. Otherwise retain the pilot's case-insensitive `MT-` metric as **provisional source-prefix QC**, not a verified complete mitochondrial gene set. Report the measured mitochondrial feature list and definition per resource; zero total counts yield missing percentage.
- The mitochondrial definition does not impose a shared scRNA/snRNA/spatial exclusion threshold. The pilot's 200-feature/25% flags remain descriptive and do not become formal filters through this standard.
- Build a per-resource coverage table distinguishing `in_source_export`, `not_in_source_export`, `off_targeted_panel`, `mapping_ambiguous` and `coverage_unknown`. An absent row in an author's export is not proof that the underlying assay did not measure the gene.
- A zero count means no observed count for an available, identified feature in that unit; it does not establish biological absence. A targeted-panel omission is unassayed for that panel. Keep coverage status separate from observed values.
- Never outer-join unavailable genes as biological zeros, and never use imputed/integrated expression to supply missing ligand, receptor or recipient-response evidence. Alignment uses documented reliable shared features; inference and program coverage use each resource's eligible original measurements.

## Required outputs before gene readiness passes

1. `feature_axes`: full small feature tables for every included library/section and hashes proving matching axes where reused.
2. `reference_manifest`: separate source/common references, exact frozen Ensembl/HGNC assets and hashes, known/unknown counting definitions and evidence links.
3. `gene_mapping_ledger`: one row per source feature with preserved identity, canonical candidates, version, HGNC fields, mapping class, evidence, decision and rule version.
4. `mapping_audit`: resolved/held counts, duplicate/conflict lists, mapping coverage, included/excluded count mass and exact count-conservation checks.
5. `gene_coverage`: resource-specific measured/exported/panel availability plus shared alignment features and coverage of each proposed marker/program or claim-specific test.
6. Derived matrix/checkpoint manifest linking source checksums, mapping rule, selected rows and numerical verification; raw checkpoints remain immutable.

These are required deliverable roles, not assertions that files with these names already exist.

## Stop gates and release decision

- Do not declare gene harmonization complete while common Ensembl/HGNC assets are unpinned, source axes are incomplete, or projection/count checks have not run.
- Stop a cross-resource gene test when its identity is ambiguous, coverage is unknown/off-panel, or source count definitions cannot support the proposed comparison. Report the failed gene/test gate; do not silently recode it as zero or exclude an entire cohort without cause.
- Stop derived-matrix release on source-axis mismatch, unexplained duplicate IDs, count duplication/loss, or missing provenance. Resolve the cause and rerun the affected narrow checks.
- A cohort may pass structural mapping with some held features only if remaining coverage supports its predeclared role. Gene readiness alone does not establish cell-state recovery, anatomical inclusion, donor independence, disease-effect identifiability or formal cohort lock.
- Record implementation status and unresolved gates before Checkpoint A acceptance. This specification does not authorize molecular ranking, change P-001/P-002, or revise frozen Stage 5 evidence.

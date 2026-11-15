# Stage 6 full-cohort QC and annotation standard

Date: 2026-11-15. Status: executable specification; full-cohort implementation and final unit lock remain pending.
Scope: native discovery-side scRNA/snRNA inputs and predefined technical assessment of reserved cohorts.
This document fixes processing rules, not validated IPF numerical cutoffs, discovered states or biological findings.

## 1. Preserve measurement and provenance

- Retain immutable source files, original feature axes, integer counts, source annotations and processing versions.
- Maintain resource -> donor -> specimen/region -> actual capture/library -> cell/nucleus identifiers.
- Qualify cell barcodes by library; an identical barcode string in different captures is not a duplicate cell.
- Record assay, chemistry, counting reference, procurement, sorting fraction and original cell-calling/QC history.
- Compute QC on original measured counts before normalization, integration or candidate-gene selection.
- Keep source-axis and harmonized-gene QC distinguishable; gene aggregation can change detected-feature counts.
- Use reference mitochondrial annotations where available; record symbol-prefix fallbacks and their limitations.
- Undefined mitochondrial percentages for zero totals remain missing, never zero.

## 2. Hard failures and cell calling

- Zero-total profiles fail expression QC; retain their identifiers and exclusion reason in the audit table.
- Negative/nonfinite counts, unexplained noninteger UMI counts or broken axes quarantine the affected input for repair.
- Ambiguous barcode joins or unresolved donor attribution quarantine the affected units; do not guess an identity.
- Resolve repeated files or duplicated measurement records before analysis; do not count them twice.
- Start from author-called/annotated cells when this is the only reliable calling boundary; disclose prior selection.
- GSE135893's 105,817 additional matrix columns are neither automatically eligible cells nor verified empty droplets.
- New calling requires confirmed per-capture raw droplets, a separately recorded method and diagnostic review.
- Never infer an empty-droplet background from pooled studies or from unassigned columns alone.

## 3. Adaptive QC flags

1. Stratify by actual capture, assay (scRNA versus snRNA) and sorting/preparation condition; summarize by donor afterward.
2. For positive profiles, calculate total counts, detected source genes, mitochondrial percentage and joint QC plots.
3. Flag low log2(total counts), low log2(detected genes) and high mitochondrial percentage at 3 scaled MADs.
4. Define scaled MAD as 1.4826 * median(abs(x - median(x))); use the specified problematic direction for each metric.
5. Export each stratum's medians, MADs, actual thresholds, flag counts and retained/flagged distributions.
6. Repeat the same rule at 4 and 5 scaled MADs as bounded QC sensitivity; do not tune against candidate results.
7. Flag high log2(counts)/log2(features) at the corresponding upper 3-MAD tail for doublet review (4/5-MAD sensitivity); do not delete using an upper count ceiling alone.

Three MADs are an initial review rule, not a validated IPF exclusion threshold. No universal 200-gene/25%-MT filter applies.
Those two pilot thresholds remain descriptive inspection flags only; do not silently promote them to formal criteria.
MAD=0, small/unstable strata, multimodality or an apparently failed whole capture triggers documented review.
If borrowing reference ranges is necessary, use comparable technical/fraction strata and record the reason and provenance.
Inspect cross-library threshold outliers: within-library adaptation cannot rescue an entirely poor-quality library.
Use distinct snRNA distributions; a low typical nuclear MT fraction does not justify requiring MT=0.

## 4. Final profile disposition

| Disposition | Required interpretation and action |
|---|---|
| PASS | No unresolved QC concern, or a documented retention decision after review; endpoint adequacy is assessed separately. |
| FAIL | Hard failure, or adjudicated severe low complexity/degradation supported by joint QC evidence; retain exclusion reasons. |
| REVIEW | Discordant flags, uncertain doublet/ambient evidence, unstable thresholds or plausible unusual biology; resolve before primary use. |

For severe low complexity, review joint low-count/low-gene flags against the degraded tail, MT/background evidence and identity coherence.
The two complexity flags are correlated observations, not independent confirmation of damage.
High MT with adequate complexity, a single low-complexity flag, or coherent rare identity requires review rather than automatic deletion.
Perform adjudication at a documented profile/group level using broad identity and QC evidence, without selecting desired disease effects.
Compare retained and excluded broad-cell composition and marker coherence; disclose any disproportionate loss of epithelial/stromal groups.
Preserve rare abnormal epithelial profiles for review; rarity, stress, mixed markers or failed reference mapping alone is not a failure rule.
Preservation does not guarantee acceptance: credible technical artifacts still fail after documented review.
Freeze actual stratum thresholds and disposition rules after full-cohort diagnostics, before formal hypothesis testing.

## 5. Doublets

- Run Scrublet on raw counts separately for each actual capture, using all available called-cell backgrounds before lineage subsetting.
- Starting settings: sim_doublet_ratio=2, n_prin_comps=30, random_state=0; reduce dimensions for undersized inputs and record changes.
- expected_doublet_rate=0.05 is a software default prior, not an observed study rate; prefer documented loading/recovery information.
- If the rate is unknown, label it unknown; use 0.05 as the diagnostic software prior with 0.02/0.10 sensitivity. These project scenarios are not estimated study doublet rates or a justified exclusion quota.
- Save observed/simulated scores, threshold, expected-rate provenance and diagnostics; review failed or unseparated score distributions.
- A predicted doublet enters REVIEW. Clear multi-identity artifacts may fail after adjudication; a score alone is not the final decision.
- Review multiple incompatible identity markers, complexity and contamination together; epithelial/stromal coexpression is not itself proof.
- Do not simulate doublets across independently sorted fractions; missing singlet parents can limit detection within sorted data.
- Transcriptomic methods poorly resolve some homotypic doublets; prior author filtering also limits what can be reassessed.
- A second detector is optional for a defined unresolved issue, not a compulsory second pipeline or majority-vote truth standard.

## 6. Ambient RNA

- Keep original counts as the common primary branch; record per-capture contamination evidence and remaining limitations.
- SoupX/CellBender require suitable raw-droplet/background inputs and diagnostics; do not force their use on filtered matrices.
- GSE135893's additional columns lack confirmed per-capture background identity; they are not a ready-made SoupX input.
- GSE159354 filtered H5 and GSE159585 called-nucleus CSV do not establish the required empty-droplet background.
- DecontX without background is an optional, explicitly limited sensitivity; absent or sorted-out source populations can bias estimation.
- If suitable raw inputs become available, predefine the correction decision using technical evidence before examining candidate effects.
- Store corrected counts in a separate layer with algorithm/version/parameters; never overwrite original measurements.
- Examine marker retention and broad identity stability after correction; low-level foreign transcripts alone do not prove a new cell state.

## 7. Hierarchical annotation

| Field | Rule |
|---|---|
| source_label/source_fraction | Preserve the author's version and collection label unchanged; sorting is not verified identity or tissue abundance. |
| broad_identity | Epithelial, stromal, endothelial, immune, other or unknown; retain major backgrounds for source assessment. |
| core_identity | Resolve relevant epithelial identities and fibroblast versus pericyte/smooth-muscle identities using multiple markers. |
| state_or_program | Data-derived, provisional until discovery-side cross-donor assessment; do not predefine final abnormal states. |
| anatomical_location | Source metadata or adjudicated spatial morphology; expression alone does not establish a true-airway location. |
| proposed_lineage_origin | Separate hypothesis with evidence status; do not infer lineage from one marker or an embedding. |

Use author labels, positive/exclusion marker combinations and donor/library distributions together; retain conflicts and confidence.
The pilot's 17-marker detectability audit is not a complete identity classifier or proof of usable fine-state recovery.

### Initial identity-review dictionary

The following project seed combinations guide cluster/local-group review, not an all-positive per-cell classifier or a disease-state whitelist.
Only use mapped, measurable features. Added fine-identity markers have not yet been audited across the full project inputs.

| Identity under review | Supporting combinations | Closest ambiguity to resolve |
|---|---|---|
| Epithelial | EPCAM, KRT8, KRT18 plus a coherent epithelial program | Low EPCAM does not exclude AT1/abnormal epithelium; mesothelium also expresses keratins |
| Basal-like | KRT5, TP63; KRT15 supports | KRT17 alone does not define basal identity; preserve KRT5-low abnormal epithelial groups |
| Club-like secretory | SCGB1A1, SCGB3A1 with secretory program | Club/mucous/SCGB3A2-associated continua; expression does not prove airway anatomy |
| Goblet/mucous | MUC5AC or MUC5B with SPDEF, AGR2 | MUC5B is not goblet-exclusive; separate secretory/submucosal-gland ambiguity |
| Multiciliated | FOXJ1 with TPPP3, PIFO and ciliary machinery | Ciliating/deuterosomal groups need not fit mature multiciliated labels |
| AT1 comparator | AGER, CAV1 with epithelial context; CLIC5/RTKN2 support | PDPN/CAV1 alone are not specific; review vascular/lymphatic and stromal alternatives |
| AT2 comparator | SFTPC, SFTPA1/2 with ABCA3, LAMP3 | Low surfactant background or SFTPB alone does not establish AT2 identity |
| Fibroblast lineage | COL1A1/2, DCN, LUM; PDGFRA supports | Collagen/ACTA2 alone is insufficient; distinguish mural, epithelial and immune programs |
| Pericyte-like | RGS5, CSPG4 with MCAM, PDGFRB | Mural/fibromyocyte boundaries can remain unresolved |
| Smooth-muscle-like | MYH11, CNN1 with TAGLN, ACTA2, MYL9 | ACTA2/TAGLN alone must not turn smooth muscle into pathological fibroblasts |
| Immune | PTPRC with myeloid TYROBP/LYZ or lymphoid CD3D/E/CD79A/MS4A1 programs | LYZ is not a pan-immune identity marker |
| Vascular endothelium | PECAM1, VWF with CLDN5/KDR and endothelial program | Mixed vascular/epithelial profiles require technical and identity review |
| Lymphatic / mesothelial guards | PROX1/FLT4/CCL21/PDPN; MSLN/WT1 with KRT8/18, respectively | Avoid calling PDPN alone AT1 or every keratin-positive group airway epithelium |

For broad identity, seek at least two mutually supporting markers plus a coherent program at cluster/local-group level and inspect the closest competing identity.
This is a project manual-review minimum, not a demand that each cell has two nonzero counts or a fixed module-score cutoff.
High confidence requires convergent identity evidence and reviewed technical alternatives; provisional retains a clear broad identity with unresolved fine identity/coverage.
Unknown/unresolved retains conflicting or insufficient evidence and its reason. Failed reference mapping alone does not establish novelty or technical failure.
Record identity-marker overlap with any later response/program endpoint; test circularity rather than using the same defining genes as independent mechanistic evidence.

Reference mapping is corroborative. Record model/version, training ancestry, scores and uncertainty; preserve unknown/review assignments.
Do not use reserved donors through an atlas, fitted label-transfer model or deconvolution reference to train discovery or tune states.
Generic published normal-identity marker knowledge can be a prior; it does not permit fitting signatures/models with reserved donors.
Healthy references can lack disease states; low-confidence mapping is not equivalent to bad QC or evidence of a novel mechanism.
Maintain versioned annotation changes; inspect unintegrated data and avoid treating a donor-specific cluster as an established identity.

## 8. Donor and endpoint eligibility

- Retain donor/specimen/library dependence; sorted fractions, repeated regions and many cells do not create independent people.
- Report post-QC broad/core identity counts per donor and capture before deciding which endpoints are estimable.
- Sparse target-cell recovery limits the affected endpoint; it need not exclude an otherwise usable donor from every analysis.
- Do not require a rare state to occur in every donor, or interpret insufficient capture as biological absence.
- For primary donor-by-identity/state aggregates, use a project starting minimum of 20 retained cells/nuclei, with predefined 10/50-cell sensitivity.
- Check count depth/informative-gene coverage as well; cell number alone does not establish an adequate pseudobulk profile.
- A primary disease contrast requires at least three independent eligible donors per group and a full-rank, estimable design.
- These are project operational lower bounds, not validated universal thresholds or power guarantees; confirm adequacy before final lock.
- Do not pool studies, sorting fractions or incompatible regions to meet a threshold; record any justified within-donor aggregation.
- A cross-study-supported core state requires donor replication within both discovery studies; otherwise retain an exploratory designation.
- Sparse/rare states remain in the annotated resource even when a particular formal contrast is not estimable.
- Preserve GSE159354 sorting limitations and GSE159585's fully confounded procurement comparison; QC does not remove those restrictions.

## 9. Required outputs before final lock

Produce the full donor/specimen/capture crosswalk; per-profile QC/flags/dispositions; per-stratum numerical threshold table; doublet/ambient audit;
annotation dictionary/version/confidence log; donor-by-identity recovery table; library/donor inclusion ledger; and bounded sensitivity summary.
Record source/code/software versions, seeds, input hashes, retained/excluded counts and unresolved REVIEW units.
Finalize the eligible donor/library/profile list only after all intended captures receive QC and annotation review and endpoint coverage is assessed.
Anatomical eligibility additionally requires the separate spatial morphology/registration review; scRNA annotation cannot satisfy that gate.
The current selected pilot does not constitute a whole-cohort QC pass, final analysis-unit lock or completed Checkpoint A.

## Primary and official sources

- [OSCA adaptive QC, batch handling and cell-type loss](https://bioconductor.org/books/release/OSCA.advanced/quality-control-redux.html).
- [Scanpy 1.11 Scrublet API and defaults](https://scanpy.readthedocs.io/en/1.11.x/api/generated/scanpy.pp.scrublet.html); [Scrublet original study and missing-parent limitations](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/).
- [SoupX maintainer input/heterogeneity guidance](https://github.com/constantAmateur/SoupX); [CellBender input and diagnostic requirements](https://cellbender.readthedocs.io/en/latest/usage/index.html).
- [DecontX official background/filtered-input limitations](https://bioconductor.org/packages/release/bioc/vignettes/decontX/inst/doc/decontX.html).
- [HLCA hierarchical annotation and mapping uncertainty](https://www.nature.com/articles/s41591-023-02327-2); [SingleR annotation diagnostics](https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html).
- [HLCA original methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10287567/) support hierarchical, context-dependent marker interpretation; the seed table above is not a claimed reproduction of its unparsed Supplementary Table 6.
- [Tsukui human lung fibroblast/mural distinctions](https://pmc.ncbi.nlm.nih.gov/articles/PMC7174390/), [human lung endothelial atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC8300155/) and [human AT1 marker study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9882278/) support the ACTA2/PDPN specificity cautions.

# Stage 6 preparation rules and cohort intake lock

Version 1.0 | Specified 2026-11-15 | Decision D-029

**Rules and source-level intake are fixed for implementation. Final post-QC donor/cell/ROI eligibility and Checkpoint A have not passed.** This is a preparation specification following the completed selected technical pilot, not a claim that the full cohorts have been cleaned or that airway niches have been demonstrated.

The user requested full-cohort cleaning, gene identity, annotation and anatomical inclusion standards before the formal cohort lock. This package implements that request within D-027's two niche opportunities. It preserves discovery of new states and mechanisms; no candidate gene list, fine state or positive niche is required for entry. No new expression processing, downloads or server jobs were performed for this specification.

## Governing standards

| Standard | Fixed now | Evidence still required |
|---|---|---|
| [Cleaning and annotation](QC_ANNOTATION_STANDARD.md) | Source counts; actual-capture/assay/fraction QC; 3 scaled-MAD review flags and 4/5 sensitivity; doublet/ambient policy; identity hierarchy and endpoint coverage | Full-cohort thresholds, dispositions, annotation confidence and donor-by-identity coverage |
| [Gene identity](GENE_STANDARD.md) | Ensembl gene stable-ID key; source identity preserved; ambiguous mappings held; duplicate-symbol and unavailable-gene rules; Ensembl 116 identity target | Frozen reference assets/hashes, complete feature axes, mapping and count-conservation audits |
| [Anatomy](ANATOMY_STANDARD.md) | All-section H&E review; technical, compartment, interface and endpoint gates; basal-facing topology and geometry; expert review | Original histology or documented sufficient substitute, expert-reviewed ROI inventory and matching comparators |
| [Exact intake](cohort_intake.json) | Source-qualified donor/library/section IDs and analytical roles; untouched external reserves | Actual final, endpoint-specific eligibility after the above gates |

These are project operational rules informed by primary research and official documentation, not IPF-specific consensus cutoffs. Source counting references remain distinct from the common gene identity reference. Histological review cannot be replaced by scRNA labels or a technically aligned image.

## Locked intake and roles

| Resource | Core candidate intake | Separate background | Role and current coverage |
|---|---|---|---|
| GSE135893 | 12 IPF / 19 libraries; 10 controls / 12 libraries | 8 other-ILD donors / 12 libraries | Discovery study 1; four libraries from four donors processed in the pilot |
| GSE159354 | 4 IPF / 20 libraries; 3 controls / 15 libraries | 3 SSc-ILD donors / 13 libraries | Discovery study 2; 16 core libraries from four donors processed; clean all 48 native libraries |
| Franzén: E-MTAB-14121 / S-BSST1410 | 4 IPF / 19 sections; 4 controls / 6 sections | No additional independent acquisition | Spatial discovery; all 25 enter review, six processed; HC post-mortem versus IPF transplant/resection is procurement-confounded |
| GSE159585 | 6 IPF / 6 deposited samples; 12 controls / 12 samples | Seven COVID donors and bulk branches outside this intake | Separate snRNA robustness; four samples processed; explant/autopsy completely confounded with disease |
| GSE136831 | Reserve 32 IPF / 28 controls | 18 COPD donors separate | Primary external transcriptomic reserve; candidate outcomes not inspected |
| GSE122960 | Reserve 4 transplant IPF / 8 controls | One early-IPF cryobiopsy and four other-ILD donors separate | Supplementary external reserve; candidate outcomes not inspected |

The 134 native intake units in the manifest include background libraries and the robustness branch; they are not 134 patients. Source-qualified donor IDs do not establish cross-source patient independence. Repeated sections and libraries remain nested within their original people. Do not sum the resources into a claimed unique-patient total.

GSE227136 remains on hold because it reuses 29/30 Habermann donors; it cannot validate the entire discovery family independently. GSE128033, GSE279404 and GSE286182 remain conditional for their documented cell-calling, donor/capture, preparation and ancestry gaps. Additional spatial resources remain extensions or claim-specific reserves pending access, panel, morphology and independence checks. No conditional resource is silently promoted into the locked core.

## GSE159354: keep discovery broad and the disease comparison interpretable

All native libraries can contribute to source-aware annotation, state discovery and coverage assessment, including unsorted samples that may recover EPCAM-low states. Different gates are separate measurement routes, not extra donors.

For a **primary within-state disease-expression contrast**, use these matched ordinary gates, subject to full QC and adequate state recovery:

- Epithelial branch: ordinary EPCAM+ from F18555, F18545 and F18536 versus N1902, N1857 and N1865.
- Fibroblast branch: Triple- from the same three IPF and three control donors, with expression-confirmed fibroblast identity; the gate itself is not a fibroblast label.
- F19609's EPCAM+SAbGAL+/- and CD90+ and N1902's CD90+ remain available for discovery and fraction-specific support. Do not reconstruct an unbiased EPCAM population by weighting the SAbGAL libraries by recovered cell numbers.
- Five unsorted libraries contribute to discovery/coverage and internal expression sensitivity. Their disease effect is not an independent validation; only three IPF and two control donors have that route.

Thus 3 versus 3 is the source-level ceiling for each matched-gate primary contrast, not a guaranteed final sample size. Do not pool different routes to rescue insufficient replication or infer tissue-wide state proportions from FACS yields. If original exhaustive-gating/sampling weights become available, record a candidate-independent protocol amendment before changing this rule. The [original study](https://insight.jci.org/articles/view/143626) describes its preparation and FACS capture; the endpoint allocation above is our project decision.

## Endpoint gates and observed feasibility warning

Use original counts for donor-aware measurement; fit within-study contrasts first where sampling differs and compare effect/uncertainty across studies. A shared representation may support matching states, but does not remove capture or procurement bias. Follow the existing [integration-method review](../stage6_integration_methods/README.md): compare unintegrated, Harmony and one justified alternative, using biological conservation, donor support and technical alignment rather than candidate significance or UMAP appearance. The best representation remains an empirical discovery-side decision.

For primary donor-by-identity/state expression aggregates, fix a project starting minimum of **20 retained cells/nuclei**, with 10/50 sensitivity, sufficient count depth and identifiable measurements. Require **at least three independent eligible donors per comparison group within the study** and an estimable design. These are operational lower bounds, not a power calculation or permission to add many covariates. Cell number does not establish power. Rare states remain annotated/exploratory when a formal contrast fails this gate; their absence from controls can require a prevalence/parent-population question rather than a within-state expression test. Spatial spots/ROIs use their own assay/geometry and donor requirements, not the 20-cell rule.

The complete deposited Habermann broad labels already show a limitation: at a 20-cell screen, only **9/12 IPF and 6/10 controls** have at least 20 author-labeled mesenchymal cells, before cleaning and before distinguishing fibroblasts from other mesenchyme. The corresponding epithelial screen gives 12/12 and 10/10. These are pre-QC source-label screens, not strict bounds on reannotation or final eligibility; fine fibroblast states may have fewer usable donors. Do not claim the entire 22-person core supports every fibroblast endpoint. The manifest stores the 10/20/50 screens and original labels for audit.

For anatomy, a control containing only alveoli can support regional mapping but cannot be used as a matched airway-interface control. Absence of a target interface restricts that endpoint, not all donor measurements. Prefer a recorded within-IPF regional endpoint to an anatomically invalid disease contrast; neither substitutes for a missing independent spatial test. A regional Visium association does not prove cell contact or signaling direction.

The retained [Franzén source description](../stage6_readiness/sources/spatial/S-BSST1410.json) also states that HC tissues were collected post-mortem and IPF tissues during transplant/resection. This design cannot separate procurement from disease by covariate adjustment or integration. Prioritize within-IPF regional/interface associations for spatial discovery; HC maps provide bounded context. A morphologically matched IPF/HC comparison remains a joint disease/procurement association, not a clean disease effect. Because this acquisition participates in discovery, its HC contrast is discovery/context evidence rather than an independent validation. GSE159585 has its own fully confounded explant/autopsy restriction.

## Final lock sequence and required record

1. Acquire/process every locked native unit, preserving existing raw files and pilot checkpoints. Resolve true-capture identities and exceptional libraries; export actual per-stratum thresholds and PASS/FAIL/REVIEW reasons. GSE135893 starts from the author-annotated cell subset, not every matrix column.
2. Acquire/hash the selected reference assets; export all source feature axes; implement the mapping ledger, coverage mask and conservation audit. Unknown source releases remain disclosed; unresolved identities cannot support cross-resource single-gene claims.
3. Apply hierarchical identity annotation with uncertainty, then tabulate cell recovery by donor, original fraction and identity. Keep major background cell classes; do not require a predicted new state to be positive for inclusion.
4. Review all 25 spatial sections with registration/masks and candidate-independent morphology; obtain qualified anatomy/pathology review and an endpoint-specific ROI/comparator ledger.
5. Resolve donor ancestry, source reuse and design estimability for the proposed endpoints. Freeze actual eligible and held/excluded units, reasons, models, representation evaluation and spatial-null plans at **Checkpoint A** before formal discovery.
6. Later, after discovery iteration, freeze state/program definitions, candidate-specific measurements and predictions at **Checkpoint B** before reserved biological outcomes. Reserved technical QC may be candidate-independent; it cannot train/tune discovery through reference reuse.

The final manifest must identify `endpoint_id`, resource version, donor/specimen/capture or section/ROI, QC/annotation/anatomy decisions, gene coverage, comparator, independent-use history, inclusion/exclusion reason, reviewer/date and source/code/reference hashes. Required gates cannot be replaced by the words "pilot passed". A failed endpoint remains `hold/ineligible`, rather than being removed from the audit trail.

**Current lock level:** `rules_and_intake_locked_final_postqc_allocation_pending`. The prepared intake contains no finally eligible profiles/ROIs, no expert sign-offs and no completed Checkpoint A. This is deliberate: the full-cohort evidence needed for those decisions does not yet exist.

## Reproduction and narrow verification

Run from the research root:

```bash
python3 analysis/stage6_specification/build_and_check_intake.py
```

The script reconstructs the exact intake from preserved metadata and pilot summaries, verifies source hashes and cardinalities, and checks that no preliminary unit has been promoted to final eligibility. It writes [validation.json](validation.json). `--write` rebuilds the intake deliberately; investigate unexpected differences first. This check validates metadata accounting, not a new expression run, gene mapping or anatomical review.

## Decision record

- **Recommendation / alignment:** Implement these full-cohort preparation standards and retain two single-cell discovery studies plus native spatial discovery, within the existing airway/airway-like epithelial–fibroblast focus.
- **Evidence:** Completed [technical pilot](../stage6_pilot/RESULTS.md), [native readiness ledger](../stage6_readiness/donor_library_ledger.json), source metadata and the linked primary/official methods.
- **Key uncertainty:** Endpoint-specific cell recovery, source-reference mapping, expert anatomical comparability and independent donors after full review.
- **Alternative:** Keep a scarce state/region exploratory or change to an estimable endpoint with an explicit amendment; do not manufacture replication by pooling incompatible units.
- **Scope / user decision:** The current request authorizes specification and intake selection. No new scope decision is required. Full data preparation remains an execution step; this document does not claim a completed final biological allocation or authorize experiments/target ranking.
- **Next verifiable action:** Execute the locked intake preparation and fill the endpoint eligibility ledger; accept or decline Checkpoint A on the observed evidence.

Earlier cohort/allocation/pilot manifests remain historical inputs with their original status; this preparation specification is the current rule source under D-029. Frozen Stage 5 assets and pending P-001/P-002 remain unchanged.

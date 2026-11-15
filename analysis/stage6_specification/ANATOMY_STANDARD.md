# Stage 6 anatomy and spatial endpoint eligibility standard

Status: rules specified on 2026-11-15; full-cohort anatomical adjudication is pending.
Purpose: identify interpretable airway/airway-like epithelial–stromal interfaces without selecting regions because a molecular candidate is positive.
This standard defines review and inclusion decisions; it does not certify an existing niche or establish signaling direction or epithelial lineage origin.

## Review universe and current evidence

- Review **all 25 Franzén sections from all eight donors**: IPF_1 (6), IPF_2 (6), IPF_3 (3), IPF_4 (4), HC_1 (2), HC_2 (1), HC_3 (2), HC_4 (1).
- Use the exact section/block/replicate mapping in `../stage6_readiness/sources/spatial/hs_visium_metadata.tsv`; do not invent a balanced repeat design.
- Only six sections have completed the selected technical pilot. Archive integrity for all members does not constitute full per-section QC or anatomical review.
- The original full-resolution H&E archive remains pending. Space Ranger low/hires images support preliminary review; they are downsampled, not original-resolution histology.
- Use original full-resolution H&E for adjudication unless a qualified reviewer explicitly documents that the available image resolves the particular compartment/interface sufficiently. Record resolution, limitations and eligible endpoints.
- Published fibrosis scores remain source annotations. Blocks/sections are cross-sectional observations, not a longitudinal progression series.
- The native S-BSST1410 description specifies HC post-mortem versus IPF transplant/resection tissue. Procurement and disease are inseparable in this acquisition; matched anatomy, QC or integration cannot identify a disease-only effect. Prioritize within-IPF regional associations; treat HC contrasts as bounded disease/procurement context.

## Separate review gates

| Gate | Required evidence | Decision |
|---|---|---|
| Technical | Exact count/barcode/coordinate join; recorded image transform and scale; credible overlay; usable tissue and artifact masks | `pass`, `hold`, `fail` |
| Anatomical | Sufficient image detail; compartment boundaries and uncertainty recorded; reviewer identity and evidence | `expert_confirmed`, `expert_uncertain`, `provisional`, `unreviewed` |
| Interface | A recognizable epithelial boundary and adjacent abluminal/basal-side stroma, with continuity and barriers assessed at assay resolution | `present`, `not_observed`, `uncertain`, `not_assessable`; `not_assessed` until review |
| Endpoint | Required anatomy, measurement, lineage coverage and comparator are available for the specified question | `eligible`, `limited`, `hold`, `ineligible` |

A failed join, unavailable image or unresolved registration places the affected section/region on hold until corrected; retain the original files and reason.
The source `in_tissue` flag indicates tissue coverage, not airway identity, cellular purity or anatomical approval.
An anatomically absent interface and an interface that cannot be assessed are different findings. Neither automatically excludes other usable endpoints from the donor.

## Candidate-independent review procedure

1. Review H&E and technical overlays before opening candidate-expression maps or inferred niche labels; record any prior molecular exposure.
2. Inventory all target-compatible anatomical interfaces and prespecified comparison compartments across every section, including negative and uncertain regions.
3. Store region polygons and artifact/barrier masks in the original image coordinate frame, with every transformation explicitly recorded.
4. Assign compartment, interface and confidence labels from morphology; retain the author's annotations separately from project annotations.
5. Obtain qualified anatomical/pathology review. Codex-assisted or nonexpert annotations remain `provisional` until reviewed; do not fabricate an expert sign-off.
6. Resolve disagreements or retain `expert_uncertain`, with an adjudication note. Freeze the candidate-independent region inventory before molecular niche selection.
7. Overlay measured RNA and lineage assignments to assess assay-specific endpoint eligibility. Do not redraw regions to maximize candidate enrichment.

Record whether diagnosis was available to reviewers; do not claim diagnosis blinding when morphology or source labels disclosed it.
Any later region revision requires its original geometry, reason, reviewer and date; preserve the prior version.

## Anatomical labels and topology

Use separate compartment and interface fields; a region may contain several recorded compartments.

| Compartment label | Inclusion meaning and boundary |
|---|---|
| `conducting_airway` / `small_airway` | Morphologically supported airway; use the more specific label only when image evidence permits |
| `airway_like_remodeled` | Bronchiolized/cystic/honeycomb-associated lining when supported; record the specific morphology and uncertainty |
| `alveolar_preserved` / `alveolar_remodeled` | Alveolar regions, kept distinct from airway interfaces |
| `fibroblastic_focus` / `dense_fibrosis` | Stromal pathology labels requiring review; neither alone establishes a particular fibroblast state |
| `vascular_perivascular` / `pleural_subpleural` | Context compartments; not substituted for the target epithelial–stromal interface |
| `lumen_acellular` / `artifact` / `indeterminate` | Exclusion/barrier or uncertainty annotations, with explicit reasons |

- Distinguish true abluminal stroma from cells across a lumen. Short two-dimensional distance across a cavity does not establish a tissue interface.
- Mark lumens, tears, folds, detached tissue, image defects and tissue edges. Do not connect disconnected tissue islands or bridge excluded gaps in a neighborhood graph.
- Record mixed regions and boundary uncertainty rather than forcing every spot/cell into a precise compartment.
- Morphology, epithelial identity, expression state and developmental origin are separate variables. Airway-like morphology or KRT5 expression does not prove airway-derived lineage.
- Visium measurements support regional mixtures/localization at their resolution; co-occurrence in a spot does not demonstrate direct cell contact.

## Endpoint-specific inclusion

| Endpoint | Required evidence | Permitted limitation |
|---|---|---|
| Regional state/program localization | Technical pass, interpretable region, relevant measured genes and lineage evidence | No direct cell-contact or signaling-direction claim |
| Airway/airway-like interface association | Expert-supported interface and measurable epithelial/stromal components; uncertainty and tissue geometry retained | Visium inference remains regional; insufficient interface resolution means hold |
| Disease–control interface comparison | The same interpretable anatomical comparison is represented in both groups, with donor-level estimability | An alveolar-only control is not matched; Franzén's procurement confounding prevents a disease-only interpretation |
| Within-IPF regional comparison | Appropriate distinct regions within IPF, with block/section dependence represented | Does not estimate a healthy comparison or temporal transition |
| Cell-neighborhood/contact-oriented test | Cell-resolved measurements, usable segmentation and panel, physically valid topology | HD bins are not automatically cells; proximity remains noncausal |
| Regional ROI program test | Valid ROI anatomy and measured RNA, with appropriate donor/comparator structure | Disconnected GeoMx ROIs do not establish continuous cell adjacency |

Apply broad-lineage coverage requirements from the cell-annotation/QC specification; do not require the candidate fine state or signal to be positive for inclusion.
Record measurable absence separately from insufficient lineage recovery, uncertain mapping and unmeasured panel genes.
Do not impose a universal minimum ROI/cell count as proof of statistical power. Record usable donors, interface extent, regions and uncertainty for each estimable endpoint.
Retain donors/regions for endpoints they can support; a failed fibroblast or airway-interface endpoint does not automatically delete valid epithelial or alveolar measurements.
Anatomical comparability does not require removing the disease-defining pathology; specify the comparison and covariates for the biological question.

## Required review record

Each record must contain:

- Identity: `resource_version`, `donor_id`, `specimen_id`, `block_id`, `section_id`, `replicate_id`, `roi_id`, `diagnosis_source`.
- Images: `image_path`, `image_sha256`, `image_type`, `image_dimensions`, `pixel_size_or_scale`, `coordinate_frame`, `transform`, `overlay_path`.
- Geometry: `roi_polygon`, `artifact_mask`, `barrier_mask`, `tissue_component_id`, `interface_geometry`, `interface_extent`, `boundary_uncertainty`.
- Morphology: `source_annotation`, `project_compartment`, `interface_type`, `interface_status`, `compartment_confidence`, `sampling_selection`.
- Review: `reviewer`, `reviewer_role`, `review_date`, `review_status`, `image_sufficiency`, `molecular_exposure`, `disagreement_or_adjudication`.
- Decisions: `technical_gate`, `endpoint_id`, `endpoint_status`, `eligible_comparator`, `lineage_coverage`, `reason`, `record_version`.

Use explicit `unknown`/`not_available` values for unresolved metadata; do not replace a missing block or specimen identity with a guessed identifier.
An endpoint cannot be finally `eligible` while a required technical/anatomical gate is unresolved. Preliminary records remain usable as an auditable review queue.

## Replication, geometry and final lock

Preserve donor → specimen/block → section → ROI → spot/cell nesting. Repeated sections and ROIs do not add independent people.
Build physical neighborhoods only within the same valid section/tissue component, retaining background populations and barrier information.
Choose physical scales appropriate to the actual assay and geometry before candidate-effect testing; record units, calibration and sensitivity scales.
Spatial nulls/comparisons must respect donor, compartment, tissue geometry and relevant abundance structure. Unrestricted label shuffling or a donor intercept alone is not sufficient spatial control.
Freeze the final endpoint-specific inclusion ledger only after all intended sections are reviewed and required gates resolved; include excluded/held records and reasons.
Until then the eight-donor/25-section list is the **candidate review universe**, not an unconditional formal airway-niche analysis cohort.

## Primary and official sources

- [Original Franzén spatial study](https://doi.org/10.1038/s41588-024-01819-2) and [original acquisition, S-BSST1410](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410); donor/section counts above use the retained native metadata.
- [10x Space Ranger spatial outputs](https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs): downsampled images, original-coordinate positions, scaling, tissue flags and registration review. Documentation version explains format semantics; it does not assert the source's software version.
- [10x description of Visium capture](https://www.10xgenomics.com/blog/answering-your-questions-about-the-visium-spatial-gene-expression-solution): spots can contain multiple cells. The endpoint eligibility and anatomy rubric above are project rules, not vendor-certified pathology thresholds.

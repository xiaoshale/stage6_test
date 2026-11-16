# Human review order and decisions

Open the offline [review page](results/review_pack/index.html). Work from measurement quality to identity, anatomy and endpoint eligibility. Save reviewed CSV copies separately; do not overwrite the generated evidence or mark an entire resource accepted from one representative plot.

## RNA: one capture at a time

Use `human_decisions/rna_capture_decisions.csv` in the package. Inspect the joint QC distributions, 3/4/5-MAD sensitivity, original labels, positive and competing identity markers, and simulated/observed doublet distributions. An outlier flag is a review queue, not an exclusion decision. Do not impose the software's doublet prior as a removal quota.

For each disposition, record reviewer/date, affected source-qualified profiles or identity group, reason and follow-up. If a group requires further examination, record a hold. Source-symbol marker scores are provisional, and CLIC5/RGS5 mapping cautions must be checked against the gene reports. A score is not an expert identity, lineage origin or disease state.

Detailed `obs_review.tsv.gz` tables remain under the server run's `rna_review/<resource>/<capture>/` folders. Use these original cell IDs when creating later reviewed masks; the portable summary does not silently exclude any profile.

## Tissue: inspect every section before candidate expression

Use `human_decisions/spatial_section_anatomy.csv`, the 25-section contact sheet and each section's original-image overview/fixed grid. Confirm compartment, basal-facing epithelium–stroma interface and image registration. Document lumens, folds, tears, edges and disconnected islands that prevent a valid neighborhood interpretation.

Inventory eligible, uncertain and unsuitable regions in `human_decisions/roi_inventory.csv`; do not choose fields because a favored candidate is expressed there. Full-resolution source H&E and JSON mask templates remain on the server. Record region geometry in the stated image coordinate system, with reviewer/date and rationale. Empty inventories/masks mean review has not occurred.

## Lock only the endpoint that the evidence supports

After QC/identity and anatomy adjudication, recalculate eligible cells/regions and independent donors per proposed endpoint. Preserve FACS fraction, assay and procurement restrictions; library/section/ROI replicates do not create extra donors. The 20-profile floor with 10/50 sensitivity is an operational starting point, not a power guarantee or global donor filter.

Resolve relevant source-reference and cross-resource donor-independence questions, then assess the design and record the final endpoint-specific allocation. This review package alone does not pass Checkpoint A or justify starting niche/target ranking.

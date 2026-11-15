# Integrating public single-cell, single-nucleus and spatial data for Stage 6

Review date: 2026-11-14. Status: methods research and a proposed project framework. No IPF integration, algorithm benchmark, matrix download, installation or server analysis was performed for this review. This is a targeted review of original research and maintainer documentation, not an exhaustive systematic review of all integration methods.

## Recommendation

Use three connected analytical components:

1. **A discovery-side shared cell-state representation**, built after source-specific quality control and metadata harmonization.
2. **Donor-aware expression and abundance analyses**, retaining original measurements and study/region provenance.
3. **Spatial evidence evaluated in its native measurement units**, linked through cell states and measured expression, with independent donors reserved for fixed prediction tests.

These components answer different questions. A common embedding identifies comparable cells; it does not itself estimate a disease effect or demonstrate a physical niche. Integration can also mean synthesizing comparable study-level effects without placing all measurements in one matrix.

The objective remains a small number of airway/airway-like epithelial–fibroblast niches with cross-donor support, defined anatomy and a specific intervention question. Preserve major cellular backgrounds for source attribution and neighborhood interpretation. Stage 5 supplies biological anchors and counterevidence, not a compulsory molecular shortlist.

## What the sources support

| Source | Direct methodological lesson | Applicability limit |
|---|---|---|
| [Luecken et al., Nature Methods 2022](https://doi.org/10.1038/s41592-021-01336-8) | Evaluate biological conservation alongside batch removal; output type and preprocessing affect comparisons. | A historical benchmark ranking does not select the best current IPF workflow. |
| [Antonsson and Melsted, Genome Research 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12315870/) | Pseudobatch experiments show that correction can alter data even when little correction is needed. Include an unintegrated baseline. | Its preference for Harmony concerns the tested calibration settings and implementations, not every rare-state or disease task. |
| [Sikkema et al., HLCA, Nature Medicine 2023](https://doi.org/10.1038/s41591-023-02327-2) | A lung-specific example combines metadata and hierarchical labels, benchmarks representations, and uses uncertain reference mapping. | The paper reports that some rare populations were not separated by Harmony/RPCA in its tests. An extended disease atlas may contain our reserved donors. |
| [Seurat official integration tutorial](https://satijalab.org/seurat/articles/integration.html) | Current integration workflows retain original layers and offer alternative shared embeddings. | A software tutorial is not a donor-aware disease-inference protocol. |
| [scvi-tools official introduction](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/quick_start/api_overview.html) and [SCVI API](https://docs.scvi-tools.org/en/latest/api/reference/scvi.model.SCVI.html#scvi.model.SCVI.setup_anndata) | Preserve counts as model input; document batch and nuisance covariates. | Covariate registration cannot identify disease separately from a completely confounded study effect. |
| [Squair et al., Nature Communications 2021](https://doi.org/10.1038/s41467-021-25960-2) | Between-replicate variability matters for condition-level differential expression. | Aggregation does not manufacture missing donors or solve a poor clinical design. |
| [Nguyen et al., Nature Communications 2023](https://doi.org/10.1038/s41467-023-37126-3) | Multi-batch DE performance depends on depth, sparsity, batch effects and the evaluation endpoint. | Disease-gene ranking performance and calibrated inference across patients are different targets; see the discussion below. |
| [cell2location original method](https://doi.org/10.1038/s41587-021-01139-4) and [official tutorial](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html) | Link reference expression signatures to measured spatial counts while modeling technical variation. | Estimated cell abundance is not an observed cell boundary or direct signaling. |
| [Li et al., spatial-deconvolution benchmark, Nature Communications 2023](https://doi.org/10.1038/s41467-023-37168-7) | Reference and measurement characteristics affect spatial deconvolution performance. | Benchmark winners are not universal IPF winners. |
| [Franzén et al., Nature Genetics 2024](https://doi.org/10.1038/s41588-024-01819-2) | A relevant pulmonary-fibrosis example separates spatial representation/mapping from raw-count donor-level disease analysis. | Its exact regressions, labels and parameters are not automatically appropriate here. |

`sources.json` contains the fuller source register, official-document links and bounded claims. Maintainer documentation establishes software input/output contracts; original research and benchmarks provide empirical evidence. Their combined use does not create a universal official IPF analysis standard.

## 1. Establish what can be compared

Start from original usable count matrices and provenance. FASTQ reprocessing is warranted when essential counting differences cannot otherwise be reconciled, rather than being a prerequisite for every public dataset.

- Maintain the donor–specimen–region–library–assay hierarchy. Resolve reused donors and specimens before adding sample totals.
- Harmonize gene identifiers, diagnosis, tissue region, treatment, sorting, technology, counting reference and original/reused status. Preserve unknown metadata explicitly.
- Determine whether an object contains counts, normalized values, residuals, imputed values or a learned representation. Do not treat these as interchangeable inputs.
- Apply quality-control principles by library and assay, with documented thresholds adapted to their distributions. A single mitochondrial-percentage or gene-count cutoff is not appropriate merely because it is convenient.
- Preserve source count layers and sample identities throughout processing. Align genes for each analysis without encoding an unmeasured targeted-panel gene as a biological zero.
- Examine whether disease, anatomical location, preparation and study can be separated statistically. A study containing only IPF and another containing only controls cannot provide a clean study-adjusted disease contrast simply through batch correction.

scRNA-seq and snRNA-seq can share an annotation framework, but nuclear/cytoplasmic capture and cell recovery differ. [Denisenko et al.](https://doi.org/10.1186/s13059-020-02048-6) and [Slyper et al.](https://doi.org/10.1038/s41591-020-0844-1) provide empirical preparation comparisons in other tissues. These support explicit modality/sampling checks, not importing their effect sizes into IPF. Analyze within-modality effects first where possible and assess whether conclusions survive stratification.

## 2. Build a shared discovery representation without forcing biological identity

Use broad lineage identification before focused epithelial and mesenchymal state analysis, retaining relevant backgrounds. Compare original labels using multiple markers, programs and anatomical information. Labels such as "basal," "basaloid" or "activated fibroblast" are not interchangeable merely because publications use similar words. Keep uncertain or genuinely unshared populations visible.

Start with a small comparison set: **unintegrated baseline, Harmony, and one justified alternative such as RPCA or scVI**. Match implementations to the working environment and input assumptions. This is a proposed economical comparison, not an empirical ranking already established on our data. These are alternative representations; stacking them sequentially is not the default plan.

Assess a concise set of properties:

- Do shared lineages align across studies where comparable populations exist?
- Are fine epithelial states, fibroblast programs and rare/unshared populations preserved?
- Are results supported by multiple donors, rather than one library or processing condition?
- Do state assignments and biological interpretations remain recognizable without integration and under a reasonable alternative?
- What are the computational cost and reproducibility of the representation?

Use [scIB's separate biological-conservation and batch-removal metrics](https://scib.readthedocs.io/en/latest/api/scib.metrics.metrics_all.html) together with marker/program and donor checks. Coarse labels cannot establish preservation of an unknown disease state, and a single combined score or attractive UMAP is insufficient. Do not choose parameters by maximizing a candidate's disease significance or its external replication.

The apparently different conclusions of the calibration study and the HLCA illustrate why this matters: avoiding artificial change and separating rare populations are related but distinct tasks. Disease or a relevant anatomical compartment must not be automatically designated as a nuisance effect to remove. Any dataset-specific population remains a biological or technical question to investigate, rather than an automatic integration failure.

## 3. Estimate disease effects using donors and suitable measurements

For the main condition-level questions, the proposed starting approach is to aggregate original counts within donor and cell type/state, retaining region when it defines the question. Fit an appropriate count or precision-weighted expression model with clinical/technical covariates that the design can actually estimate. Repeated regions from the same person require paired/hierarchical handling or an explicitly justified aggregation.

[Zimmerman et al.](https://doi.org/10.1038/s41467-021-21038-1) supports explicit correlation modeling and finds aggregation conservative in its settings. [Hafner et al. 2025](https://doi.org/10.1093/bib/bbaf397) evaluates nested single-study and atlas scenarios, including calibration and runtime. These reinforce a donor-aware choice rather than a universal pseudobulk rule. [The official dreamlet tutorial](https://diseaseneurogenomics.github.io/dreamlet/articles/dreamlet.html) describes sample-by-cell-type count aggregation and precision-weighted mixed models; this is an option for repeated measurements, not an assurance that a small dataset can estimate many covariates.

Two complementary routes are available:

| Route | When reasonable | Required checks |
|---|---|---|
| Joint donor-level model including study | Diagnoses, sampling and measurements are sufficiently comparable, with estimable study/disease contrasts | Design rank, confounding, residual variation, influential donors, source-specific direction and sensitivity |
| Within-study contrasts followed by effect synthesis | Tissue, chemistry or selection differs substantially, but the same biological contrast can be estimated within studies | Common state/endpoint definition, effect direction and uncertainty, heterogeneity and donor independence; formal meta-analysis only when sufficiently supported |

Summaries should report effect sizes and uncertainty, with multiple-testing handling where applicable. An intersection of significant-gene lists discards information about power and effect size. A large number of cells does not create a large number of patients. Study-level synthesis likewise does not make overlapping patients independent.

There is a genuine literature qualification. Nguyen et al. benchmarked multi-batch DE workflows, reported limited pseudobulk sensitivity in some settings, and favored cell-level approaches for disease-gene prioritization. That is a reason to assess sensitivity and model specification, not evidence that cells from one lung are independent patient replicates. For the condition-comparison branch, our target is a reproducible between-donor disease effect; donor-aware inference is therefore the primary proposal. Other model-specific DE procedures may be valid with their own assumptions and calibration. We are not asserting that all corrected-expression DE is prohibited or that pseudobulk is always most powerful.

Other core niche questions concern regional variation, state–response association or neighborhood enrichment within IPF. Define their own estimands and comparisons while retaining donor dependence; the project is not restricted to an IPF-versus-control differential-expression screen.

Separate expression from abundance. A shared state can support a within-state expression contrast; an IPF-restricted state absent from controls cannot automatically support the same within-state case–control test. Its prevalence or abundance and its relationship to a broader parent population are different questions. GSE159354's sorted fractions also require explicit sampling limitations: their observed proportions are not absolute tissue abundance, and gating can change which states are captured.

[muscat's original multi-sample study](https://doi.org/10.1038/s41467-020-19894-4) is useful for distinguishing cell-population definition, within-population state expression and abundance questions. For our niche hypotheses, report these separately instead of treating a change in any one as proof of the others.

## 4. Connect single-cell and spatial measurements through states and observed expression

For Visium-like data, use discovery-side cell states with original reference counts to estimate reference expression signatures, then map them to each spatial library's measured counts using a justified method such as cell2location. Its tutorial uses untransformed spatial counts and accommodates technical variation in reference-signature estimation. A Harmony/scVI coordinate is not a gene-expression signature for spatial deconvolution.

Preserve each section/core's physical coordinates and pathology. A joint transcriptomic embedding does not create physical neighbors across disconnected sections, cores or patients. Review mapped states against registered morphology and measured markers, and examine uncertainty. Regional programs and within-section associations should then be summarized or modeled with donor structure preserved.

For Xenium/CosMx, first establish the actual measured panel and usable segmentation. Evaluate shared cell states and directly measured outputs; an absent panel gene is untestable. A common measured marker/program subset may enable a specific cross-platform comparison, while broad spatial data retain the wider discovery space. Imputed off-panel expression is a prediction, not independent measured confirmation.

[Xenium output documentation](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs) distinguishes transcript locations, morphology and segmentation masks/polygons. [Space Ranger documentation](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/segmented-outputs) distinguishes binned and segmented outputs. Current software capabilities do not prove that an older public archive contains reliable cell-level segmentation. GeoMx isolated ROIs support regional expression comparisons, not contiguous cell-neighborhood reconstruction.

For the core niche question, separate three endpoints: state/region presence, local state-pair enrichment, and recipient response. Account for anatomy, local cell abundance and donor. Use suitable spatial null comparisons within the actual tissue geometry. Avoid circular tests in which the same genes both define a cell/state or deconvolution estimate and constitute the claimed response. These are analysis-design requirements to specify before formal testing; no numerical spatial threshold has been selected in this review.

## 5. Protect external testing and retain the scientific objective

Under the current provisional allocation, discovery-side alignment would use GSE135893 and native GSE159354, with GSE128033 if informative after provenance checks. Franzén/Mayr broad spatial measurements and panel-covered Vannan neighborhoods can refine discovery. These roles can change after bounded resource screening and readiness checks, before reserved outcomes are used.

Keep GSE136831 and appropriate independent spatial donors out of reference construction, candidate selection and parameter tuning when claiming their later external validation. Mayr's integrated PF-ILD reference includes reserved Adams/Reyfman data, so use original spatial counts with a discovery-only reference rather than that integrated single-cell object. An exact HLCA model/object version also needs a constituent-donor audit; the healthy core and extended disease atlas have different provenance.

After discovery, fix the state/program definitions, region, candidate, direction and endpoint before applying the predefined mapping/testing procedure to reserved data. Report negative, inconclusive and unmeasurable predictions distinctly. A study that already contributed to discovery cannot regain external-validation status after a role switch.

Broad abnormal-epithelial/myofibroblast localization is established prior work, including [Mayr et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC11313858/). The intended contribution is a more specific state–response–anatomy–regulation link with an intervention question. Spatial RNA association alone does not establish signaling direction, necessity or repair-support failure; these remain later functional questions.

## Next verifiable work

Complete the bounded resource gap check and the donor/library/panel/histology ledger. Then use discovery-side native IPF data for complete-file acquisition, parsing and technical QC; document counts and annotation versions, parameter provenance and computational needs. Assess whether cross-donor target coverage supports the proposed comparisons before selecting the main representation and locking the formal analysis specification.

This recommendation does not change the approved biological scope, select a final molecule, or claim that any integration method has passed on our IPF data. Algorithm choice and numerical criteria remain empirical decisions for a later, authorized data-analysis step.

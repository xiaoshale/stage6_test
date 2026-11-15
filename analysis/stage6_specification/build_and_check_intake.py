"""Build/check the metadata intake snapshot; this does not clean expression data."""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
INPUTS = [
    "analysis/stage6_readiness/donor_library_ledger.json",
    "analysis/stage6_readiness/sources/GSE159585_snRNAseq_metadata.csv.gz",
    "analysis/stage6_readiness/sources/spatial/hs_visium_metadata.tsv",
    "analysis/stage6_readiness/sources/spatial/S-BSST1410.json",
] + [
    f"analysis/stage6_pilot/results/pilot_001/{name}/summary.json"
    for name in ("GSE135893", "GSE159354", "GSE159585", "Franzen")
]


def read_json(path):
    return json.loads(path.read_text())


def build(root):
    ledger = read_json(root / INPUTS[0])["cohorts"]
    pilots = {
        name: read_json(root / f"analysis/stage6_pilot/results/pilot_001/{name}/summary.json")
        for name in ("GSE135893", "GSE159354", "GSE159585", "Franzen")
    }
    resources = {}

    def unit(resource, donor, library, diagnosis, group, pilot, **details):
        return {
            "unit_uid": f"{resource}::{library}",
            "donor_uid": f"{resource}::{donor}",
            "donor_id": donor,
            "library_or_section_id": library,
            "source_diagnosis": diagnosis,
            "analysis_group": group,
            "technical_processing": "pilot_pass" if pilot else "not_processed",
            "full_cohort_qc_disposition": "not_adjudicated",
            "project_annotation_status": "not_reviewed",
            "endpoint_eligibility": "not_adjudicated",
            **details,
        }

    def add(resource, role, modality, units, limits):
        counts = {}
        for group in sorted({u["analysis_group"] for u in units}):
            subset = [u for u in units if u["analysis_group"] == group]
            counts[group] = {
                "source_donors": len({u["donor_uid"] for u in subset}),
                "units": len(subset),
                "pilot_processed_units": sum(u["technical_processing"] == "pilot_pass" for u in subset),
            }
        resources[resource] = {
            "role": role,
            "modality": modality,
            "intake_status": "scope_locked_postqc_pending",
            "cross_source_independence": "not_fully_resolved",
            "counts": counts,
            "limits": limits,
            "units": sorted(units, key=lambda u: u["unit_uid"]),
        }

    units = []
    for row in ledger["GSE135893"]["library_rows"]:
        diagnosis = row["source_diagnosis"]
        units.append(unit(
            "GSE135893", row["donor_id"], row["library_id"], diagnosis,
            diagnosis if diagnosis in ("IPF", "Control") else "other_ILD",
            row["donor_id"] in pilots["GSE135893"]["selected_donors"],
            specimen_id=row["specimen_id"], region_id=row["region_id"],
            facs_fraction=row["facs_fraction"],
            capture_identity="orig.ident library proxy; specimen/capture crosswalk unresolved",
            gsm_crosswalk_status=row["geo_crosswalk_status"],
            source_broad_label_counts=row["provided_broad_lineage_cell_counts"],
        ))
    add("GSE135893", "primary_discovery_study_1", "scRNA", units, [
        "Only author-annotated called-cell subset is eligible; extra matrix columns are not certified empty droplets or cells.",
        "Other-ILD diagnosis discrepancy and missing exact region/fraction mapping remain documented.",
        "Source mesenchymal counts are not verified fibroblast counts.",
    ])

    pilot_gsms = {r["gsm"] for r in pilots["GSE159354"]["libraries"]}
    units = []
    for row in ledger["GSE159354"]["library_rows"]:
        diagnosis = row["diagnosis"]
        fraction = row["facs_fraction"]
        units.append(unit(
            "GSE159354", row["donor_id"], row["library_id"], diagnosis,
            diagnosis if diagnosis in ("IPF", "Control") else "other_ILD",
            row["gsm"] in pilot_gsms,
            specimen_id=row["specimen_id"], region_id=row["region_id"],
            facs_fraction=fraction, submitter_library_id=row["submitter_library_id"],
            primary_comparable_route=(
                "none" if diagnosis not in ("IPF", "Control") else
                "epithelial" if fraction == "EPCAM+" else
                "stromal" if fraction == "Triple-" else "none"
            ),
            capture_identity="separate submitted native H5 library; retain capture provenance",
            discovery_state_search_eligible=True,
            unsorted_disease_effect_sensitivity_only=fraction == "Unsorted",
        ))
    add("GSE159354", "primary_discovery_study_2", "scRNA", units, [
        "Clean all 48 native libraries; core IPF/control intake has 35 libraries including five unsorted.",
        "Fraction-specific expression; no absolute tissue proportions from gated libraries.",
        "SAbGAL-split/CD90 routes and unsorted libraries are not exchangeable with ordinary EPCAM+/Triple- routes.",
    ])

    with gzip.open(root / INPUTS[1], "rt") as handle:
        nuclei = defaultdict(list)
        for row in csv.DictReader(handle, delimiter=";"):
            if row["condition"] in ("IPF", "control"):
                nuclei[(row["patient_ID"], row["sample_ID"], row["condition"])].append(row["cellType"])
    units = []
    for (donor, library, diagnosis), labels in nuclei.items():
        units.append(unit(
            "GSE159585", donor, library, diagnosis,
            "Control" if diagnosis == "control" else "IPF",
            donor in pilots["GSE159585"]["donors"],
            specimen_id=None, region_id=None, facs_fraction=None,
            procurement="autopsy" if diagnosis == "control" else "explant",
            capture_identity="deposited sample_ID; verify physical capture before doublet modeling",
            source_broad_label_counts=dict(sorted(Counter(labels).items())),
        ))
    add("GSE159585", "separate_procurement_confounded_robustness", "snRNA", units, [
        "Disease and procurement are completely confounded; not a clean third discovery or independent disease-validation cohort.",
        "Seven COVID donors and bulk branches excluded from this core intake.",
    ])

    pilot_sections = {r["section_id"] for r in pilots["Franzen"]["sections"]}
    with (root / INPUTS[2]).open() as handle:
        sections = list(csv.DictReader(handle, delimiter="\t"))
    units = []
    for row in sections:
        units.append(unit(
            "Franzen", row["subject_alias"], row["sample_id"], row["condition"],
            "Control" if row["condition"] == "control" else "IPF",
            row["sample_id"] in pilot_sections,
            source_accessions=["E-MTAB-14121", "S-BSST1410"],
            source_tissue_alias=row["tissue_alias"],
            source_region_label=row["B_tissue_selection"],
            source_replicate_label=row["replicate"],
            source_fibrosis_score=row["fibrotic_extent_score_by_pathologist_0-3"],
            source_procurement="post_mortem" if row["condition"] == "control" else "transplant_or_resection",
            procurement_evidence="S-BSST1410 description; disease/procurement inseparable in this acquisition",
            expert_anatomy_status="unreviewed", interface_status="not_assessed",
            roi_review_records=[],
        ))
    add("Franzen", "primary_spatial_discovery", "Visium", units, [
        "One acquisition represented by two accessions; eight donors, 25 sections, unbalanced repeats.",
        "Original full-resolution H&E and expert review pending; technical images do not certify airway interfaces.",
        "HC post-mortem versus IPF transplant/resection procurement is confounded with disease; matched anatomy does not identify a disease-only effect.",
        "B1/B2/B3 and repeated sections are not independent people or longitudinal stages.",
    ])

    external = {}
    adams = ledger["GSE136831"]
    external["GSE136831"] = {
        "role": "reserved_primary_transcriptomic_validation",
        "candidate_outcomes_inspected": False,
        "status": "reserved_not_finally_eligible",
        "donors_by_group": {
            group: sorted(k for k, v in adams["donors"].items() if v["diagnoses"] == [group])
            for group in ("IPF", "Control", "COPD")
        },
        "library_ids_by_group": {
            group: sorted(k for k, v in adams["libraries"].items() if v["diagnoses"] == [group])
            for group in ("IPF", "Control", "COPD")
        },
        "independence_status": "source ancestry partly audited; cross-source donor independence not fully resolved",
        "limits": ["Eleven control donors reused from Raredon; preserve ancestry, not new replication."],
    }
    reyfman = defaultdict(list)
    for row in ledger["GSE122960"]["library_rows"]:
        group = (
            "early_IPF_cryobiopsy" if "early" in row["specimen_branch"] else
            "Control" if row["diagnosis"] == "Donor" else
            "IPF" if row["diagnosis"] == "Idiopathic pulmonary fibrosis" else "other_ILD"
        )
        reyfman[group].append({
            "donor_label": row["donor_label"], "gsm": row["gsm"],
            "source_diagnosis": row["diagnosis"], "specimen_branch": row["specimen_branch"],
        })
    external["GSE122960"] = {
        "role": "reserved_supplementary_transcriptomic_validation",
        "candidate_outcomes_inspected": False,
        "status": "reserved_not_finally_eligible",
        "donors_by_branch": dict(sorted(reyfman.items())),
        "independence_status": "cross-source donor identity remains to be verified",
        "limits": ["Single early-IPF cryobiopsy is separate and cannot establish a disease-stage trajectory."],
    }

    coverage = {}
    for resource, labels in (("GSE135893", ("Epithelial", "Mesenchymal")), ("GSE159585", ("epithelial", "stroma"))):
        coverage[resource] = {}
        for label in labels:
            donor_counts = Counter()
            for u in resources[resource]["units"]:
                donor_counts[(u["analysis_group"], u["donor_id"])] += u["source_broad_label_counts"].get(label, 0)
            coverage[resource][label] = {
                group: {
                    str(threshold): sum(n >= threshold for (g, _), n in donor_counts.items() if g == group)
                    for threshold in (10, 20, 50)
                }
                for group in ("IPF", "Control")
            }

    return {
        "version": "1.0", "specified_on": "2026-11-15",
        "decision": "D-029",
        "status": "rules_and_intake_locked_final_postqc_allocation_pending",
        "checkpoint_A_passed": False,
        "formal_discovery_started": False,
        "independence_boundary": "Counts are source-qualified donors; no pooled unique-patient total is asserted.",
        "source_artifacts": [
            {"path": p, "sha256": hashlib.sha256((root / p).read_bytes()).hexdigest()}
            for p in INPUTS
        ],
        "resources": resources,
        "source_label_coverage_screen": {
            "boundary": "Pre-QC author broad labels, aggregated by donor; not verified fibroblast/fine-state eligibility or evidence of absence.",
            "donors_at_least_n_source_labelled_profiles": coverage,
        },
        "reserved_external": external,
        "conditional_resources": {
            "GSE227136": "hold: 29/30 Habermann donors overlap; no whole-family independent validation",
            "GSE128033": "hold: cell calling and additional sample/capture mapping",
            "GSE279404": "hold: donor/specimen count and airway-removal boundary",
            "GSE286182": "hold: native/PCLS separation and Yale/Adams independence",
            "other_spatial_resources": "claim-specific extension/reserve only after access, donor, panel, segmentation and anatomy gates",
        },
        "final_lock_required": [
            "Full native unit acquisition, numeric QC and adjudicated cell/library dispositions",
            "Frozen Ensembl/HGNC assets and audited mapping/feature-coverage ledgers",
            "Hierarchical annotation and donor-by-identity endpoint coverage",
            "Expert-reviewed section/ROI anatomy and matched endpoint comparators",
            "Donor independence, estimable models and discovery-representation assessment plan",
            "Versioned endpoint-specific final manifest and Checkpoint A decision",
        ],
    }


def check(data):
    checks = []

    def require(name, condition):
        checks.append({"check": name, "passed": bool(condition)})

    expected = {
        "GSE135893": {"IPF": (12, 19, 2), "Control": (10, 12, 2), "other_ILD": (8, 12, 0)},
        "GSE159354": {"IPF": (4, 20, 8), "Control": (3, 15, 8), "other_ILD": (3, 13, 0)},
        "GSE159585": {"IPF": (6, 6, 2), "Control": (12, 12, 2)},
        "Franzen": {"IPF": (4, 19, 3), "Control": (4, 6, 3)},
    }
    all_uids = []
    for resource, groups in expected.items():
        info = data["resources"][resource]
        for group, values in groups.items():
            observed = info["counts"][group]
            require(f"{resource}/{group} donors, units, processed pilot units", tuple(observed[k] for k in ("source_donors", "units", "pilot_processed_units")) == values)
        units = info["units"]
        all_uids.extend(u["unit_uid"] for u in units)
        require(f"{resource} qualified IDs preserve source identity", all(u["unit_uid"] == f"{resource}::{u['library_or_section_id']}" and u["donor_uid"] == f"{resource}::{u['donor_id']}" for u in units))
        require(f"{resource} annotations remain unreviewed", all(u["project_annotation_status"] == "not_reviewed" for u in units))
        require(f"{resource} no premature endpoint lock", all(u["endpoint_eligibility"] == "not_adjudicated" for u in units))
        require(f"{resource} no pilot flag promoted to formal QC", all(u["full_cohort_qc_disposition"] == "not_adjudicated" for u in units))
    require("134 native intake units", len(all_uids) == 134)
    require("globally unique source-qualified unit IDs", len(set(all_uids)) == len(all_uids))
    require("Checkpoint A remains pending", data["checkpoint_A_passed"] is False)
    require("no formal discovery claimed", data["formal_discovery_started"] is False)
    dep = data["resources"]["GSE159354"]["units"]
    require("five unsorted disease-effect sensitivity libraries retained separately", sum(u["unsorted_disease_effect_sensitivity_only"] for u in dep) == 5)
    require("all native 159354 libraries can contribute to state search", all(u["discovery_state_search_eligible"] for u in dep))
    require("159354 unsorted sensitivity flag matches each source fraction", all(u["unsorted_disease_effect_sensitivity_only"] == (u["facs_fraction"] == "Unsorted") for u in dep))
    require("159354 primary routes match core diagnosis and exact fraction", all(u["primary_comparable_route"] == ("none" if u["analysis_group"] not in ("IPF", "Control") else "epithelial" if u["facs_fraction"] == "EPCAM+" else "stromal" if u["facs_fraction"] == "Triple-" else "none") for u in dep))
    for route in ("epithelial", "stromal"):
        for group in ("IPF", "Control"):
            donors = {u["donor_id"] for u in dep if u["primary_comparable_route"] == route and u["analysis_group"] == group}
            expected_donors = {"F18555", "F18545", "F18536"} if group == "IPF" else {"N1902", "N1857", "N1865"}
            require(f"159354 matched {route}/{group} exact three source donors", donors == expected_donors)
    external = data["reserved_external"]
    for group, count in (("IPF", 32), ("Control", 28), ("COPD", 18)):
        require(f"136831 reserved {group} donor count", len(external["GSE136831"]["donors_by_group"][group]) == count)
    for group, count in (("IPF", 4), ("Control", 8), ("early_IPF_cryobiopsy", 1), ("other_ILD", 4)):
        require(f"122960 reserved {group} donor count", len(external["GSE122960"]["donors_by_branch"][group]) == count)
    require("reserved candidate outcomes untouched", all(x["candidate_outcomes_inspected"] is False for x in external.values()))
    require("no expert anatomy or interface certification fabricated", all(u["expert_anatomy_status"] == "unreviewed" and u["interface_status"] == "not_assessed" and not u["roi_review_records"] for u in data["resources"]["Franzen"]["units"]))
    require("Franzen source procurement confounding retained", all(u["source_procurement"] == ("post_mortem" if u["analysis_group"] == "Control" else "transplant_or_resection") for u in data["resources"]["Franzen"]["units"]))
    coverage = data["source_label_coverage_screen"]["donors_at_least_n_source_labelled_profiles"]["GSE135893"]
    require("Habermann pre-QC mesenchymal 20-cell screen", coverage["Mesenchymal"]["IPF"]["20"] == 9 and coverage["Mesenchymal"]["Control"]["20"] == 6)
    require("Habermann pre-QC epithelial 20-cell screen", coverage["Epithelial"]["IPF"]["20"] == 12 and coverage["Epithelial"]["Control"]["20"] == 10)
    if not all(c["passed"] for c in checks):
        raise ValueError([c for c in checks if not c["passed"]])
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", type=Path, default=PROJECT)
    parser.add_argument("--write", action="store_true", help="Write the intake snapshot before checking it.")
    args = parser.parse_args()
    expected = build(args.evidence_root)
    checks = check(expected)
    target = HERE / "cohort_intake.json"
    if args.write:
        target.write_text(json.dumps(expected, indent=2, ensure_ascii=False) + "\n")
    if read_json(target) != expected:
        raise ValueError("Stored intake differs from metadata/source hashes; investigate before updating.")
    checks.append({"check": "stored intake exactly reproduces from hashed sources", "passed": True})
    result = {"scope": "Metadata intake integrity only; no expression cleaning or anatomical adjudication", "checks_passed": len(checks), "checks": checks}
    (HERE / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Passed {len(checks)} metadata/intake checks; final post-QC allocation remains pending.")


if __name__ == "__main__":
    main()

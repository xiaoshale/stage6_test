"""Check review boundaries on a tiny synthetic capture; run in the server runtime."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import rna_review as review


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    genes = review.MARKERS + ["MT-ND1", "RP11-34P13.3"]
    patterns = [
        "EPCAM KRT8 KRT18 KRT19 KRT17",
        "DCN LUM COL1A1 COL1A2 PDGFRA",
        "MYH11 CNN1 TAGLN ACTA2 MYL9",
        "PROX1 FLT4 CCL21 PDPN",
        "MSLN WT1 KRT8 KRT18",
        "EPCAM KRT8 KRT18 KRT19 KRT17 DCN LUM COL1A1 COL1A2 PDGFRA",
        "ACTA2", "PDPN", "", "KRT5 TP63 KRT15",
    ]
    values = np.zeros((len(patterns), len(genes)), dtype=np.int64)
    for row, pattern in enumerate(patterns):
        for gene in pattern.split():
            values[row, genes.index(gene)] = 20
        if pattern:
            values[row, genes.index("MT-ND1")] = 1
    names = [f"fixture::{i}" for i in range(len(patterns))]
    obs = pd.DataFrame({"donor_id": "fixture_donor", "library_id": "fixture_capture",
                        "source_population": ["Epithelial", "Mesenchymal", "Mesenchymal", "Endothelial", "Other", "Other", "Other", "Other", "Other", "Epithelial"],
                        "source_celltype": ["source_" + str(i) for i in range(len(patterns))],
                        "original_note": "preserve punctuation and original metadata"}, index=names)
    var = pd.DataFrame({"source_symbol": genes}, index=[f"fixture_gene_{i}" for i in range(len(genes))])
    source = args.output / "source.h5ad"
    ad.AnnData(sparse.csr_matrix(values), obs=obs, var=var).write_h5ad(source)
    original_hash = review.sha256(source)
    record = {"path": str(source.resolve()), "resource": "SYNTHETIC_FIXTURE", "library_id": "fixture_capture",
              "donor_id": "fixture_donor", "assay": "scRNA", "fraction": "EPCAM+",
              "capture_identity_status": "proxy_unresolved", "source_label_columns": ["source_population", "source_celltype"]}
    summary = review.process_record(record, args.output / "review", review.sha256(Path(review.__file__)))
    exported = pd.read_csv(args.output / "review/obs_review.tsv.gz", sep="\t", index_col=0, keep_default_na=False)
    checks = {}
    checks["all_profile_ids_and_order_preserved"] = exported.index.tolist() == names
    checks["original_file_unchanged"] = review.sha256(source) == original_hash
    checks["all_original_metadata_retained"] = all(exported[c].astype(str).tolist() == obs[c].astype(str).tolist() for c in obs)
    checks["source_symbol_field_recognized"] = summary["symbol_field"] == "source_symbol"
    checks["source_population_comparison_available"] = summary["source_comparison_field"] == "source_population"
    checks["capture_gate_does_not_assign_identity"] = exported.iloc[1]["review_provisional_core_identity"] == "fibroblast_lineage"
    checks["smooth_muscle_not_forced_to_fibroblast"] = exported.iloc[2]["review_provisional_core_identity"] == "smooth_muscle_like"
    checks["lymphatic_and_mesothelial_guards"] = exported.iloc[3]["review_provisional_core_identity"] == "lymphatic_guard" and exported.iloc[4]["review_provisional_core_identity"] == "mesothelial_guard"
    checks["mixed_identity_is_unresolved"] = exported.iloc[5]["review_provisional_broad_identity"] == "unknown"
    checks["ACTA2_or_PDPN_alone_not_identity"] = all(exported.iloc[i]["review_provisional_broad_identity"] == "unknown" for i in (6, 7))
    checks["zero_profile_retained_and_hard_flagged"] = exported.iloc[8]["review_automated_disposition"] == "HARD_FAIL_ZERO_COUNTS_RETAINED" and summary["hard_fail_zero_profiles"] == 1
    checks["identities_and_eligibility_unadjudicated"] = set(exported["review_expert_identity_status"]) == {"not_adjudicated"} and set(exported["review_endpoint_eligibility"]) == {"not_adjudicated"}
    checks["small_capture_doublets_missing_not_negative"] = summary["scrublet"]["unscored_profiles"] == len(patterns) and exported["review_doublet_flag_prior_0.05"].eq("").all()
    coverage = pd.read_csv(args.output / "review/identity_coverage.tsv", sep="\t")
    checks["unscored_doublet_denominator_explicit"] = coverage["doublet_unscored_profiles"].sum() == len(patterns) and coverage["no_3mad_flags_and_scored_doublet_negative"].sum() == 0
    axis = pd.read_csv(args.output / "review/source_feature_axis.tsv.gz", sep="\t")
    checks["dotted_source_symbol_preserved"] = "RP11-34P13.3" in axis["source_symbol"].values
    base = np.array([0., .1, .5, 1.])
    checks["conditional_prior_odds_identity_and_order"] = np.allclose(review.prior_score(base, .05), base) and np.all(review.prior_score(base, .02) <= base) and np.all(review.prior_score(base, .10) >= base)
    missing = np.zeros(len(review.MARKERS), dtype=bool)
    group = review.marker_groups(pd.Series(["unknown"]), np.zeros((1, len(missing))), np.zeros((1, len(missing)), dtype=bool), missing, "fixture")
    checks["missing_axis_feature_is_coverage_unknown"] = all(r["measurement"] == "not_in_source_export_coverage_unknown" and r["fraction_nonzero"] is None for r in group)
    bad = ad.read_h5ad(source)
    bad.X = bad.X.astype(float)
    bad.X.data[0] = 1.5
    invalid = args.output / "invalid_fractional_counts.h5ad"
    bad.write_h5ad(invalid)
    rejected = False
    try:
        review.process_record({**record, "path": str(invalid.resolve())}, args.output / "invalid_review", "fixture")
    except ValueError as error:
        rejected = "nonnegative integers" in str(error)
    checks["fractional_counts_quarantined"] = rejected
    checks = {name: bool(value) for name, value in checks.items()}
    result = {"fixture": "Synthetic review boundary checks; not biological validation", "passed": int(sum(checks.values())), "total": len(checks), "checks": checks}
    (args.output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    if not all(checks.values()):
        raise AssertionError([k for k, v in checks.items() if not v])
    print(json.dumps(result))


if __name__ == "__main__":
    main()

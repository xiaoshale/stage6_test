"""Read selected native sorted libraries; retain raw counts without integration."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from qc_spatial import read_counts


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    m = json.loads(args.manifest.read_text())
    cohort = next(r for r in m["primary_data_group"] if r["accession"] == "GSE159354")
    libraries = [(r, group) for group in ["primary_libraries", "full_spectrum_libraries"] for r in cohort[group]]
    if len(libraries) != len({r["gsm"] for r, _ in libraries}) or len(libraries) != 16:
        raise ValueError("Expected sixteen distinct core/background libraries")
    started = time.monotonic()
    rows = []
    reference_genes = None
    for source, role in libraries:
        gsm = source["gsm"]
        url, = source["count_or_component_urls"]
        path = args.raw_dir / url.rsplit("/", 1)[1]
        print(datetime.now(timezone.utc).isoformat(), gsm, "reading complete H5", flush=True)
        counts, barcodes, var = read_counts(path.read_bytes())
        counts = counts.astype(np.int64, copy=False)
        counts.eliminate_zeros()
        totals = np.asarray(counts.sum(axis=1)).ravel()
        genes = np.diff(counts.indptr)
        if "gene_symbol" not in var:
            raise ValueError("Missing gene symbols needed for defined mitochondrial metric")
        mt = var["gene_symbol"].str.upper().str.startswith("MT-").to_numpy()
        if not mt.any():
            raise ValueError("No MT- symbols; do not assume zero mitochondrial signal")
        mitochondrial = np.asarray(counts[:, mt].sum(axis=1)).ravel()
        pct = np.divide(100.0 * mitochondrial, totals, out=np.full(len(totals), np.nan), where=totals > 0)
        obs = pd.DataFrame({"source_barcode": barcodes, "donor_id": source["donor_id"], "diagnosis": source["diagnosis"], "library_id": gsm, "facs_fraction": source["facs_fraction"], "pilot_role": role, "total_counts": totals, "n_genes_by_counts": genes, "total_counts_mt": mitochondrial, "pct_counts_mt": pct, "flag_zero_counts": totals == 0, "flag_low_genes_lt200": genes < 200, "flag_high_mt_gt25pct": pct > 25}, index=[gsm + "::" + b for b in barcodes])
        var["mt"] = mt
        obj = ad.AnnData(X=counts, obs=obs, var=var)
        obj.uns["technical_pilot"] = {"source_json": json.dumps(source), "count_file_sha256": sha256(path), "manifest_sha256": sha256(args.manifest), "counts_storage": "X is the original integer matrix; no normalization or QC filtering.", "facs_boundary": "FACS fractions are collection conditions, not verified cell identities or tissue proportions.", "qc_flags": "<200 detected genes, >25% MT-prefix counts or zero counts; descriptive only, not formal inclusion rules."}
        out = args.output_dir / (gsm + "_raw_qc.h5ad")
        obj.write_h5ad(out, compression="gzip")
        reload = ad.read_h5ad(out)
        if reload.shape != obj.shape or not reload.obs_names.equals(obj.obs_names) or not reload.var_names.equals(obj.var_names) or not sparse.issparse(reload.X) or (reload.X != counts).nnz:
            raise ValueError("H5AD raw count/identifier roundtrip mismatch")
        pd.testing.assert_frame_equal(reload.obs, obj.obs, check_dtype=False, check_categorical=False)
        pd.testing.assert_frame_equal(reload.var, obj.var, check_dtype=False, check_categorical=False)
        obj.obs.to_csv(args.output_dir / (gsm + "_cell_qc.csv"))
        same_axis = True if reference_genes is None else var.index.equals(reference_genes)
        if reference_genes is None:
            reference_genes = var.index.copy()
        rows.append({"gsm": gsm, "donor_id": source["donor_id"], "diagnosis": source["diagnosis"], "facs_fraction": source["facs_fraction"], "cells_retained": obj.n_obs, "genes_retained": obj.n_vars, "nonzero_entries": int(counts.nnz), "total_counts": int(totals.sum()), "median_counts": float(np.median(totals)), "median_detected_genes": float(np.median(genes)), "median_pct_mt": float(np.nanmedian(pct)) if np.isfinite(pct).any() else None, "flag_zero_counts": int((totals == 0).sum()), "flag_low_genes_lt200": int((genes < 200).sum()), "flag_high_mt_gt25pct": int((pct > 25).sum()), "gene_axis_equals_first_library": same_axis, "source_count_sha256": sha256(path), "checkpoint_sha256": sha256(out), "roundtrip_passed": True})
        pd.DataFrame(rows).to_csv(args.output_dir / "library_qc.csv", index=False)
        print(gsm, "PASS", obj.shape, flush=True)
    summary = {"accession": "GSE159354", "status": "technical_qc_passed", "libraries": rows, "unique_donors": len({r["donor_id"] for r in rows}), "total_library_cells_retained": sum(r["cells_retained"] for r in rows), "all_source_filtered_cells_and_genes_retained": True, "cell_identity_boundary": "Cell Ranger-filtered source barcodes retained; FACS EPCAM/Triple-negative/CD31/CD45 labels do not establish biological cell states or sufficient fibroblast recovery.", "comparability_boundary": "No absolute abundance inference across sorted fractions, no batch integration or cohort disease tests. Each library has a separate checkpoint; differing gene axes are not silently zero-filled.", "formal_biological_analysis_ready": False, "elapsed_seconds": time.monotonic()-started, "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024, "script_sha256": sha256(Path(__file__)), "qc_filtering": "None", "not_performed": ["Normalization", "Clustering", "Differential expression", "Niche or target analysis"]}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()

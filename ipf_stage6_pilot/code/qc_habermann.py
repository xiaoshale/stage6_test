#!/usr/bin/env python3
"""Audit GSE135893 original counts and retain the manifest's annotated pilot cells.

This technical pilot never filters cells/genes, normalizes counts, discovers cell
states, clusters, calculates disease effects or tests niches/targets. Author
annotations remain explicitly source-provided, rather than new annotations.
"""

import argparse
import gc
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import time


FILES = {
    "counts": "GSE135893_matrix.mtx.gz",
    "genes": "GSE135893_genes.tsv.gz",
    "barcodes": "GSE135893_barcodes.tsv.gz",
    "metadata": "GSE135893_IPF_metadata.csv.gz",
}


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_selection(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    candidates = [entry for entry in manifest["primary_data_group"]
                  if entry["accession"] == "GSE135893"]
    if len(candidates) != 1:
        raise ValueError("Manifest must contain one GSE135893 primary-data entry")
    entry = candidates[0]
    donors = entry["selected_donors"]
    donor_ids = [donor["donor_id"] for donor in donors]
    if len(donor_ids) != 4 or len(set(donor_ids)) != 4:
        raise ValueError("The approved pilot requires four unique source donors")
    if sorted(donor["diagnosis"] for donor in donors) != ["Control", "Control", "IPF", "IPF"]:
        raise ValueError("The approved pilot requires two IPF and two control donors")
    for donor in donors:
        if not donor["library_ids"]:
            raise ValueError(f"Missing source library mapping for {donor['donor_id']}")
    return entry


def log(message):
    print(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), message, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-genes-flag", type=int, default=200)
    parser.add_argument("--max-mito-pct-flag", type=float, default=25.0)
    args = parser.parse_args()
    if args.threads < 1 or args.min_genes_flag < 0 or not 0 <= args.max_mito_pct_flag <= 100:
        parser.error("Invalid thread limit or technical flag threshold")
    raw_dir = args.raw_dir.resolve(strict=True)
    manifest_path = args.manifest.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        parser.error("--output-dir must be absent or empty; previous outputs are preserved")
    inputs = {role: (raw_dir / filename).resolve(strict=True) for role, filename in FILES.items()}
    selection = load_selection(manifest_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[name] = str(args.threads)
    started = time.perf_counter()

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import io, sparse
    from threadpoolctl import threadpool_limits

    log("Reading complete source axes and author metadata")
    genes = pd.read_csv(inputs["genes"], sep="\t", header=None, dtype=str, keep_default_na=False)
    barcode_table = pd.read_csv(inputs["barcodes"], sep="\t", header=None, dtype=str,
                                keep_default_na=False)
    if genes.shape[1] != 1 or barcode_table.shape[1] != 1:
        raise ValueError("Expected the original one-column gene-symbol and barcode axes")
    gene_symbols = genes.iloc[:, 0]
    barcodes = pd.Index(barcode_table.iloc[:, 0], name="source_barcode")
    metadata = pd.read_csv(inputs["metadata"], index_col=0, keep_default_na=False)
    metadata.index = metadata.index.astype(str)
    if not barcodes.is_unique or not metadata.index.is_unique:
        raise ValueError("Source matrix and metadata barcode identifiers must be unique")
    if gene_symbols.eq("").any() or barcodes.isin([""]).any():
        raise ValueError("Source axes contain empty identifiers")
    required = {"Sample_Name", "Diagnosis", "orig.ident", "population", "nCount_RNA", "nFeature_RNA"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Source metadata lacks columns: {sorted(required - set(metadata.columns))}")
    all_metadata_positions = barcodes.get_indexer(metadata.index)
    if (all_metadata_positions < 0).any():
        raise ValueError("Some author-annotated barcodes are absent from the original matrix axis")

    requested_donors = [donor["donor_id"] for donor in selection["selected_donors"]]
    selected_metadata = metadata.loc[metadata["Sample_Name"].isin(requested_donors)]
    for donor in selection["selected_donors"]:
        rows = selected_metadata.loc[selected_metadata["Sample_Name"].eq(donor["donor_id"])]
        if rows.empty or set(rows["Diagnosis"]) != {donor["diagnosis"]}:
            raise ValueError(f"Missing donor or conflicting diagnosis: {donor['donor_id']}")
        if set(rows["orig.ident"]) != set(donor["library_ids"]):
            raise ValueError(f"Conflicting source library mapping: {donor['donor_id']}")
        expected_broad = donor.get("provided_broad_lineage_cell_counts")
        if expected_broad is not None and rows["population"].value_counts().to_dict() != expected_broad:
            raise ValueError(f"Source broad-label counts differ from manifest: {donor['donor_id']}")
    # Selection is exclusively by exact author metadata barcode, never inferred
    # from a library prefix. Matrix columns lacking author annotation are audited
    # globally, but are not reclassified as called cells or added to this subset.
    selected_positions = np.sort(barcodes.get_indexer(selected_metadata.index))
    selected_barcodes = barcodes.take(selected_positions)
    selected_metadata = selected_metadata.reindex(selected_barcodes)
    if selected_metadata["Sample_Name"].isna().any():
        raise ValueError("Barcode join failed for the selected source-annotated subset")

    log("Hashing original inputs and reading the complete sparse MTX")
    provenance = {role: {"file": str(path), "bytes": path.stat().st_size,
                          "sha256": file_sha256(path)} for role, path in inputs.items()}
    with threadpool_limits(limits=args.threads):
        with gzip.open(inputs["counts"], "rb") as handle:
            global_counts = io.mmread(handle)
            if handle.read(1):
                raise ValueError("Unexpected trailing payload after the Matrix Market data")
    if not sparse.issparse(global_counts):
        raise ValueError("The source Matrix Market payload is not sparse")
    if global_counts.shape != (len(gene_symbols), len(barcodes)):
        raise ValueError("Full MTX shape disagrees with its gene/barcode axes")
    expected_header = selection.get("full_matrix_header", {})
    actual_header = {"genes": global_counts.shape[0], "barcode_columns": global_counts.shape[1],
                     "nonzero_entries": int(global_counts.nnz)}
    for key, actual in actual_header.items():
        if key in expected_header and expected_header[key] != actual:
            raise ValueError(f"Full MTX {key} differs from the checked manifest header")
    values = global_counts.data
    if not np.issubdtype(values.dtype, np.integer):
        raise ValueError("The original counts must have an integer sparse data type")
    if values.size and values.min() < 0:
        raise ValueError("The original counts contain negative entries")
    source_dtype = str(values.dtype)
    if not values.size or values.max() <= np.iinfo(np.int32).max:
        global_counts.data = values.astype(np.int32, copy=False)
    del values
    log(f"Selecting {len(selected_barcodes)} exact annotated barcodes; no QC filtering")
    # A complete sparse read plus CSC conversion fits the intended 32 GB budget;
    # no dense global matrix or cell-by-gene copy is constructed.
    global_csc = global_counts.tocsc()
    del global_counts
    gc.collect()
    counts = global_csc[:, selected_positions].T.tocsr()
    del global_csc
    gc.collect()
    counts.sum_duplicates()
    counts.eliminate_zeros()
    counts.sort_indices()

    obs = selected_metadata.add_prefix("source_").copy()
    obs["source_barcode"] = selected_barcodes.to_numpy()
    obs["donor_id"] = obs["source_Sample_Name"]
    obs["diagnosis"] = obs["source_Diagnosis"]
    obs["library_id"] = obs["source_orig.ident"]
    obs["source_broad_annotation"] = obs["source_population"]
    obs["source_metadata_present"] = True
    obs.index = pd.Index(["GSE135893::" + barcode for barcode in selected_barcodes], name="cell_id")
    var = pd.DataFrame({"source_gene_symbol": gene_symbols.to_numpy(),
                        "source_gene_row_1based": np.arange(1, len(gene_symbols) + 1)},
                       index=pd.Index(gene_symbols, name="gene_symbol_unique"))
    var["mt"] = var["source_gene_symbol"].str.upper().str.startswith("MT-")
    if not var["mt"].any():
        raise ValueError("No MT-prefixed source symbols: mitochondrial metric cannot be assumed zero")
    total_counts = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    n_genes = np.diff(counts.indptr)
    mito_counts = np.asarray(counts[:, var["mt"].to_numpy()].sum(axis=1)).ravel().astype(np.int64)
    pct_mito = np.full(len(total_counts), np.nan, dtype=float)
    np.divide(100.0 * mito_counts, total_counts, out=pct_mito, where=total_counts > 0)
    obs["total_counts"] = total_counts
    obs["n_genes_by_counts"] = n_genes
    obs["total_counts_mt"] = mito_counts
    obs["pct_counts_mt"] = pct_mito
    obs["flag_low_genes"] = n_genes < args.min_genes_flag
    obs["flag_high_mito"] = pct_mito > args.max_mito_pct_flag
    obs["flag_zero_counts"] = total_counts == 0
    obs["flag_any_technical"] = obs[["flag_low_genes", "flag_high_mito", "flag_zero_counts"]].any(axis=1)
    obs["source_nCount_RNA_matches_computed"] = pd.to_numeric(obs["source_nCount_RNA"]).eq(total_counts)
    obs["source_nFeature_RNA_matches_computed"] = pd.to_numeric(obs["source_nFeature_RNA"]).eq(n_genes)
    var["total_counts"] = np.asarray(counts.sum(axis=0)).ravel().astype(np.int64)
    var["n_cells_by_counts"] = np.asarray((counts > 0).sum(axis=0)).ravel().astype(np.int64)
    duplicate_gene_symbols = int(var.index.duplicated().sum())
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X
    thresholds = {"n_genes_strictly_below": args.min_genes_flag,
                  "pct_counts_mt_strictly_above": args.max_mito_pct_flag,
                  "zero_counts": True, "cells_or_genes_filtered": False}
    adata.uns["technical_pilot"] = {
        "accession": "GSE135893", "counts_storage": "X and layers/counts contain identical original integer counts",
        "selection": "Exact source-annotated barcodes for the four manifest donors; all selected cells and genes retained",
        "annotation_status": "Source-provided broad and cell-type annotations copied without reannotation",
        "flag_thresholds_json": json.dumps(thresholds),
        "selected_manifest_entry_json": json.dumps(selection),
        "manifest_sha256": file_sha256(manifest_path),
        "source_provenance_json": json.dumps(provenance),
        "no_biological_analysis": True,
    }
    adata.obs.to_csv(output_dir / "QC.csv", index_label="cell_id")
    grouping = ["donor_id", "diagnosis", "library_id"]
    donor_qc = adata.obs.groupby(grouping, observed=True, dropna=False).agg(
        cells=("source_barcode", "size"), median_total_counts=("total_counts", "median"),
        median_detected_genes=("n_genes_by_counts", "median"), median_pct_counts_mt=("pct_counts_mt", "median"),
        flagged_low_genes=("flag_low_genes", "sum"), flagged_high_mito=("flag_high_mito", "sum"),
        flagged_zero_counts=("flag_zero_counts", "sum"), flagged_any_technical=("flag_any_technical", "sum"))
    donor_qc.to_csv(output_dir / "donor_library_QC.csv")
    broad_qc = adata.obs.groupby(grouping + ["source_broad_annotation"], observed=True, dropna=False).agg(
        cells=("source_barcode", "size"), median_total_counts=("total_counts", "median"),
        median_detected_genes=("n_genes_by_counts", "median"), median_pct_counts_mt=("pct_counts_mt", "median"),
        flagged_any_technical=("flag_any_technical", "sum"))
    broad_qc.to_csv(output_dir / "source_broad_annotation_QC.csv")
    log("Writing and independently reloading sparse counts.h5ad")
    output_path = output_dir / "counts.h5ad"
    adata.write_h5ad(output_path, compression="gzip")
    reloaded = ad.read_h5ad(output_path)
    checks = {
        "shape": reloaded.shape == adata.shape,
        "cell_ids": reloaded.obs_names.equals(adata.obs_names),
        "gene_ids": reloaded.var_names.equals(adata.var_names),
        "X_original_counts": sparse.issparse(reloaded.X) and (reloaded.X != adata.X).nnz == 0,
        "counts_layer": sparse.issparse(reloaded.layers["counts"]) and (reloaded.layers["counts"] != adata.X).nnz == 0,
        "all_selected_annotated_cells_retained": adata.n_obs == len(selected_metadata),
        "all_source_genes_retained": adata.n_vars == len(gene_symbols),
        "no_embeddings_or_neighbor_graphs": not reloaded.obsm and not reloaded.obsp,
    }
    pd.testing.assert_frame_equal(reloaded.obs, adata.obs, check_dtype=False, check_categorical=False)
    pd.testing.assert_frame_equal(reloaded.var, adata.var, check_dtype=False, check_categorical=False)
    checks["all_obs_and_var_values"] = True
    if not all(checks.values()):
        raise ValueError(f"Sparse count round-trip verification failed: {checks}")
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_rss_bytes = peak_rss if platform.system() == "Darwin" else peak_rss * 1024
    summary = {
        "accession": "GSE135893", "technical_execution_complete": True, "formal_biological_analysis_ready": False,
        "selected_donors": requested_donors, "selected_shape_cells_genes": list(adata.shape),
        "selected_nnz": int(adata.X.nnz), "selected_total_counts": int(total_counts.sum()),
        "selected_sparse_dtype": str(adata.X.dtype), "source_sparse_dtype": source_dtype,
        "global_counts_axes_audit": {**actual_header, "complete_sparse_payload_read": True,
                                     "source_metadata_rows": len(metadata),
                                     "all_metadata_barcodes_map_to_unique_matrix_columns": True,
                                     "matrix_columns_without_this_author_annotation": len(barcodes) - len(metadata),
                                     "unannotated_columns_included_in_pilot": 0,
                                     "matrix_axis_not_assumed_to_be_called_cells": True},
        "mitochondrial_gene_definition": "Case-insensitive MT- prefix in original gene-symbol axis",
        "mitochondrial_gene_rows": int(var["mt"].sum()), "duplicate_source_gene_symbol_rows": duplicate_gene_symbols,
        "flag_thresholds": thresholds,
        "technical_flag_counts": {name: int(adata.obs[name].sum()) for name in
                                  ["flag_low_genes", "flag_high_mito", "flag_zero_counts", "flag_any_technical"]},
        "source_metric_disagreements": {"nCount_RNA": int((~obs["source_nCount_RNA_matches_computed"]).sum()),
                                         "nFeature_RNA": int((~obs["source_nFeature_RNA_matches_computed"]).sum())},
        "source_broad_annotation_boundary": "Author labels only, not verified fibroblast states, tissue abundance, spatial proximity or new annotation",
        "reload_validation": checks,
        "original_inputs": provenance,
        "manifest": {"file": str(manifest_path), "sha256": file_sha256(manifest_path)},
        "script_sha256": file_sha256(Path(__file__).resolve()),
        "versions": {name: importlib.metadata.version(name) for name in
                     ["anndata", "numpy", "scipy", "pandas", "h5py", "threadpoolctl"]},
        "python": platform.python_version(), "threads": args.threads,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "peak_rss_bytes": int(peak_rss_bytes),
        "counts_h5ad": {"file": output_path.name, "sha256": file_sha256(output_path)},
        "not_performed": ["Cell or gene QC filtering", "Normalization or log transform", "Cell-state discovery",
                          "PCA or clustering", "Differential expression", "Niche or target analysis"],
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"PASS: {adata.n_obs} source-annotated cells, {adata.n_vars} genes; all cells retained")


if __name__ == "__main__":
    main()

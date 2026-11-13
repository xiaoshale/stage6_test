#!/usr/bin/env python3
"""Run a technical PBMC3k smoke test; this is not an IPF analysis.

The QC cutoffs and clustering settings below are demonstration settings, not
Stage 6 biological criteria. Clusters are unannotated technical outputs.
"""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import time


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    input_path = args.input.resolve(strict=True)
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        parser.error("--output-dir must be absent or empty; previous runs are preserved")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set limits before scientific libraries initialize their thread pools.
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS",
    ):
        os.environ[name] = str(args.threads)
    os.environ["MPLBACKEND"] = "Agg"
    started = time.perf_counter()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numba
    import numpy as np
    import scanpy as sc
    from scipy import sparse

    numba.set_num_threads(args.threads)
    sc.settings.n_jobs = args.threads
    sc.settings.verbosity = 2
    input_hash = file_sha256(input_path)
    adata = sc.read_h5ad(input_path)
    original_shape = list(adata.shape)
    if not adata.obs_names.is_unique:
        raise ValueError("Input cell identifiers must be unique")
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Input must contain finite, nonnegative raw counts")
    if not np.allclose(values, np.rint(values), rtol=0, atol=1e-8):
        raise ValueError("Input values are not integer-like raw counts")

    adata.var["source_gene_name"] = adata.var_names.astype(str)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True,
    )
    adata.uns["technical_smoke_test"] = {
        "expected_dataset": "PBMC3k public 10x Genomics example",
        "interpretation": "Technical test only; not lung or IPF data",
        "source_input_sha256": input_hash,
        "seed": args.seed,
        "threads": args.threads,
    }
    # This preserves all input cells, genes, and counts before any filtering.
    adata.write_h5ad(output_dir / "raw_counts_checkpoint.h5ad", compression="gzip")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, column, label in zip(
        axes,
        ("n_genes_by_counts", "total_counts", "pct_counts_mt"),
        ("Detected genes per cell", "Counts per cell", "Mitochondrial counts (%)"),
    ):
        ax.hist(adata.obs[column], bins=45, color="#367F9D", edgecolor="white")
        ax.set_xlabel(label)
        ax.set_ylabel("Cells")
    axes[0].axvline(200, color="#BA4A3B", linestyle="--")
    axes[2].axvline(20, color="#BA4A3B", linestyle="--")
    fig.suptitle("Technical smoke test | PBMC3k | Not IPF\nPre-filter QC; tutorial-only cutoffs")
    fig.tight_layout()
    fig.savefig(output_dir / "qc.png", dpi=160)
    plt.close(fig)

    # These explicit demonstration thresholds are not patient-analysis criteria.
    keep = (adata.obs["n_genes_by_counts"] >= 200) & (adata.obs["pct_counts_mt"] < 20)
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    filtered_shape = list(adata.shape)
    if min(adata.shape) < 4:
        raise ValueError("Too few cells or genes remain for this smoke test")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat", n_top_genes=min(2000, adata.n_vars),
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    n_pcs = min(30, adata.n_obs - 1, n_hvg - 1)
    if n_pcs < 2:
        raise ValueError("Too few variable genes for a useful PCA smoke test")
    sc.pp.pca(
        adata, n_comps=n_pcs, mask_var="highly_variable",
        svd_solver="arpack", random_state=args.seed,
    )
    sc.pp.neighbors(
        adata, n_neighbors=min(10, adata.n_obs - 1),
        n_pcs=n_pcs, random_state=args.seed,
    )
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(
        adata, resolution=0.5, random_state=args.seed,
        flavor="igraph", directed=False, n_iterations=2, key_added="leiden",
    )
    if not np.isfinite(adata.obsm["X_pca"]).all():
        raise ValueError("PCA contains non-finite values")
    if not np.isfinite(adata.obsm["X_umap"]).all():
        raise ValueError("UMAP contains non-finite values")

    fig, ax = plt.subplots(figsize=(7.4, 5.5))
    categories = list(adata.obs["leiden"].cat.categories)
    palette = plt.get_cmap("tab20")
    for index, category in enumerate(categories):
        selected = adata.obs["leiden"].eq(category).to_numpy()
        xy = adata.obsm["X_umap"][selected]
        ax.scatter(xy[:, 0], xy[:, 1], s=5, color=palette(index % 20),
                   label=f"Cluster {category}", linewidths=0, rasterized=True)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Technical smoke test | PBMC3k | Not IPF\nUnannotated Leiden clusters")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, markerscale=2)
    fig.tight_layout()
    fig.savefig(output_dir / "umap.png", dpi=180)
    plt.close(fig)

    columns = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "leiden"]
    adata.obs[columns].to_csv(output_dir / "cell_qc_clusters.csv", index_label="cell_id")
    final_path = output_dir / "processed.h5ad"
    adata.write_h5ad(final_path, compression="gzip")
    reloaded = sc.read_h5ad(final_path)
    counts_before, counts_after = adata.layers["counts"], reloaded.layers["counts"]
    if sparse.issparse(counts_before) or sparse.issparse(counts_after):
        counts_match = (sparse.csr_matrix(counts_before) != sparse.csr_matrix(counts_after)).nnz == 0
    else:
        counts_match = np.array_equal(counts_before, counts_after)
    validation = {
        "raw_counts_finite_nonnegative_integer_like": True,
        "input_cell_ids_unique": True,
        "roundtrip_shape": reloaded.shape == adata.shape,
        "roundtrip_cell_ids": reloaded.obs_names.equals(adata.obs_names),
        "roundtrip_gene_ids": reloaded.var_names.equals(adata.var_names),
        "roundtrip_raw_counts": bool(counts_match),
        "roundtrip_cluster_labels": np.array_equal(
            reloaded.obs["leiden"].astype(str), adata.obs["leiden"].astype(str),
        ),
        "roundtrip_pca": np.array_equal(reloaded.obsm["X_pca"], adata.obsm["X_pca"]),
        "roundtrip_umap": np.array_equal(reloaded.obsm["X_umap"], adata.obsm["X_umap"]),
        "embeddings_finite": bool(
            np.isfinite(reloaded.obsm["X_pca"]).all()
            and np.isfinite(reloaded.obsm["X_umap"]).all()
        ),
    }
    if not all(validation.values()):
        raise RuntimeError(f"Saved output validation failed: {validation}")

    packages = [
        "scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib",
        "scikit-learn", "umap-learn", "numba", "igraph", "h5py",
    ]
    summary = {
        "status": "PASS",
        "purpose": "Technical environment and execution test; no IPF or biological conclusions",
        "expected_dataset": "PBMC3k (blood-cell example, not a lung cohort)",
        "input": str(input_path),
        "input_sha256": input_hash,
        "acquisition_provenance": "See accompanying download record; this script does not infer provenance from expression",
        "shape_before_filtering": original_shape,
        "shape_after_filtering": filtered_shape,
        "highly_variable_genes": n_hvg,
        "pca_components": n_pcs,
        "cluster_counts": {
            str(key): int(value)
            for key, value in adata.obs["leiden"].value_counts(sort=False).items()
        },
        "demonstration_parameters": {
            "min_genes_per_cell": 200, "max_mito_percent_exclusive": 20,
            "min_cells_per_gene": 3, "normalization_target": 10000,
            "hvg_flavor": "seurat", "hvg_requested": min(2000, adata.n_vars),
            "neighbors": min(10, adata.n_obs - 1), "leiden_resolution": 0.5,
            "leiden_flavor": "igraph", "leiden_iterations": 2,
            "threads": args.threads, "seed": args.seed,
        },
        "validation": validation,
        "wall_seconds": round(time.perf_counter() - started, 3),
        "peak_rss_MiB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 3)
        if platform.system() == "Linux" else None,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name) for name in packages},
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    artifacts = sorted(path for path in output_dir.iterdir() if path.is_file())
    success = {
        "status": "PASS", "summary": summary_path.name,
        "output_sha256": {path.name: file_sha256(path) for path in artifacts},
    }
    (output_dir / "SUCCESS.json").write_text(json.dumps(success, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

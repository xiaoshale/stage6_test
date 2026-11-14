"""Stream the native GSE159585 CSV into an unfiltered technical-pilot H5AD.

Only the manifest-selected donors are retained. The source CSV is gene by cell;
one complete gene row is parsed at a time, without a global dense expression
table. Sparse CSC columns correspond directly to successive source gene rows.
"""

from __future__ import annotations

import argparse
from array import array
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import platform
import resource
import sys
import time


ACCESSION = "GSE159585"
COUNTS_FILENAME = "GSE159585_snRNAseq_raw_counts.csv.gz"
METADATA_FILENAME = "GSE159585_snRNAseq_metadata.csv.gz"


def log(stage: str, **values) -> None:
    print(json.dumps({"stage": stage, **values}), flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_samples(manifest: dict) -> dict:
    component = manifest["technical_robustness_component"]
    if component["accession"] != ACCESSION:
        raise ValueError("The robustness component must identify GSE159585.")
    samples = {}
    for sample in component["selected_samples"]:
        donor = sample["patient_ID"]
        if donor in samples:
            raise ValueError(f"Duplicate manifest patient_ID: {donor}")
        procurement = [
            value.split(":", 1)[1].strip()
            for value in sample["GEO_characteristics"]
            if value.lower().startswith("sample type:")
        ]
        if len(procurement) != 1:
            raise ValueError(f"Missing or ambiguous procurement for {donor}")
        samples[donor] = {**sample, "procurement": procurement[0]}
    if not samples:
        raise ValueError("No manifest-selected samples.")
    return samples


def read_metadata(path: Path, samples: dict) -> tuple[set, dict]:
    all_barcodes = set()
    selected = {}
    required = {"", "sample_ID", "patient_ID", "condition", "cellType", "subcluster"}
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Unexpected metadata columns: {reader.fieldnames}")
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(row[field] is None for field in required):
                raise ValueError(f"Malformed metadata row {line_number}")
            barcode = row[""]
            if not barcode or barcode in all_barcodes:
                raise ValueError(f"Blank or duplicate metadata barcode at row {line_number}")
            all_barcodes.add(barcode)
            donor = row["patient_ID"]
            if donor in samples:
                sample = samples[donor]
                if row["condition"] != sample["condition"]:
                    raise ValueError(f"Condition mismatch for {donor}")
                if row["sample_ID"] not in sample["metadata_sample_IDs"]:
                    raise ValueError(f"Sample ID mismatch for {donor}: {row['sample_ID']}")
                selected[barcode] = row
    present = {row["patient_ID"] for row in selected.values()}
    if present != set(samples):
        raise ValueError(f"Selected donors missing from metadata: {sorted(set(samples) - present)}")
    return all_barcodes, selected


def sparse_fingerprint(matrix) -> str:
    """Hash every stored integer and sparse index without a large bytes copy."""
    digest = hashlib.sha256(str(matrix.shape).encode())
    for values in (matrix.data, matrix.indices, matrix.indptr):
        digest.update(str(values.dtype).encode())
        digest.update(memoryview(values))
    return digest.hexdigest()


def stream_counts(path: Path, all_metadata_barcodes: set, metadata: dict):
    import numpy as np
    from scipy import sparse

    stored_values = array("q")
    stored_indices = array("i")
    column_pointers = array("q", [0])
    genes = []
    genes_nonzero_cells = []
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if len(header) < 2 or header[0] != "":
            raise ValueError("Expected an empty gene-index header followed by cell IDs.")
        global_barcodes = header[1:]
        global_set = set(global_barcodes)
        if len(global_set) != len(global_barcodes):
            raise ValueError("Duplicate cell IDs in the original count header.")
        if global_set != all_metadata_barcodes:
            raise ValueError(
                "Count/metadata global barcode sets differ: "
                f"count_only={sorted(global_set - all_metadata_barcodes)[:10]}, "
                f"metadata_only={sorted(all_metadata_barcodes - global_set)[:10]}"
            )
        selected_columns = [
            index for index, barcode in enumerate(header) if index and barcode in metadata
        ]
        barcodes = [header[index] for index in selected_columns]
        if len(barcodes) != len(metadata):
            raise ValueError("Selected metadata cells do not map one-to-one to count columns.")
        n_cells = len(barcodes)
        total_counts = np.zeros(n_cells, dtype=np.int64)
        detected_genes = np.zeros(n_cells, dtype=np.int32)
        mt_counts = np.zeros(n_cells, dtype=np.int64)
        expected_width = len(header)
        log("count_header_validated", global_cells=len(global_barcodes), selected_cells=n_cells)

        for gene_number, row in enumerate(reader, start=1):
            if len(row) != expected_width or not row[0]:
                raise ValueError(f"Malformed count row {gene_number}: {len(row)} fields")
            gene = row[0]
            try:
                values = np.fromiter(
                    (int(row[index]) for index in selected_columns),
                    dtype=np.int64,
                    count=n_cells,
                )
            except (ValueError, OverflowError) as error:
                raise ValueError(f"Non-integer or overflowing selected count at gene {gene}") from error
            if np.any(values < 0):
                raise ValueError(f"Negative selected count at gene {gene}")
            nonzero = np.flatnonzero(values)
            stored_values.extend(values[nonzero].tolist())
            stored_indices.extend(nonzero.tolist())
            column_pointers.append(len(stored_values))
            genes.append(gene)
            genes_nonzero_cells.append(len(nonzero))
            total_counts += values
            detected_genes += values > 0
            if gene.upper().startswith("MT-"):
                mt_counts += values
            if gene_number % 2000 == 0:
                log("count_rows", genes=gene_number, selected_nonzero=len(stored_values))

        # Iterating the gzip reader to EOF also checks its trailer/CRC; a
        # truncated gzip raises before this point even if its last row is valid.
        if not genes:
            raise ValueError("The count file has no gene rows.")

    matrix = sparse.csc_matrix(
        (
            np.frombuffer(stored_values, dtype=np.int64),
            np.frombuffer(stored_indices, dtype=np.int32),
            np.frombuffer(column_pointers, dtype=np.int64),
        ),
        shape=(n_cells, len(genes)),
    )
    matrix.check_format(full_check=True)
    if np.any(total_counts < 0):
        raise ValueError("Integer overflow in selected cell totals.")
    log("count_eof", genes=len(genes), selected_cells=n_cells, selected_nonzero=matrix.nnz)
    return matrix, genes, barcodes, total_counts, detected_genes, mt_counts, genes_nonzero_cells


def write_table(path: Path, fields: list, rows: list) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(args, summary: dict) -> None:
    import anndata as ad
    import numpy as np

    counts_path = args.raw_dir / COUNTS_FILENAME
    metadata_path = args.raw_dir / METADATA_FILENAME
    manifest = json.loads(args.manifest.read_text())
    samples = selected_samples(manifest)
    all_barcodes, metadata = read_metadata(metadata_path, samples)
    log("metadata", global_cells=len(all_barcodes), selected_cells=len(metadata), donors=sorted(samples))
    summary["input_files"] = {
        "counts": {"path": str(counts_path), "compressed_bytes": counts_path.stat().st_size},
        "metadata": {"path": str(metadata_path), "sha256": sha256_file(metadata_path)},
        "manifest": {"path": str(args.manifest), "sha256": sha256_file(args.manifest)},
    }
    matrix, genes, barcodes, totals, detected, mt, gene_cells = stream_counts(
        counts_path, all_barcodes, metadata
    )
    summary["full_count_file_read_to_eof"] = True
    summary["gzip_trailer_crc_passed"] = True
    summary["input_files"]["counts"]["sha256"] = sha256_file(counts_path)
    pct_mt = np.divide(100.0 * mt, totals, out=np.full(len(totals), np.nan), where=totals != 0)
    obs = {
        "source_cell_id": barcodes,
        "donor_id": [metadata[barcode]["patient_ID"] for barcode in barcodes],
        "source_sample_id": [metadata[barcode]["sample_ID"] for barcode in barcodes],
        "diagnosis": [metadata[barcode]["condition"] for barcode in barcodes],
        "source_coarse_cell_type": [metadata[barcode]["cellType"] for barcode in barcodes],
        "source_subcluster": [metadata[barcode]["subcluster"] for barcode in barcodes],
        "gsm": [samples[metadata[barcode]["patient_ID"]]["GSM"] for barcode in barcodes],
        "procurement": [samples[metadata[barcode]["patient_ID"]]["procurement"] for barcode in barcodes],
        "total_counts": totals,
        "n_genes_by_counts": detected,
        "total_counts_mt": mt,
        "pct_counts_mt": pct_mt,
        "flag_zero_counts": totals == 0,
        "flag_low_genes_lt200": detected < 200,
        "flag_high_mt_gt25pct": pct_mt > 25.0,
    }
    adata = ad.AnnData(
        X=matrix,
        obs=obs,
        var={"source_gene_name": genes, "n_cells_by_counts": gene_cells},
    )
    adata.obs_names = barcodes
    adata.var_names = genes
    duplicate_gene_rows = len(genes) - len(set(genes))
    if duplicate_gene_rows:
        adata.var_names_make_unique()
    flag_definitions = {
        "flag_zero_counts": "Selected cell has zero total source counts.",
        "flag_low_genes_lt200": "Fewer than 200 genes with nonzero source counts; descriptive flag only.",
        "flag_high_mt_gt25pct": "More than 25% of counts in source gene labels beginning MT-; descriptive flag only. Mitochondrial percentage is NaN when total counts are zero.",
    }
    confounding = (
        "All six native IPF subjects are explants and all twelve non-COVID controls are autopsies. "
        "Disease and procurement are completely confounded in this cohort; integration or a "
        "covariate cannot separate their effects. This is a separate technical/cross-assay "
        "robustness component, not equally clean third disease-discovery evidence. "
        "Cross-study donor independence remains unverified."
    )
    adata.uns["technical_pilot"] = {
        "accession": ACCESSION,
        "counts": "Original nonnegative integer source counts in X; all selected cells and gene rows retained.",
        "flags": flag_definitions,
        "source_annotations": "Published coarse/subcluster labels preserved; no new cell annotation inferred.",
        "confounding": confounding,
        "duplicate_gene_rows": duplicate_gene_rows,
        "analysis_boundary": "Technical QC only; no differential expression, PCA, clustering, niches or target inference.",
    }
    before_hash = sparse_fingerprint(adata.X)
    output_h5ad = args.output_dir / "GSE159585_raw_qc.h5ad"
    log("write_h5ad", cells=adata.n_obs, genes=adata.n_vars, sparse_format=adata.X.format)
    adata.write_h5ad(output_h5ad, compression="gzip")
    log("reload_h5ad")
    reloaded = ad.read_h5ad(output_h5ad)
    if reloaded.shape != adata.shape or reloaded.X.format != "csc":
        raise ValueError("Reloaded matrix shape or sparse format differs.")
    if sparse_fingerprint(reloaded.X) != before_hash:
        raise ValueError("Reload changed original sparse integer values or indices.")
    if list(reloaded.obs_names) != barcodes or list(reloaded.var["source_gene_name"]) != genes:
        raise ValueError("Reload changed the source barcode or gene axis.")
    for column in ["donor_id", "source_sample_id", "diagnosis", "source_coarse_cell_type", "source_subcluster", "gsm", "procurement"]:
        if list(reloaded.obs[column]) != list(adata.obs[column]):
            raise ValueError(f"Reload changed source metadata column {column}")
    if not np.array_equal(np.asarray(reloaded.X.sum(axis=1)).ravel(), totals):
        raise ValueError("Reloaded cell totals differ from the streamed counts.")
    cell_rows = []
    for barcode, row in reloaded.obs.iterrows():
        cell_rows.append({
            "cell_id": barcode,
            **{name: row[name] for name in [
                "donor_id", "source_sample_id", "diagnosis", "gsm", "procurement",
                "source_coarse_cell_type", "total_counts", "n_genes_by_counts",
                "pct_counts_mt", "flag_zero_counts", "flag_low_genes_lt200", "flag_high_mt_gt25pct",
            ]},
        })
    write_table(args.output_dir / "cell_qc.tsv", list(cell_rows[0]), cell_rows)
    donor_rows = []
    donor_detail = {}
    donor_ids = np.asarray(obs["donor_id"])
    for donor, sample in samples.items():
        selected = donor_ids == donor
        lineage_counts = dict(Counter(np.asarray(obs["source_coarse_cell_type"])[selected]))
        defined_mt_percentages = pct_mt[selected][np.isfinite(pct_mt[selected])]
        donor_row = {
            "donor_id": donor,
            "diagnosis": sample["condition"],
            "gsm": sample["GSM"],
            "procurement": sample["procurement"],
            "retained_source_cells": int(selected.sum()),
            "median_total_counts": float(np.median(totals[selected])),
            "median_detected_genes": float(np.median(detected[selected])),
            "median_pct_counts_mt": float(np.median(defined_mt_percentages)) if len(defined_mt_percentages) else None,
            "flag_zero_counts_cells": int(np.sum(obs["flag_zero_counts"][selected])),
            "flag_low_genes_lt200_cells": int(np.sum(obs["flag_low_genes_lt200"][selected])),
            "flag_high_mt_gt25pct_cells": int(np.sum(obs["flag_high_mt_gt25pct"][selected])),
        }
        donor_rows.append(donor_row)
        donor_detail[donor] = {**donor_row, "published_coarse_annotation_counts": lineage_counts}
    write_table(args.output_dir / "donor_qc.tsv", list(donor_rows[0]), donor_rows)
    summary.update({
        "status": "technical_qc_passed",
        "source_global_metadata_cells": len(all_barcodes),
        "selected_cells_retained": len(barcodes),
        "gene_rows_retained": len(genes),
        "duplicate_source_gene_rows_retained": duplicate_gene_rows,
        "sparse_nonzero_counts": int(matrix.nnz),
        "count_dtype": str(matrix.dtype),
        "sparse_format": "csc",
        "integer_validation": "Every selected count in every gene row parsed as signed int64, nonnegative; all global rows checked for exact header width. Non-selected count values were not numerically inspected.",
        "barcode_join": "Complete count-header and metadata barcode sets are equal; selected cells joined by exact source cell ID, in source count-header order.",
        "reload_checks": "Shape, source gene/barcode axes, all original sparse data/indices, cell totals and donor/source labels match exactly.",
        "sparse_content_sha256": before_hash,
        "qc_flag_definitions": flag_definitions,
        "qc_filtering": "No selected cell, gene row or count removed or normalized. Flags are descriptive, not validated snRNA inclusion rules.",
        "donors": donor_detail,
        "procurement_confounding": confounding,
        "outputs": {
            "h5ad": str(output_h5ad),
            "cell_qc": str(args.output_dir / "cell_qc.tsv"),
            "donor_qc": str(args.output_dir / "donor_qc.tsv"),
        },
    })
    log("complete", selected_cells=len(barcodes), genes=len(genes), donors=len(samples))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename in ["summary.json", "GSE159585_raw_qc.h5ad", "cell_qc.tsv", "donor_qc.tsv"]:
        if (args.output_dir / filename).exists():
            raise FileExistsError(f"Use a fresh output directory; refusing to overwrite {filename}")
    start = time.monotonic()
    summary = {
        "accession": ACCESSION,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "started",
        "full_count_file_read_to_eof": False,
        "gzip_trailer_crc_passed": False,
        "analysis_boundary": "Technical parsing and descriptive QC only; no expression-based biological tests.",
    }
    try:
        run(args, summary)
    except Exception as error:
        summary.update({"status": "failed", "error_type": type(error).__name__, "error": str(error)})
        raise
    finally:
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        summary["peak_rss_bytes"] = int(rss if sys.platform == "darwin" else rss * 1024)
        summary["elapsed_seconds"] = time.monotonic() - start
        summary["completed_utc"] = datetime.now(timezone.utc).isoformat()
        summary["python"] = platform.python_version()
        (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()

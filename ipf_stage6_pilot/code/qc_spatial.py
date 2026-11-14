#!/usr/bin/env python3
"""Parse the six declared Franzen sections without filtering or biological analysis.

X stores the author's raw integer counts for all capture positions, including
off-tissue positions. Author in_tissue and filtered-H5 membership are flags only.
The input ZIP is never extracted wholesale. Reading each selected member to EOF
also invokes zipfile's CRC verification for that member.
"""

import argparse
import hashlib
import io
import json
import re
import traceback
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import anndata as ad
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse


POSITION_COLUMNS = [
    "barcode", "in_tissue", "array_row", "array_col",
    "pxl_row_in_fullres", "pxl_col_in_fullres",
]
REQUIRED_FILES = [
    "raw_feature_bc_matrix.h5", "filtered_feature_bc_matrix.h5",
    "spatial/tissue_positions_list.csv", "spatial/scalefactors_json.json",
    "spatial/tissue_lowres_image.png", "spatial/tissue_hires_image.png",
]


def progress(message):
    print(f"[{datetime.now(timezone.utc).isoformat()}] {message}", flush=True)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def strings(values):
    return np.asarray([
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ], dtype=str)


def read_counts(payload):
    """Read 10x feature-by-barcode CSC and return barcode-by-feature CSR."""
    with h5py.File(io.BytesIO(payload), "r") as handle:
        require("matrix" in handle, "Expected a 10x /matrix HDF5 group")
        group = handle["matrix"]
        shape = tuple(int(value) for value in group["shape"][()])
        values = group["data"][()]
        indices = group["indices"][()]
        indptr = group["indptr"][()]
        barcodes = strings(group["barcodes"][()])
        features = group["features"]
        gene_ids = strings(features["id"][()])
        var = pd.DataFrame({"gene_id": gene_ids}, index=gene_ids)
        for source, target in [
            ("name", "gene_symbol"), ("feature_type", "feature_type"),
            ("genome", "genome"),
        ]:
            if source in features:
                var[target] = strings(features[source][()])
    require(len(shape) == 2, "Count matrix must have two dimensions")
    require(shape == (len(gene_ids), len(barcodes)), "H5 feature/barcode shape mismatch")
    require(len(gene_ids) > 0 and len(barcodes) > 0, "Empty feature or barcode axis")
    require(len(set(barcodes)) == len(barcodes), "Duplicate count barcodes")
    require(len(set(gene_ids)) == len(gene_ids), "Duplicate feature IDs; no automatic renaming")
    require(all(gene_ids) and all(barcodes), "Blank gene IDs or barcodes")
    require(np.isfinite(values).all() and (values >= 0).all(), "Invalid negative/nonfinite counts")
    require(np.equal(values, np.floor(values)).all(), "Noninteger source counts")
    require(len(indptr) == shape[1] + 1 and indptr[0] == 0, "Invalid CSC pointer array")
    require(indptr[-1] == len(values) == len(indices), "CSC arrays have inconsistent lengths")
    require((np.diff(indptr) >= 0).all(), "CSC pointers are not monotonic")
    require(((indices >= 0) & (indices < shape[0])).all(), "CSC feature index outside shape")
    matrix = sparse.csc_matrix((values, indices, indptr), shape=shape).T.tocsr()
    matrix.sort_indices()
    require(matrix.has_canonical_format, "Source sparse matrix has duplicate feature entries")
    return matrix, pd.Index(barcodes, name="barcode"), var


def read_positions(payload):
    frame = pd.read_csv(io.BytesIO(payload), header=None, dtype=str)
    require(frame.shape[1] == 6, "Expected six tissue_positions_list.csv columns")
    if frame.iloc[0, 0] == "barcode":
        frame = frame.iloc[1:].copy()
    frame.columns = POSITION_COLUMNS
    require(not frame.isna().any().any(), "Missing values in source positions")
    require(frame["barcode"].is_unique, "Duplicate position barcodes")
    for column in POSITION_COLUMNS[1:]:
        numeric = pd.to_numeric(frame[column], errors="raise")
        require(np.isfinite(numeric).all(), f"Nonfinite {column}")
        require(np.equal(numeric, np.floor(numeric)).all(), f"Noninteger {column}")
        frame[column] = numeric.astype(np.int64)
    require(frame["in_tissue"].isin([0, 1]).all(), "in_tissue must be author-supplied 0 or 1")
    require((frame[POSITION_COLUMNS[2:]] >= 0).all().all(), "Negative source coordinate")
    return frame.set_index("barcode", drop=False)


def read_image(payload):
    with Image.open(io.BytesIO(payload)) as image:
        image.load()
        require(image.mode in ("RGB", "RGBA"), f"Unexpected tissue image mode: {image.mode}")
        return np.asarray(image).copy()


def selected_metadata(manifest, metadata_path):
    groups = [
        group for group in manifest["primary_data_group"]
        if "E-MTAB-14121" in group.get("accession", "")
    ]
    require(len(groups) == 1, "Manifest must contain one Franzen primary data group")
    sections = groups[0]["selected_sections"]
    require(len(sections) == 6, "This declared pilot expects six sections")
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str, keep_default_na=False)
    require(metadata["sample_id"].is_unique, "Duplicate section IDs in source metadata")
    metadata = metadata.set_index("sample_id", drop=False)
    ids = [section["sample_id"] for section in sections]
    require(len(set(ids)) == len(ids), "Duplicate selected section IDs")
    require(len({section["subject_alias"] for section in sections}) == 6,
            "Declared six-section pilot must retain six distinct source donors")
    for section in sections:
        section_id = section["sample_id"]
        require(re.fullmatch(r"[A-Za-z0-9_-]+", section_id) is not None, "Unsafe section identifier")
        require(section_id in metadata.index, f"Section missing from metadata: {section_id}")
        source = metadata.loc[section_id]
        for key, value in section.items():
            require(key in source.index and str(source[key]) == str(value),
                    f"Manifest/source metadata mismatch: {section_id}/{key}")
    return sections, metadata


def read_member(archive, name, ledger):
    path = PurePosixPath(name)
    require(not path.is_absolute() and ".." not in path.parts, "Unsafe ZIP member name")
    matches = [item for item in archive.infolist() if item.filename == name]
    require(len(matches) == 1, f"Expected one exact ZIP member: {name}")
    info = matches[0]
    require(not info.is_dir(), f"Expected a file ZIP member: {name}")
    with archive.open(info, "r") as handle:
        payload = handle.read()
    require(len(payload) == info.file_size, f"Incomplete ZIP member read: {name}")
    ledger.append({
        "member": name, "uncompressed_bytes": info.file_size,
        "compressed_bytes": info.compress_size, "crc32": f"{info.CRC:08x}",
        "crc_verified_on_read": True, "sha256": hashlib.sha256(payload).hexdigest(),
    })
    return payload


def matrix_equal(left, right):
    return left.shape == right.shape and (left != right).nnz == 0


def quantiles(values):
    return {str(q): float(np.quantile(values, q)) for q in (0, 0.25, 0.5, 0.75, 1)}


def image_checks(positions, images, factors):
    for key in ["tissue_lowres_scalef", "tissue_hires_scalef", "spot_diameter_fullres"]:
        require(key in factors, f"Missing spatial scale: {key}")
        require(np.isfinite(factors[key]) and factors[key] > 0, f"Invalid spatial scale: {key}")
    low, high = images["lowres"], images["hires"]
    low_scale, high_scale = factors["tissue_lowres_scalef"], factors["tissue_hires_scalef"]
    dimension_error = np.abs(np.asarray(low.shape[:2]) - np.asarray(high.shape[:2]) * low_scale / high_scale)
    require((dimension_error <= 1.5).all(), "Low/hires image dimensions disagree with scale factors")
    checks = {
        "lowres_shape": list(low.shape), "hires_shape": list(high.shape),
        "dimension_mapping_error_in_lowres_pixels_y_x": dimension_error.tolist(),
        "dimension_mapping_tolerance_lowres_pixels": 1.5,
    }
    xy = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    in_tissue = positions["in_tissue"].to_numpy() == 1
    require(in_tissue.any(), "No author in-tissue capture positions")
    for resolution, image in images.items():
        scale = factors[f"tissue_{resolution}_scalef"]
        scaled = xy * scale
        inside = ((scaled[:, 0] >= 0) & (scaled[:, 0] < image.shape[1]) &
                  (scaled[:, 1] >= 0) & (scaled[:, 1] < image.shape[0]))
        require(inside[in_tissue].all(), f"Author in-tissue centers fall outside {resolution} image")
        checks[resolution] = {
            "all_position_centers_inside_image": int(inside.sum()),
            "all_position_centers_outside_image": int((~inside).sum()),
            "in_tissue_centers_inside_image": int(inside[in_tissue].sum()),
            "implied_fullres_height_width": [image.shape[0] / scale, image.shape[1] / scale],
            "scaled_spot_diameter_pixels": float(factors["spot_diameter_fullres"] * scale),
        }
    return checks


def overlay(path, positions, images, factors, title):
    figure, axes = plt.subplots(1, 2, figsize=(12, 6))
    flags = positions["in_tissue"].to_numpy() == 1
    for axis, resolution in zip(axes, ("lowres", "hires")):
        image = images[resolution]
        scale = factors[f"tissue_{resolution}_scalef"]
        xy = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy() * scale
        axis.imshow(image, origin="upper")
        axis.scatter(xy[~flags, 0], xy[~flags, 1], s=1, c="deepskyblue", alpha=0.35,
                     linewidths=0, label="Author off-tissue capture position")
        axis.scatter(xy[flags, 0], xy[flags, 1], s=2, c="yellow", alpha=0.45,
                     linewidths=0, label="Author in-tissue capture position")
        axis.set_xlim(-0.5, image.shape[1] - 0.5)
        axis.set_ylim(image.shape[0] - 0.5, -0.5)
        axis.set_title(f"{resolution}: source coordinates x scale")
        axis.axis("off")
    axes[0].legend(loc="lower left", fontsize=6)
    figure.suptitle(title + "\nTechnical registration check; positions are not cells", fontsize=10)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def process_section(archive, section, source_metadata, output):
    section_id = section["sample_id"]
    progress(f"{section_id}: reading selected ZIP members")
    folder = output / section_id
    folder.mkdir()
    members = []
    raw, barcodes, var = read_counts(read_member(archive, f"{section_id}/{REQUIRED_FILES[0]}", members))
    filtered, filtered_barcodes, filtered_var = read_counts(
        read_member(archive, f"{section_id}/{REQUIRED_FILES[1]}", members))
    positions = read_positions(read_member(archive, f"{section_id}/{REQUIRED_FILES[2]}", members))
    factors = json.loads(read_member(archive, f"{section_id}/{REQUIRED_FILES[3]}", members))
    images = {
        "lowres": read_image(read_member(archive, f"{section_id}/{REQUIRED_FILES[4]}", members)),
        "hires": read_image(read_member(archive, f"{section_id}/{REQUIRED_FILES[5]}", members)),
    }
    require(len(positions) == 4992, f"{section_id}: expected 4,992 capture positions")
    require(set(barcodes) == set(positions.index), "Raw H5 and source position barcode sets differ")
    require(var.equals(filtered_var), "Raw/filtered H5 feature identities or order differ")
    indexer = barcodes.get_indexer(filtered_barcodes)
    require((indexer >= 0).all(), "Author filtered barcode missing from raw H5")
    require(matrix_equal(raw[indexer], filtered), "Raw and author filtered counts disagree")
    aligned = positions.loc[barcodes].copy()
    aligned["author_filtered_barcode"] = barcodes.isin(filtered_barcodes)
    filtered_outside = int(((aligned["in_tissue"] == 0) & aligned["author_filtered_barcode"]).sum())
    require(filtered_outside == 0, "Author filtered H5 contains an author off-tissue barcode")
    checks = image_checks(positions, images, factors)
    progress(f"{section_id}: verifying raw counts, position joins and image scale")
    counts = np.asarray(raw.sum(axis=1, dtype=np.int64)).ravel()
    positive = raw > 0
    obs = aligned.copy()
    for key, value in source_metadata.items():
        obs[key] = str(value)
    obs["n_counts"] = counts
    obs["n_features_by_counts"] = np.asarray(positive.sum(axis=1)).ravel().astype(np.int64)
    var["n_capture_positions_by_counts"] = np.asarray(positive.sum(axis=0)).ravel().astype(np.int64)
    var["total_counts"] = np.asarray(raw.sum(axis=0, dtype=np.int64)).ravel()
    if "gene_symbol" in var:
        mito = var["gene_symbol"].str.upper().str.startswith("MT-").to_numpy()
        mito_counts = np.asarray(raw[:, mito].sum(axis=1, dtype=np.int64)).ravel()
        obs["pct_counts_mito"] = np.divide(mito_counts * 100.0, counts,
                                           out=np.full(len(counts), np.nan), where=(counts != 0) & mito.any())
    obs.index = pd.Index([f"{section_id}:{barcode}" for barcode in barcodes], name="capture_position_id")
    data = ad.AnnData(X=raw, obs=obs, var=var)
    data.obsm["spatial"] = obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    data.uns["spatial"] = {section_id: {"images": images, "scalefactors": factors}}
    data.uns["source_metadata"] = {key: str(value) for key, value in source_metadata.items()}
    data.uns["technical_provenance"] = {
        "accession": "E-MTAB-14121 / S-BSST1410",
        "count_source_member": f"{section_id}/raw_feature_bc_matrix.h5",
        "coordinate_order": "obsm spatial is x=fullres column, y=fullres row",
        "data_scope": "All author raw capture positions and all features; no filtering or normalization",
        "unit": "Capture position, not cell; off-tissue positions are retained for technical QC only",
        "author_filtered_counts_verified_against_raw": True,
    }
    target = folder / "raw_counts_spatial.h5ad"
    data.write_h5ad(target, compression="gzip")
    obs.to_csv(folder / "capture_position_qc.tsv", sep="\t", index=True)
    overlay(folder / "coordinate_overlay.png", positions, images, factors,
            f"{section_id} | {section['subject_alias']} | {section['condition']}")
    progress(f"{section_id}: reloading H5AD and checking count/coordinate/image round trip")
    restored = ad.read_h5ad(target)
    require(sparse.issparse(restored.X), "Reloaded X is not sparse")
    require(restored.X.dtype == raw.dtype and matrix_equal(restored.X, raw), "H5AD count round trip failed")
    require(restored.obs_names.equals(data.obs_names) and restored.var_names.equals(data.var_names),
            "H5AD barcode/gene identity round trip failed")
    require(np.array_equal(restored.obsm["spatial"], data.obsm["spatial"]), "Coordinate round trip failed")
    for column in ["barcode", "in_tissue", "author_filtered_barcode", "subject_alias", "condition"]:
        require(np.array_equal(restored.obs[column].to_numpy(), data.obs[column].to_numpy()),
                f"Source flag/identity round trip failed: {column}")
    for resolution in images:
        require(np.array_equal(restored.uns["spatial"][section_id]["images"][resolution], images[resolution]),
                f"Image round trip failed: {resolution}")
    for key, value in factors.items():
        require(restored.uns["spatial"][section_id]["scalefactors"][key] == value,
                f"Scale factor round trip failed: {key}")
    in_tissue = obs["in_tissue"].to_numpy() == 1
    warnings = []
    if int(in_tissue.sum()) != len(filtered_barcodes):
        warnings.append("Author in_tissue and filtered-H5 barcode sets differ; both flags and all raw counts retained")
    summary = {
        "section_id": section_id, "donor_id": section["subject_alias"],
        "diagnosis": section["condition"], "source_metadata": dict(source_metadata),
        "capture_positions": data.n_obs, "features": data.n_vars,
        "author_in_tissue_positions": int(in_tissue.sum()),
        "author_off_tissue_positions": int((~in_tissue).sum()),
        "author_filtered_h5_positions": len(filtered_barcodes),
        "all_position_and_raw_barcode_sets_equal": True,
        "filtered_count_submatrix_equal": True,
        "total_umis_all_capture_positions": int(counts.sum()),
        "total_umis_author_in_tissue": int(counts[in_tissue].sum()),
        "zero_count_positions": int((counts == 0).sum()),
        "counts_per_author_in_tissue_position_quantiles": quantiles(counts[in_tissue]),
        "features_per_author_in_tissue_position_quantiles": quantiles(obs.loc[in_tissue, "n_features_by_counts"]),
        "unique_gene_ids": bool(var.index.is_unique),
        "ensembl_gene_id_count": sum(bool(re.fullmatch(r"ENSG\d+(?:\.\d+)?", value)) for value in var.index),
        "duplicate_gene_symbol_count": int(var["gene_symbol"].duplicated().sum()) if "gene_symbol" in var else 0,
        "feature_type_counts": var["feature_type"].value_counts().to_dict() if "feature_type" in var else {},
        "image_checks": checks, "selected_member_checks": members,
        "reload_verified": True, "h5ad": str(target),
        "coordinate_overlay": str(folder / "coordinate_overlay.png"),
        "capture_position_qc": str(folder / "capture_position_qc.tsv"),
        "warnings": warnings, "status": "passed_technical_checks",
    }
    write_json(folder / "section_summary.json", summary)
    progress(f"{section_id}: passed; {data.n_obs} capture positions, {int(in_tissue.sum())} author in-tissue")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    require(not any(args.output_dir.iterdir()), "Output directory must be empty; use a unique run directory")
    archive_path = args.raw_dir / "hs_visium_spaceranger_output.zip"
    metadata_path = args.raw_dir / "hs_visium_metadata.tsv"
    result = {
        "status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "accession": "E-MTAB-14121 / S-BSST1410", "sections": [],
        "archive": str(archive_path), "manifest": str(args.manifest),
        "count_policy": "Raw H5 counts for all 4,992 capture positions; all features retained; no normalization or filtering",
        "analysis_boundary": "Technical parsing, donor/source joins and image-coordinate agreement only; no deconvolution, spatial association or target inference",
    }
    try:
        manifest = json.loads(args.manifest.read_text())
        sections, metadata = selected_metadata(manifest, metadata_path)
        result["archive_bytes"] = archive_path.stat().st_size
        result["metadata_sha256"] = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
        with zipfile.ZipFile(archive_path) as archive:
            for section in sections:
                row = process_section(archive, section, metadata.loc[section["sample_id"]].to_dict(), args.output_dir)
                result["sections"].append(row)
                write_json(args.output_dir / "summary.json", result)
        flat = [{key: row[key] for key in [
            "section_id", "donor_id", "diagnosis", "capture_positions", "features",
            "author_in_tissue_positions", "author_off_tissue_positions", "author_filtered_h5_positions",
            "total_umis_all_capture_positions", "total_umis_author_in_tissue", "zero_count_positions",
            "reload_verified", "status",
        ]} for row in result["sections"]]
        qc = pd.DataFrame(flat)
        qc.to_csv(args.output_dir / "section_donor_qc.tsv", sep="\t", index=False)
        result["section_donor_qc"] = str(args.output_dir / "section_donor_qc.tsv")
        result["biological_donors"] = int(qc["donor_id"].nunique())
        result["sections_completed"] = len(qc)
        result["status"] = "passed_technical_checks"
        result["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(args.output_dir / "summary.json", result)
        progress("Franzen six-section technical QC completed")
    except Exception as error:
        result["status"] = "failed"
        result["error"] = str(error)
        result["traceback"] = traceback.format_exc()
        result["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(args.output_dir / "summary.json", result)
        raise


if __name__ == "__main__":
    main()

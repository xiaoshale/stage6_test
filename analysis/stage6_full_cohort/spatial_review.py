#!/usr/bin/env python3
"""Full Franzén technical QC and candidate-independent human histology review pack.

Reuses the validated pilot parser; no normalization, ROI selection by RNA,
deconvolution, anatomical certification or disease-effect analysis is performed.
"""

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import re
import traceback
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

from PIL import Image, ImageDraw, ImageFont


EXPECTED_DONORS = {
    "HC_1": 2, "HC_2": 1, "HC_3": 2, "HC_4": 1,
    "IPF_1": 6, "IPF_2": 6, "IPF_3": 3, "IPF_4": 4,
}
REVIEW_FIELDS = [
    "section_id", "donor_id", "diagnosis", "block_id", "region_label", "replicate",
    "procurement", "technical_status", "source_image", "source_image_sha256",
    "source_image_width", "source_image_height", "review_image", "image_scale",
    "histology_mapping_status", "histology_dimension_check", "image_sufficiency",
    "compartment", "interface_status", "mask_status", "expert_status", "reviewer",
    "reviewer_role", "review_date", "molecular_exposure", "endpoint_eligibility", "notes",
]
ROI_FIELDS = [
    "roi_id", "section_id", "donor_id", "block_id", "coordinate_frame", "roi_polygon",
    "artifact_mask", "barrier_mask", "tissue_component_id", "interface_geometry",
    "interface_extent", "boundary_uncertainty", "source_annotation", "project_compartment",
    "interface_type", "interface_status", "reviewer", "reviewer_role", "review_date",
    "expert_status", "molecular_exposure", "endpoint_id", "endpoint_eligibility", "reason",
]
FEATURE_FIELDS = [
    "source_row_1based", "source_gene_id", "source_symbol", "genome", "feature_type",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export_feature_axis(features, path):
    """Preserve the native row order using the common gene-mapper input schema."""
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(FEATURE_FIELDS)
        count = 0
        for count, (_, feature) in enumerate(features.iterrows(), 1):
            writer.writerow([count, feature["gene_id"], feature.get("gene_symbol", ""),
                             feature.get("genome", ""), feature.get("feature_type", "")])
    return {"path": str(path), "sha256": sha256(path), "source_features": count,
            "schema": FEATURE_FIELDS, "source_row_order_preserved": True}


def load_qc(path):
    spec = importlib.util.spec_from_file_location("stage6_pilot_spatial_qc", path)
    require(spec is not None and spec.loader is not None, "Cannot load the pilot spatial module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metadata_rows(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(len(rows) == 25, "Expected the native 25-section acquisition")
    require(len({r["sample_id"] for r in rows}) == 25, "Duplicate section IDs")
    require(Counter(r["subject_alias"] for r in rows) == EXPECTED_DONORS,
            "Native donor/section counts differ from the fixed intake")
    require(len({r["tissue_alias"] for r in rows}) == 16, "Expected 16 source tissue blocks")
    require(Counter(r["condition"] for r in rows) == {"IPF": 19, "control": 6},
            "Native diagnosis/section counts differ")
    for row in rows:
        require(re.fullmatch(r"[A-Za-z0-9_-]+", row["sample_id"]) is not None,
                "Unsafe section identifier")
        require((row["condition"] == "IPF") == row["subject_alias"].startswith("IPF_"),
                "Source donor diagnosis mismatch")
    return rows


def histology_members(archive, rows):
    """Map exact section/sample-name tokens only; never guess by donor or block."""
    images = [i for i in archive.infolist()
              if not i.is_dir() and "__MACOSX" not in PurePosixPath(i.filename).parts
              and not PurePosixPath(i.filename).name.startswith("._")
              and PurePosixPath(i.filename).suffix.lower() in (".jpg", ".jpeg", ".png", ".tif", ".tiff")]
    result = {}
    for row in rows:
        candidates = []
        for info in images:
            name = str(PurePosixPath(info.filename))
            if any(re.search(r"(?<![A-Za-z0-9])" + re.escape(token) + r"(?![A-Za-z0-9])", name)
                   for token in (row["sample_id"], row["sample_name"])):
                candidates.append(info.filename)
        result[row["sample_id"]] = {
            "status": "unique_exact_token" if len(candidates) == 1 else
                      "missing_exact_token" if not candidates else "ambiguous_exact_token",
            "candidates": candidates,
        }
    chosen = [r["candidates"][0] for r in result.values() if r["status"] == "unique_exact_token"]
    require(len(chosen) == len(set(chosen)), "One histology image matched multiple sections")
    return result


def extract_image(archive, member, target):
    path = PurePosixPath(member)
    require(not path.is_absolute() and ".." not in path.parts, "Unsafe histology ZIP member")
    infos = [i for i in archive.infolist() if i.filename == member]
    require(len(infos) == 1, "Duplicate histology ZIP member")
    info = infos[0]
    digest = hashlib.sha256()
    with archive.open(info) as source, target.open("wb") as destination:
        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
            destination.write(chunk)
    require(target.stat().st_size == info.file_size, "Incomplete original histology extraction")
    return {"member": member, "bytes": info.file_size, "crc32": f"{info.CRC:08x}",
            "crc_verified_on_read": True, "sha256": digest.hexdigest(),
            "original_bytes_preserved": True, "server_only_original_path": str(target)}


def font(size):
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size)
    except OSError:
        return ImageFont.load_default()


def open_histology(path):
    # Bound the decoded pixels explicitly; original JPEG bytes remain unchanged.
    old_limit = Image.MAX_IMAGE_PIXELS
    try:
        Image.MAX_IMAGE_PIXELS = None
        with Image.open(path) as source:
            require(source.width * source.height <= 300_000_000,
                    "Histology exceeds the declared 300-million-pixel review decode limit")
            orientation = int(source.getexif().get(274, 1))
            image = source.convert("RGB")
    finally:
        Image.MAX_IMAGE_PIXELS = old_limit
    return image, orientation


def review_images(image, folder, title, provenance, side=2500):
    """Whole field and a fixed 3x3 grid; neither uses RNA nor a predicted ROI."""
    preview = image.copy()
    preview.thumbnail((side, side))
    scale_x, scale_y = preview.width / image.width, preview.height / image.height
    canvas = Image.new("RGB", (preview.width, preview.height + 100), "white")
    canvas.paste(preview, (0, 100))
    draw = ImageDraw.Draw(canvas)
    draw.text((12, 8), title, font=font(24), fill="black")
    draw.text((12, 42), provenance, font=font(18), fill="black")
    draw.text((12, 70), "H&E review pending; no anatomical label assigned", font=font(18), fill="black")
    canvas.save(folder / "tissue_overview.jpg", quality=90)
    tiles = []
    grid = Image.new("RGB", (3000, 3120), "white")
    for row in range(3):
        for col in range(3):
            box = (col * image.width // 3, row * image.height // 3,
                   (col + 1) * image.width // 3, (row + 1) * image.height // 3)
            tile = image.crop(box)
            tile.thumbnail((1000, 960))
            origin = (col * 1000, row * 1040)
            grid.paste(tile, (origin[0], origin[1] + 60))
            ImageDraw.Draw(grid).text(origin, f"R{row + 1}C{col + 1} source px {box}",
                                      font=font(18), fill="black")
            tiles.append({"row": row + 1, "column": col + 1, "source_pixel_box_xyxy": list(box),
                          "display_width_height": list(tile.size)})
    grid.save(folder / "fixed_grid_zoom.jpg", quality=90)
    write_json(folder / "review_image_provenance.json", {
        "title": title, "source": provenance, "source_width_height": list(image.size),
        "overview_content_width_height": list(preview.size),
        "overview_content_offset_xy": [0, 100], "overview_scale_xy": [scale_x, scale_y],
        "zoom_grid": tiles, "selection": "Uniform 3x3 image grid, independent of measured RNA",
        "expert_review_status": "unreviewed", "not_a_cell_segmentation": True,
    })
    return [scale_x, scale_y]


def diagnostic_plots(qc, data, folder, row, limits):
    import numpy as np
    import matplotlib.pyplot as plt

    sid = row["sample_id"]
    image = data.uns["spatial"][sid]["images"]["hires"]
    factor = float(data.uns["spatial"][sid]["scalefactors"]["tissue_hires_scalef"])
    xy = data.obsm["spatial"] * factor
    specifications = [
        ("raw_umi.png", np.log1p(data.obs["n_counts"].to_numpy()), "log1p raw UMI", limits[0]),
        ("detected_genes.png", data.obs["n_features_by_counts"].to_numpy(), "Detected source features", limits[1]),
        ("author_in_tissue_mask.png", data.obs["in_tissue"].to_numpy(), "Author in_tissue flag (0/1)", 1),
    ]
    for filename, values, label, upper in specifications:
        figure, axis = plt.subplots(figsize=(8, 8))
        axis.imshow(image, origin="upper", alpha=0.65)
        points = axis.scatter(xy[:, 0], xy[:, 1], c=values, s=5, linewidths=0,
                              cmap="viridis", vmin=0, vmax=max(float(upper), 1), alpha=0.85)
        axis.set_xlim(-0.5, image.shape[1] - 0.5)
        axis.set_ylim(image.shape[0] - 0.5, -0.5)
        axis.axis("off")
        axis.set_title(f"{sid} | {row['subject_alias']} | {row['condition']}\n"
                       "All raw capture positions; not cells or anatomical certification", fontsize=10)
        figure.colorbar(points, ax=axis, shrink=0.7, label=label + " (shared scale)")
        figure.tight_layout()
        figure.savefig(folder / filename, dpi=140)
        plt.close(figure)


def contact_sheet(rows, output):
    sheet = Image.new("RGB", (2500, 2900), "white")
    for index, row in enumerate(rows):
        with Image.open(output / row["sample_id"] / "tissue_overview.jpg") as image:
            image.thumbnail((490, 490))
            x, y = (index % 5) * 500, (index // 5) * 580
            sheet.paste(image, (x, y))
        ImageDraw.Draw(sheet).text((x + 5, y + 495),
                                  f"{row['sample_id']}\n{row['subject_alias']} {row['B_tissue_selection']} rep{row['replicate']}\nExpert: unreviewed",
                                  font=font(17), fill="black")
    sheet.save(output / "human_review" / "all25_histology_contact_sheet.jpg", quality=90)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--histology-archive", type=Path)
    parser.add_argument("--pilot-qc-module", type=Path)
    args = parser.parse_args()
    candidates = [Path(__file__).with_name("qc_spatial.py"),
                  Path(__file__).resolve().parents[1] / "stage6_pilot" / "qc_spatial.py"]
    qc_path = args.pilot_qc_module or next((p for p in candidates if p.exists()), candidates[0])
    require(qc_path.is_file(), "Supply --pilot-qc-module pointing to the unchanged pilot qc_spatial.py")
    qc = load_qc(qc_path)
    rows = metadata_rows(args.metadata)
    args.output.mkdir(parents=True, exist_ok=True)
    require(not any(args.output.iterdir()), "Output must be empty; use a unique run directory")
    review_dir = args.output / "human_review"
    axes_dir = args.output / "feature_axes"
    review_dir.mkdir()
    axes_dir.mkdir()
    result = {
        "status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "resource": "E-MTAB-14121 / S-BSST1410", "sections": [],
        "candidate_independent": True, "expert_review_status": "unreviewed",
        "formal_anatomical_eligibility": "pending_human_review",
        "procurement_limit": "HC post-mortem versus IPF transplant/resection; no disease-only effect identifiable",
        "count_policy": "Original raw sparse counts, all positions/features; author flags retained; no filtering or normalization",
        "sources": [{"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)}
                    for p in [args.archive, args.metadata, qc_path, Path(__file__)]] +
                   ([{"path": str(args.histology_archive), "bytes": args.histology_archive.stat().st_size,
                      "sha256": sha256(args.histology_archive)}] if args.histology_archive else []),
    }
    histology = zipfile.ZipFile(args.histology_archive) if args.histology_archive else None
    try:
        mapping = histology_members(histology, rows) if histology else {
            r["sample_id"]: {"status": "original_histology_not_supplied", "candidates": []} for r in rows}
        write_json(review_dir / "histology_member_mapping.json", mapping)
        with zipfile.ZipFile(args.archive) as archive:
            require(archive.testzip() is None, "Space Ranger archive CRC check failed")
            result["whole_spaceranger_archive_crc_passed"] = True
            for row in rows:
                sid = row["sample_id"]
                section = qc.process_section(archive, row, row, args.output)
                result["sections"].append(section)
                write_json(args.output / "summary.json", result)
        import numpy as np

        limits = (max(math.log1p(s["counts_per_author_in_tissue_position_quantiles"]["1"]) for s in result["sections"]),
                  max(s["features_per_author_in_tissue_position_quantiles"]["1"] for s in result["sections"]))
        review_rows = []
        for row, summary in zip(rows, result["sections"]):
            sid = row["sample_id"]
            folder = args.output / sid
            data = qc.ad.read_h5ad(folder / "raw_counts_spatial.h5ad")
            # Original feature rows are exported after count and H5AD roundtrip checks.
            summary["feature_axis"] = export_feature_axis(
                data.var, axes_dir / f"{sid}.features.tsv.gz")
            require(summary["feature_axis"]["source_features"] == data.n_vars,
                    "Feature export does not preserve every native row")
            diagnostic_plots(qc, data, folder, row, limits)
            native = data.uns["spatial"][sid]
            image = Image.fromarray(np.asarray(native["images"]["hires"])).convert("RGB")
            provenance = "Space Ranger hires fallback; original H&E pending"
            original = None
            dimension_status = "not_assessable_original_not_supplied_or_unmapped"
            if mapping[sid]["status"] == "unique_exact_token":
                member = mapping[sid]["candidates"][0]
                original_dir = folder / "original_histology"
                original_dir.mkdir()
                target = original_dir / (sid + PurePosixPath(member).suffix.lower())
                original = extract_image(histology, member, target)
                image, orientation = open_histology(target)
                original["width_height"] = list(image.size)
                original["exif_orientation"] = orientation
                scale = float(native["scalefactors"]["tissue_hires_scalef"])
                expected = np.asarray(native["images"]["hires"].shape[:2][::-1]) / scale
                error = np.abs(np.asarray(image.size) - expected)
                dimension_status = "dimensions_consistent_visual_review_pending" if (
                    (error <= 1.5 / scale).all() and orientation == 1) else "hold_dimensions_or_exif_orientation"
                original["dimension_error_fullres_pixels_xy"] = error.tolist()
                original["dimension_status"] = dimension_status
                provenance = f"Original H&E {image.width}x{image.height}px; {dimension_status}"
            scale_xy = review_images(image, folder,
                                     f"{sid} | {row['subject_alias']} | {row['condition']}", provenance)
            summary["histology"] = {"mapping": mapping[sid], "original": original,
                                    "dimension_status": dimension_status, "expert_status": "unreviewed"}
            write_json(folder / "section_summary.json", summary)
            write_json(folder / "anatomy_mask_template.json", {
                "section_id": sid, "coordinate_frame": "original Space Ranger input image pixels x/y",
                "source_image_dimensions_confirmed": dimension_status.startswith("dimensions_consistent"),
                "artifact_polygons": [], "barrier_polygons": [], "compartment_polygons": [],
                "interface_polylines": [], "review_status": "unreviewed",
                "empty_arrays_mean": "No annotations have been made; not absence of artifacts or interfaces",
                "endpoint_eligibility": "not_adjudicated",
            })
            record = {key: "" for key in REVIEW_FIELDS}
            record.update({
                "section_id": sid, "donor_id": row["subject_alias"], "diagnosis": row["condition"],
                "block_id": row["tissue_alias"], "region_label": row["B_tissue_selection"], "replicate": row["replicate"],
                "procurement": "post_mortem" if row["condition"] == "control" else "transplant_or_resection",
                "technical_status": summary["status"],
                "source_image": original["server_only_original_path"] if original else "Space Ranger hires in H5AD",
                "source_image_sha256": original["sha256"] if original else next(
                    m["sha256"] for m in summary["selected_member_checks"] if m["member"].endswith("tissue_hires_image.png")),
                "source_image_width": image.width, "source_image_height": image.height,
                "review_image": str(folder / "tissue_overview.jpg"), "image_scale": json.dumps(scale_xy),
                "histology_mapping_status": mapping[sid]["status"], "histology_dimension_check": dimension_status,
                "image_sufficiency": "unreviewed", "compartment": "unreviewed", "interface_status": "not_assessed",
                "mask_status": "unreviewed_empty_template", "expert_status": "unreviewed",
                "molecular_exposure": "Technical UMI/features available; no candidate-expression ROI selection",
                "endpoint_eligibility": "not_adjudicated",
            })
            review_rows.append(record)
            qc.progress(f"{sid}: candidate-independent image review material generated; expert status unreviewed")
            del data, image
        with (review_dir / "section_review.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=REVIEW_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(review_rows)
        with (review_dir / "roi_inventory.tsv").open("w") as handle:
            csv.writer(handle, delimiter="\t").writerow(ROI_FIELDS)
        contact_sheet(rows, args.output)
        (review_dir / "README.md").write_text(
            "# Human spatial anatomy review queue\n\n"
            "All 25 sections are included. Every expert/anatomical/endpoint decision remains unreviewed.\n"
            "Review H&E first; inventory all target-compatible and comparison regions independently of candidate RNA.\n"
            "The fixed zoom grid covers the whole image and does not identify an airway or a niche.\n"
            "Original histology bytes stay in server-only original_histology folders. Review JPEGs are resized derivatives.\n"
            "Check source dimensions, registration, basal-side interfaces, lumens, tears and folds before marking regions.\n"
            "Empty mask arrays mean unreviewed, not artifact-free. Do not connect tissue across cavities or disconnected sections.\n"
            "Visium positions are mixtures, not cells. Technical UMI/feature plots do not certify anatomy.\n"
            "HC is post-mortem and IPF transplant/resection; matched anatomy cannot identify a disease-only effect.\n"
            "Use donor/block/section/ROI nesting; prioritize within-IPF regional questions.\n\n"
            "Source: https://www.ebi.ac.uk/biostudies/studies/S-BSST1410\n"
            "Image semantics: https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs\n")
        result.update({
            "status": "passed_technical_checks_human_review_pending", "sections_completed": len(rows),
            "biological_donors": len(EXPECTED_DONORS), "native_features_exported": True,
            "original_histology_sections": sum(r["histology_mapping_status"] == "unique_exact_token" for r in review_rows),
            "original_histology_dimension_holds": sum(r["histology_dimension_check"].startswith("hold") for r in review_rows),
            "technical_plot_scale": {"log1p_raw_umi_max": limits[0], "detected_features_max": limits[1],
                                     "basis": "Shared maxima across author in-tissue positions; off-tissue values may saturate"},
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        })
        write_json(args.output / "summary.json", result)
        artifacts = [{"path": str(p.relative_to(args.output)), "bytes": p.stat().st_size, "sha256": sha256(p)}
                     for p in sorted(args.output.rglob("*")) if p.is_file()]
        write_json(args.output / "artifact_sha256.json", {"created_after_count_checks": True, "artifacts": artifacts})
        qc.progress("All 25 native sections passed technical QC; human anatomy/eligibility review remains pending")
    except Exception as error:
        result.update({"status": "failed", "error": str(error), "traceback": traceback.format_exc()})
        write_json(args.output / "summary.json", result)
        raise
    finally:
        if histology is not None:
            histology.close()


if __name__ == "__main__":
    main()

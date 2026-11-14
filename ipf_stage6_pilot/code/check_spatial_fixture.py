#!/usr/bin/env python3
"""Exercise spatial parsing and failure handling with a small synthetic ZIP.

The fixture keeps 4,992 capture positions to exercise the declared assay guard,
but has only three features, four nonzero positions and very small images.
Nothing is downloaded or deleted. Use a fresh --work-dir for every invocation.
"""

import argparse
import csv
import io
import json
import subprocess
import sys
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from PIL import Image
from scipy import sparse


def h5_bytes(matrix, barcodes):
    counts = sparse.csc_matrix(matrix.T)
    buffer = io.BytesIO()
    with h5py.File(buffer, "w") as handle:
        group = handle.create_group("matrix")
        for key, values in {
            "data": counts.data, "indices": counts.indices,
            "indptr": counts.indptr, "shape": counts.shape,
            "barcodes": np.asarray(barcodes, dtype="S"),
        }.items():
            group.create_dataset(key, data=values)
        features = group.create_group("features")
        for key, values in {
            "id": ["ENSG00000000001", "ENSG00000000002", "ENSG00000000003"],
            "name": ["KRT17", "COL1A1", "MT-ND1"],
            "feature_type": ["Gene Expression"] * 3,
            "genome": ["GRCh38"] * 3,
        }.items():
            features.create_dataset(key, data=np.asarray(values, dtype="S"))
    return buffer.getvalue()


def png_bytes(image):
    buffer = io.BytesIO()
    image.save(buffer, format="PNG")
    return buffer.getvalue()


def build_inputs(folder, bad_filtered=False):
    folder.mkdir()
    barcodes = [f"BC{index:05d}-1" for index in range(4992)]
    raw = np.zeros((4992, 3), dtype=np.int32)
    raw[0] = [1, 0, 0]
    raw[1] = [0, 2, 3]
    raw[17] = [4, 5, 0]
    raw[-1] = [0, 0, 11]  # Off-tissue signal must be retained.
    selected = [0, 1, 17]
    filtered = raw[selected].copy()
    if bad_filtered:
        filtered[0, 0] = 99
    positions_buffer = io.StringIO()
    writer = csv.writer(positions_buffer, lineterminator="\n")
    for index, barcode in enumerate(barcodes):
        row, column = divmod(index, 78)
        writer.writerow([barcode, int(index in selected), row, column,
                         20 + row * 3, 10 + column * 2])
    source_image = np.zeros((250, 200, 3), dtype=np.uint8)
    source_image[:, :, 0] = np.arange(200)[None, :]
    source_image[:, :, 1] = np.arange(250)[:, None]
    source_image[:, :, 2] = 100
    image = Image.fromarray(source_image)
    low = image.resize((40, 50))
    high = image.resize((100, 125))
    factors = {"tissue_lowres_scalef": 0.2, "tissue_hires_scalef": 0.5,
               "spot_diameter_fullres": 4.0, "fiducial_diameter_fullres": 6.0}
    sections = [
        {"sample_id": f"FIXTURE-{index}", "sample_name": f"fixture_{index}",
         "subject_alias": f"D{index}", "condition": "IPF" if index <= 3 else "control",
         "tissue": "lung"}
        for index in range(1, 7)
    ]
    with (folder / "hs_visium_metadata.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sections[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(sections)
    manifest = folder / "pilot_manifest.json"
    manifest.write_text(json.dumps({"primary_data_group": [{
        "accession": "E-MTAB-14121 / S-BSST1410", "selected_sections": sections,
    }]}, indent=2) + "\n")
    payloads = {
        "raw_feature_bc_matrix.h5": h5_bytes(raw, barcodes),
        "filtered_feature_bc_matrix.h5": h5_bytes(filtered, [barcodes[index] for index in selected]),
        "spatial/tissue_positions_list.csv": positions_buffer.getvalue().encode(),
        "spatial/scalefactors_json.json": json.dumps(factors).encode(),
        "spatial/tissue_lowres_image.png": png_bytes(low),
        "spatial/tissue_hires_image.png": png_bytes(high),
    }
    with zipfile.ZipFile(folder / "hs_visium_spaceranger_output.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for section in sections:
            for name, payload in payloads.items():
                archive.writestr(f"{section['sample_id']}/{name}", payload)
    return manifest, raw


def run_case(script, raw_dir, manifest, output, log_prefix):
    run = subprocess.run([
        sys.executable, str(script), "--raw-dir", str(raw_dir),
        "--manifest", str(manifest), "--output-dir", str(output),
    ], capture_output=True, text=True)
    log_prefix.with_suffix(".stdout.log").write_text(run.stdout)
    log_prefix.with_suffix(".stderr.log").write_text(run.stderr)
    return run


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    check(not any(args.work_dir.iterdir()), "Fixture work directory must be empty")
    script = Path(__file__).with_name("qc_spatial.py")
    print("Building small six-section fixture with asymmetric x/y coordinates", flush=True)
    good_dir = args.work_dir / "good_inputs"
    manifest, expected = build_inputs(good_dir)
    good_output = args.work_dir / "good_output"
    good = run_case(script, good_dir, manifest, good_output, args.work_dir / "good")
    check(good.returncode == 0, "Good fixture failed: " + good.stderr[-2000:])
    summary = json.loads((good_output / "summary.json").read_text())
    check(summary["status"] == "passed_technical_checks", "Good fixture status is not passed")
    check(summary["sections_completed"] == 6 and summary["biological_donors"] == 6,
          "Section/donor totals differ from the fixture")
    for index in range(1, 7):
        section_id = f"FIXTURE-{index}"
        data = ad.read_h5ad(good_output / section_id / "raw_counts_spatial.h5ad")
        check(sparse.issparse(data.X), "Stored counts became dense")
        check(data.shape == (4992, 3), "A raw spot or feature was dropped")
        check(np.array_equal(data.X.toarray(), expected), "Stored raw counts differ from hand-set counts")
        check(np.array_equal(data.obsm["spatial"][[0, 1, 17, -1]],
                             [[10, 20], [12, 20], [44, 20], [164, 209]]),
              "Stored x/y orientation or barcode ordering is incorrect")
        check(data.obs["in_tissue"].to_numpy().tolist().count(1) == 3, "Author tissue flags changed")
        check(np.flatnonzero(data.obs["author_filtered_barcode"]).tolist() == [0, 1, 17],
              "Author filtered flags changed")
        check(data.X[-1].sum() == 11 and data.obs.iloc[-1]["in_tissue"] == 0,
              "Off-tissue counts were filtered")
        check(data.X[2].sum() == 0, "Zero-count capture position was not preserved")
        check(data.uns["spatial"][section_id]["images"]["lowres"].shape == (50, 40, 3),
              "Low-resolution image shape changed")
        check(data.uns["spatial"][section_id]["scalefactors"]["tissue_hires_scalef"] == 0.5,
              "Image scale changed")
    for section in summary["sections"]:
        data = ad.read_h5ad(good_output / section["section_id"] / "raw_counts_spatial.h5ad")
        assert data.obs.loc[data.obs["n_counts"] == 0, "pct_counts_mito"].isna().all()
    print("Good fixture passed: raw counts, flags, x/y, images, scales and H5AD reload", flush=True)
    bad_dir = args.work_dir / "bad_filtered_inputs"
    bad_manifest, _ = build_inputs(bad_dir, bad_filtered=True)
    bad_output = args.work_dir / "bad_filtered_output"
    bad = run_case(script, bad_dir, bad_manifest, bad_output, args.work_dir / "bad_filtered")
    check(bad.returncode != 0, "Inconsistent filtered counts were incorrectly accepted")
    failed = json.loads((bad_output / "summary.json").read_text())
    check(failed["status"] == "failed" and "Raw and author filtered counts disagree" in failed["error"],
          "Bad fixture failed for an unexpected reason")
    report = {
        "status": "passed", "good_fixture_six_section_round_trip": True,
        "asymmetric_xy_verified": True, "raw_off_tissue_and_zero_count_positions_preserved": True,
        "bad_filtered_counts_rejected": True, "good_summary": str(good_output / "summary.json"),
        "expected_failure_summary": str(bad_output / "summary.json"),
        "scope": "Synthetic technical fixture only; no evidence about the real dataset or biological interpretation",
    }
    (args.work_dir / "fixture_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Bad filtered-count fixture rejected as expected; fixture_checks.json written", flush=True)


if __name__ == "__main__":
    main()

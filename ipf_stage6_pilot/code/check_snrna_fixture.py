"""Exercise the snRNA CLI with an independently specified tiny source matrix.

Run in the existing scientific Python environment. The deliberately dense
expected matrix below is only 5 by 3; the production reader remains sparse.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np


SCRIPT = Path(__file__).with_name("qc_snrna.py")
COUNTS = "GSE159585_snRNAseq_raw_counts.csv.gz"
METADATA = "GSE159585_snRNAseq_metadata.csv.gz"

CELL_A = "LUNG04_IPF_CELL_A_1"
CELL_ZERO = "LUNG04_IPF_CELL_ZERO_1"
CELL_B = "LUNG03_IPF_CELL_B_1"
CELL_C = "LUNG01_HC_CELL_C_1"
CELL_D = "LUNG010_HC_CELL_D_1"
CELL_COVID = "LUNG19_COVID-19_CELL_E_1"
HEADER = ["", CELL_COVID, CELL_D, CELL_ZERO, CELL_C, CELL_B, CELL_A]
GENE_ROWS = [
    ["GENE_A", 9, 1, 0, 0, 2, 4],
    ["MT-CO1", 9, 0, 0, 5, 1, 0],
    ["GENE_A", 9, 3, 0, 0, 0, 2],
]


def create_case(directory: Path, failure: str | None = None) -> tuple[Path, Path, Path]:
    directory.mkdir(parents=True, exist_ok=False)
    raw = directory / "raw"
    raw.mkdir()
    metadata_rows = [
        [CELL_A, "LUNG04", "IPF1", "IPF", "epithelial", "source_A"],
        [CELL_B, "LUNG03", "IPF2", "IPF", "stroma", "source_B"],
        [CELL_COVID, "LUNG19", "COV-1", "COVID-19", "immune", "excluded"],
        [CELL_C, "LUNG01", "Non-COV1", "control", "endothelial", "source_C"],
        [CELL_D, "LUNG010", "Non-COV2", "control", "immune", "source_D"],
        [CELL_ZERO, "LUNG04", "IPF1", "IPF", "epithelial", "source_zero"],
    ]
    if failure == "missing_metadata_cell":
        # Keep every selected donor present, so this specifically exercises
        # the complete global-header join, including the unselected COVID cell.
        metadata_rows = [row for row in metadata_rows if row[0] != CELL_COVID]
    with gzip.open(raw / METADATA, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter=";")
        writer.writerow(["", "sample_ID", "patient_ID", "condition", "cellType", "subcluster"])
        writer.writerows(metadata_rows)
    gene_rows = [list(row) for row in GENE_ROWS]
    if failure == "fractional_count":
        gene_rows[0][2] = "0.5"
    with gzip.open(raw / COUNTS, "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerows(gene_rows)
    selected = []
    for donor, source_sample, condition, procurement, gsm in [
        ("IPF1", "LUNG04", "IPF", "explant", "GSM_TEST_IPF1"),
        ("IPF2", "LUNG03", "IPF", "explant", "GSM_TEST_IPF2"),
        ("Non-COV1", "LUNG01", "control", "autopsy", "GSM_TEST_C1"),
        ("Non-COV2", "LUNG010", "control", "autopsy", "GSM_TEST_C2"),
    ]:
        selected.append({
            "patient_ID": donor,
            "condition": condition,
            "metadata_sample_IDs": [source_sample],
            "GSM": gsm,
            "GEO_characteristics": [f"sample type: {procurement}"],
        })
    manifest = directory / "manifest.json"
    manifest.write_text(json.dumps({
        "technical_robustness_component": {
            "accession": "GSE159585", "selected_samples": selected,
        },
    }))
    return raw, manifest, directory / "out"


def invoke(raw: Path, manifest: Path, output: Path):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--raw-dir", str(raw), "--manifest", str(manifest), "--output-dir", str(output)],
        capture_output=True,
        text=True,
        check=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    raw, manifest, output = create_case(args.work_dir / "valid")
    result = invoke(raw, manifest, output)
    if result.returncode:
        raise AssertionError(f"Valid fixture failed:\n{result.stdout}\n{result.stderr}")
    data = ad.read_h5ad(output / "GSE159585_raw_qc.h5ad")
    expected = np.array([
        [1, 0, 3],  # Non-COV2, exact LUNG010 source ID retained.
        [0, 0, 0],  # Zero-count IPF1 cell is flagged, not removed.
        [0, 5, 0],  # Non-COV1.
        [2, 1, 0],  # IPF2.
        [4, 0, 2],  # IPF1.
    ], dtype=np.int64)
    assert data.X.format == "csc" and data.X.dtype == np.int64
    np.testing.assert_array_equal(data.X.toarray(), expected)
    assert list(data.obs_names) == [CELL_D, CELL_ZERO, CELL_C, CELL_B, CELL_A]
    assert CELL_COVID not in data.obs_names
    assert list(data.obs["donor_id"]) == ["Non-COV2", "IPF1", "Non-COV1", "IPF2", "IPF1"]
    assert list(data.obs["source_coarse_cell_type"]) == ["immune", "epithelial", "endothelial", "stroma", "epithelial"]
    assert data.obs.loc[CELL_D, "source_sample_id"] == "LUNG010"
    assert list(data.var["source_gene_name"]) == ["GENE_A", "MT-CO1", "GENE_A"]
    assert data.var_names.is_unique
    np.testing.assert_array_equal(data.obs["total_counts"], [4, 0, 5, 3, 6])
    np.testing.assert_array_equal(data.obs["n_genes_by_counts"], [2, 0, 1, 2, 2])
    assert bool(data.obs.loc[CELL_ZERO, "flag_zero_counts"])
    assert np.isnan(data.obs.loc[CELL_ZERO, "pct_counts_mt"])
    assert bool(data.obs.loc[CELL_C, "flag_high_mt_gt25pct"])
    summary = json.loads((output / "summary.json").read_text())
    assert summary["status"] == "technical_qc_passed"
    assert summary["full_count_file_read_to_eof"] and summary["gzip_trailer_crc_passed"]
    assert summary["source_global_metadata_cells"] == 6
    assert summary["selected_cells_retained"] == 5
    assert summary["gene_rows_retained"] == 3
    assert summary["duplicate_source_gene_rows_retained"] == 1
    assert summary["donors"]["IPF1"]["retained_source_cells"] == 2
    assert "confounded" in summary["procurement_confounding"]
    assert (output / "cell_qc.tsv").is_file() and (output / "donor_qc.tsv").is_file()

    failed_cases = []
    for failure, expected_error in [
        ("missing_metadata_cell", "global barcode sets differ"),
        ("fractional_count", "Non-integer"),
    ]:
        raw, manifest, output = create_case(args.work_dir / failure, failure)
        result = invoke(raw, manifest, output)
        assert result.returncode != 0, f"Invalid fixture accepted: {failure}"
        summary = json.loads((output / "summary.json").read_text())
        assert summary["status"] == "failed" and expected_error in summary["error"]
        assert not summary["full_count_file_read_to_eof"]
        assert not (output / "GSE159585_raw_qc.h5ad").exists()
        failed_cases.append(failure)
    print(json.dumps({
        "status": "passed", "valid_matrix_shape": [5, 3],
        "failure_cases_rejected": failed_cases, "work_dir": str(args.work_dir),
    }), flush=True)


if __name__ == "__main__":
    main()

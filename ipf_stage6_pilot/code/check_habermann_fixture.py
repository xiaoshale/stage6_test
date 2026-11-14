#!/usr/bin/env python3
"""Run a hand-specified, deliberately scrambled GSE135893 technical fixture.

The five-column source includes one unannotated same-library barcode that must
stay excluded, four manifest donors, duplicate gene symbols and a retained
zero-count source-annotated cell. This fixture is not biological data.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--qc-script", type=Path,
                        default=Path(__file__).with_name("qc_habermann.py"))
    args = parser.parse_args()
    work_dir = args.work_dir.resolve()
    if work_dir.exists() and any(work_dir.iterdir()):
        parser.error("--work-dir must be absent or empty; earlier results are preserved")
    raw_dir = work_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    texts = {
        "GSE135893_genes.tsv.gz": "MT-CO1\nGAPDH\nGAPDH\n",
        "GSE135893_barcodes.tsv.gz": (
            "F01157_CCC\nF00431_UNANNOTATED\nF01214_BBB\nF00431_AAA\nF01174_DDD\n"),
        # Columns are C, excluded U, B, A, D. D is explicitly all zero.
        "GSE135893_matrix.mtx.gz": (
            "%%MatrixMarket matrix coordinate integer general\n% hand-specified fixture\n"
            "3 5 7\n2 1 7\n1 2 99\n1 3 1\n2 3 2\n3 3 3\n1 4 2\n2 4 4\n"),
    }
    for filename, text in texts.items():
        with gzip.open(raw_dir / filename, "wt") as handle:
            handle.write(text)
    # Metadata order A,D,C,B deliberately disagrees with matrix order C,U,B,A,D.
    rows = [
        ["F00431_AAA", "F00431", 6, 2, "IPF", "TILD001", "NTI", "Disease", "Epithelial"],
        ["F01174_DDD", "F01174", 0, 0, "Control", "THD0002", "DNA", "Control", "Mesenchymal"],
        ["F01157_CCC", "F01157", 7, 1, "Control", "THD0001", "DNA", "Control", "Mesenchymal"],
        ["F01214_BBB", "F01214", 6, 3, "IPF", "TILD010", "NTI", "Disease", "Epithelial"],
    ]
    with gzip.open(raw_dir / "GSE135893_IPF_metadata.csv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "orig.ident", "nCount_RNA", "nFeature_RNA", "Diagnosis",
                         "Sample_Name", "Sample_Source", "Status", "population"])
        writer.writerows(rows)
    manifest = {"primary_data_group": [{
        "accession": "GSE135893",
        "full_matrix_header": {"genes": 3, "barcode_columns": 5, "nonzero_entries": 7},
        "selected_donors": [
            {"donor_id": "TILD001", "diagnosis": "IPF", "library_ids": ["F00431"],
             "provided_broad_lineage_cell_counts": {"Epithelial": 1}},
            {"donor_id": "TILD010", "diagnosis": "IPF", "library_ids": ["F01214"],
             "provided_broad_lineage_cell_counts": {"Epithelial": 1}},
            {"donor_id": "THD0001", "diagnosis": "Control", "library_ids": ["F01157"],
             "provided_broad_lineage_cell_counts": {"Mesenchymal": 1}},
            {"donor_id": "THD0002", "diagnosis": "Control", "library_ids": ["F01174"],
             "provided_broad_lineage_cell_counts": {"Mesenchymal": 1}},
        ],
    }]}
    manifest_path = work_dir / "fixture_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    output_dir = work_dir / "qc"
    subprocess.run([sys.executable, str(args.qc_script.resolve(strict=True)),
                    "--raw-dir", str(raw_dir), "--manifest", str(manifest_path),
                    "--output-dir", str(output_dir), "--threads", "4"], check=True)

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    result = ad.read_h5ad(output_dir / "counts.h5ad")
    summary = json.loads((output_dir / "summary.json").read_text())
    # Expected retained C,B,A,D rows are specified independently of QC code.
    expected_counts = np.array([[0, 7, 0], [1, 2, 3], [2, 4, 0], [0, 0, 0]])
    assert result.shape == (4, 3)
    assert sparse.issparse(result.X)
    np.testing.assert_array_equal(result.X.toarray(), expected_counts)
    np.testing.assert_array_equal(result.layers["counts"].toarray(), expected_counts)
    assert result.obs["source_barcode"].tolist() == ["F01157_CCC", "F01214_BBB", "F00431_AAA", "F01174_DDD"]
    assert result.obs["donor_id"].tolist() == ["THD0001", "TILD010", "TILD001", "THD0002"]
    assert result.obs["library_id"].tolist() == ["F01157", "F01214", "F00431", "F01174"]
    assert result.obs["diagnosis"].tolist() == ["Control", "IPF", "IPF", "Control"]
    assert result.obs["source_broad_annotation"].tolist() == ["Mesenchymal", "Epithelial", "Epithelial", "Mesenchymal"]
    assert result.var["source_gene_symbol"].tolist() == ["MT-CO1", "GAPDH", "GAPDH"]
    assert result.var_names.is_unique
    assert result.obs["total_counts"].tolist() == [7, 6, 6, 0]
    assert result.obs["n_genes_by_counts"].tolist() == [1, 3, 2, 0]
    assert result.obs["total_counts_mt"].tolist() == [0, 1, 2, 0]
    np.testing.assert_allclose(result.obs["pct_counts_mt"], [0, 100 / 6, 100 / 3, np.nan], equal_nan=True)
    assert result.obs["flag_low_genes"].tolist() == [True, True, True, True]
    assert result.obs["flag_high_mito"].tolist() == [False, False, True, False]
    assert result.obs["flag_zero_counts"].tolist() == [False, False, False, True]
    assert result.obs["source_nCount_RNA_matches_computed"].all()
    assert result.obs["source_nFeature_RNA_matches_computed"].all()
    assert summary["selected_total_counts"] == 19
    assert summary["selected_nnz"] == 6
    assert summary["global_counts_axes_audit"]["matrix_columns_without_this_author_annotation"] == 1
    assert summary["global_counts_axes_audit"]["unannotated_columns_included_in_pilot"] == 0
    assert summary["duplicate_source_gene_symbol_rows"] == 1
    assert summary["technical_flag_counts"]["flag_high_mito"] == 1
    assert summary["technical_flag_counts"]["flag_zero_counts"] == 1
    assert all(summary["reload_validation"].values())
    csv_result = pd.read_csv(output_dir / "QC.csv", index_col=0)
    assert csv_result.shape[0] == 4
    assert csv_result["source_barcode"].tolist() == result.obs["source_barcode"].tolist()
    assert pd.read_csv(output_dir / "source_broad_annotation_QC.csv")["cells"].sum() == 4
    assert not result.obsm and not result.obsp
    verdict = {"fixture": "PASS", "matrix_order_and_exact_annotation_join": True,
               "same_library_unannotated_barcode_excluded": True, "four_donor_library_mappings": True,
               "original_sparse_counts_roundtrip": True, "mitochondrial_and_zero_count_flags": True,
               "all_selected_cells_including_zero_count_cell_retained": True,
               "duplicate_gene_symbols_preserved": True}
    (work_dir / "fixture_validation.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()

"""Audit fixed marker detectability in existing raw-count checkpoints, read-only."""
import argparse
import csv
import json
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

PANEL = {
    "epithelial": "EPCAM KRT8 KRT18 KRT5 SCGB1A1 SFTPC".split(),
    "stromal": "DCN LUM COL1A1 COL1A2 PDGFRA".split(),
    "immune": "PTPRC LYZ CD3D".split(),
    "endothelial": "PECAM1 VWF CLDN5".split(),
}
SOURCES = [
    ("GSE135893", "counts.h5ad", "source_gene_symbol"),
    ("GSE159354", "*_raw_qc.h5ad", "gene_symbol"),
    ("GSE159585", "GSE159585_raw_qc.h5ad", "source_gene_name"),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    sources, output_rows = [], 0
    columns = ["resource", "source_file", "scope", "donor_id", "library_id",
               "source_coarse_field", "source_coarse_label", "marker_set", "marker",
               "present_in_source_var", "measurement_status", "source_symbol_rows", "duplicate_symbol_rows",
               "source_cells", "measured_cells", "nonzero_cells", "fraction_nonzero"]
    with (args.output_dir / "marker_capture.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for resource, pattern, symbol_field in SOURCES:
            files = sorted((args.run_dir / resource).glob(pattern))
            if not files:
                raise FileNotFoundError(f"No {resource}/{pattern} checkpoints")
            for path in files:
                data = ad.read_h5ad(path, backed="r")
                try:
                    if symbol_field not in data.var or "donor_id" not in data.obs:
                        raise ValueError(f"Missing source symbol/donor field in {path}")
                    library = "library_id" if "library_id" in data.obs else "gsm"
                    if library not in data.obs:
                        raise ValueError(f"Missing source library field in {path}")
                    coarse = next((k for k in ["source_broad_annotation", "source_coarse_cell_type"]
                                   if k in data.obs), None)
                    symbols = data.var[symbol_field].astype(str).to_numpy()
                    rows = {marker: np.flatnonzero(symbols == marker)
                            for markers in PANEL.values() for marker in markers}
                    present = np.array(sorted({int(i) for indices in rows.values() for i in indices}))
                    block = data.X[:, present] if len(present) else None
                    if block is not None and not sparse.issparse(block):
                        raise ValueError(f"Expected original sparse X in {path}")
                    groups = []
                    for scope, fields in [("donor_library", ["donor_id", library]),
                                          ("source_coarse_annotation", ["donor_id", library, coarse])]:
                        if None in fields:
                            continue
                        for key, indices in data.obs.groupby(fields, observed=True, sort=False,
                                                             dropna=False).indices.items():
                            groups.append((scope, key, indices))
                    for marker_set, markers in PANEL.items():
                        for marker in markers:
                            indices = rows[marker]
                            detected = None
                            if len(indices):
                                detected = np.asarray((block[:, np.searchsorted(present, indices)] > 0)
                                                      .sum(axis=1)).ravel() > 0
                            for scope, key, cells in groups:
                                measured = len(cells) if len(indices) else 0
                                nonzero = int(detected[cells].sum()) if measured else None
                                writer.writerow(dict(zip(columns, [
                                    resource, str(path), scope, str(key[0]), str(key[1]),
                                    coarse or "", str(key[2]) if len(key) > 2 else "",
                                    marker_set, marker, bool(len(indices)), "measured" if len(indices) else "unmeasured",
                                    len(indices), max(0, len(indices) - 1), len(cells), measured,
                                    nonzero, nonzero / measured if measured else None])))
                                output_rows += 1
                    sources.append({"resource": resource, "file": str(path), "cells": data.n_obs,
                                    "gene_rows": data.n_vars, "source_symbol_field": symbol_field,
                                    "library_field": library, "source_coarse_field": coarse})
                    print(json.dumps({"resource": resource, "file": path.name,
                                      "cells": data.n_obs, "status": "audited"}), flush=True)
                finally:
                    data.file.close()
    summary = {"status": "completed", "fixed_marker_sets": PANEL, "marker_count": 17, "source_checkpoints": sources,
               "output_rows": output_rows, "output": "marker_capture.tsv",
               "boundary": "Technical detectability only. Conventional marker-set labels are not cell assignments. "
                           "No filtering, normalization, disease comparison, state annotation or biological hypothesis. "
                           "All existing cells retained. Missing source symbols are unmeasured, not zero expression. "
                           "Duplicate symbol rows use detection union only; original X is never merged or rewritten. "
                           "Coarse labels, when available, are unchanged author labels; FACS fractions are not cell identities."}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()

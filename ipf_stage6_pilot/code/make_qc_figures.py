"""Plot descriptive raw-count distributions from the completed QC tables."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--run-dir", type=Path, required=True)
args = p.parse_args()
out = args.run_dir / "figures"
out.mkdir(exist_ok=False)
tables = {
    "GSE135893": (pd.read_csv(args.run_dir / "GSE135893/QC.csv"), "source-annotated cells", "n_genes_by_counts", "total_counts"),
    "GSE159354": (pd.concat([pd.read_csv(f) for f in sorted((args.run_dir / "GSE159354").glob("*_cell_qc.csv"))], ignore_index=True), "cells from sixteen sorted libraries", "n_genes_by_counts", "total_counts"),
    "GSE159585": (pd.read_csv(args.run_dir / "GSE159585/cell_qc.tsv", sep="\t"), "source nuclei (procurement-confounded cohort)", "n_genes_by_counts", "total_counts"),
    "Franzen": (pd.concat([pd.read_csv(f, sep="\t") for f in sorted((args.run_dir / "Franzen").glob("*/capture_position_qc.tsv"))], ignore_index=True), "author in-tissue capture positions from six sections", "n_features_by_counts", "n_counts"),
}
for resource, (frame, unit, features, counts) in tables.items():
    if resource == "Franzen":
        frame = frame.loc[frame["in_tissue"] == 1]
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    for axis, field, label in [(axes[0], counts, "Original total counts"), (axes[1], features, "Detected source features")]:
        values = pd.to_numeric(frame[field], errors="raise").to_numpy()
        axis.hist(np.log10(values + 1), bins=60, color="#365f80", edgecolor="white", linewidth=0.3)
        axis.set_xlabel(f"log10({label.lower()} + 1)")
        axis.set_ylabel("Measured units")
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_title(label, fontsize=10)
    figure.suptitle(f"{resource}: descriptive technical QC", fontsize=13)
    figure.text(0.5, 0.015, f"{len(frame):,} {unit}. Units are nested within donors/libraries; no disease-effect test.", ha="center", fontsize=7)
    figure.tight_layout(rect=(0, 0.05, 1, 0.94))
    figure.savefig(out / f"{resource}_technical_qc.png", dpi=160)
    plt.close(figure)
    print(resource, "QC figure saved", flush=True)

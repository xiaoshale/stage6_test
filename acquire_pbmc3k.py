#!/usr/bin/env python3
"""Acquire Scanpy's referenced raw PBMC3k H5AD as a technical control only."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import urllib.request


SOURCE_URL = "https://falexwolf.de/data/pbmc3k_raw.h5ad"
SOURCE_DEFINITION_URL = (
    "https://raw.githubusercontent.com/scverse/scanpy/1.11.5/"
    "src/scanpy/datasets/_datasets.py"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        parser.error("--output-dir must be absent or empty; previous downloads are preserved")
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "2"

    import anndata as ad
    h5ad_path = output_dir / "pbmc3k_raw.h5ad"
    digest = hashlib.sha256()
    downloaded_bytes = 0
    with urllib.request.urlopen(SOURCE_URL, timeout=90) as response, h5ad_path.open("xb") as target:
        expected_bytes = response.headers.get("Content-Length")
        for block in iter(lambda: response.read(1024 * 1024), b""):
            target.write(block)
            digest.update(block)
            downloaded_bytes += len(block)
    downloaded_at = datetime.now(timezone.utc).isoformat()
    if expected_bytes is not None and downloaded_bytes != int(expected_bytes):
        raise ValueError("Downloaded H5AD size does not match Content-Length")

    # Read in backed read-only mode; the downloaded H5AD bytes remain unchanged.
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        shape = list(adata.shape)
        if tuple(shape) != (2700, 32738):
            raise ValueError(f"Unexpected PBMC3k count shape: {shape}")
        if not adata.obs_names.is_unique:
            raise ValueError("PBMC3k cell barcodes must be unique")
    finally:
        adata.file.close()
    record = {
        "source_url": SOURCE_URL,
        "source_kind": "Scanpy 1.11.5 reference H5AD derived from 10x PBMC3k",
        "source_definition_url": SOURCE_DEFINITION_URL,
        "downloaded_at_utc": downloaded_at,
        "purpose": "Technical PBMC3k control only; no IPF conclusions or donor replication claim",
        "downloaded_bytes": downloaded_bytes,
        "h5ad_file": h5ad_path.name,
        "h5ad_sha256": digest.hexdigest(),
        "raw_shape_cells_genes": shape,
    }
    (output_dir / "acquisition.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()

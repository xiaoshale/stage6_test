#!/usr/bin/env bash
set -euo pipefail

test_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
run_id="${1:?Provide a unique run identifier}"
if [[ ! "$run_id" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "Run identifier must contain only letters, numbers, underscores or hyphens." >&2
    exit 1
fi
mkdir -p "$test_root/logs" "$test_root/tmp" "$test_root/cache" "$test_root/runs"
export TMPDIR="$test_root/tmp"
export MPLCONFIGDIR="$test_root/cache/matplotlib"
export NUMBA_CACHE_DIR="$test_root/cache/numba"
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export NUMBA_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
# Linux virtual-address-space ceiling: 16 GiB. This is not a production sizing rule.
ulimit -v 16777216

if [[ ! -f "$test_root/input/acquisition.json" ]]; then
    "$test_root/.venv/bin/python" "$test_root/acquire_pbmc3k.py" \
        --output-dir "$test_root/input"
fi
"$test_root/.venv/bin/python" - "$test_root/input" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

input_dir = Path(sys.argv[1])
record = json.loads((input_dir / "acquisition.json").read_text())
digest = hashlib.sha256()
with (input_dir / "pbmc3k_raw.h5ad").open("rb") as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
if digest.hexdigest() != record["h5ad_sha256"]:
    raise SystemExit("Input H5AD does not match the acquisition record; analysis stopped.")
print("Input acquisition SHA-256 verified.", flush=True)
PY

run_dir="$test_root/runs/$run_id"
if [[ -e "$run_dir" ]]; then
    echo "Refusing to replace an existing run: $run_dir" >&2
    exit 1
fi
"$test_root/.venv/bin/python" -u "$test_root/run_smoke.py" \
    --input "$test_root/input/pbmc3k_raw.h5ad" \
    --output-dir "$run_dir" --threads 2 --seed 0

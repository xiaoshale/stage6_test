#!/usr/bin/env bash
set -euo pipefail

test_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
base_python="/home/data/t120631/miniconda3/envs/litreview/bin/python"
mkdir -p "$test_root/logs" "$test_root/tmp" "$test_root/cache"
export TMPDIR="$test_root/tmp"
export PIP_CACHE_DIR="$test_root/cache/pip"
export MPLCONFIGDIR="$test_root/cache/matplotlib"
export NUMBA_CACHE_DIR="$test_root/cache/numba"
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2

if [[ -e "$test_root/.venv" ]]; then
    echo "Refusing to replace an existing test environment." >&2
    exit 1
fi
"$base_python" -m venv "$test_root/.venv"
"$test_root/.venv/bin/python" -m pip install \
    --index-url https://pypi.org/simple --only-binary=:all: \
    -r "$test_root/requirements.txt"
"$test_root/.venv/bin/python" -m pip check
"$test_root/.venv/bin/python" -m pip freeze > "$test_root/environment.freeze.txt"
"$test_root/.venv/bin/python" -c \
    'import sys, scanpy, anndata, scipy, igraph; print("Environment imports PASS:", sys.version.split()[0], "scanpy", scanpy.__version__, "anndata", anndata.__version__, "scipy", scipy.__version__, "igraph", igraph.__version__)'

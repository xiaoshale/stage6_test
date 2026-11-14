"""Create only the explicitly authorized pilot directory; reuse the tested runtime."""
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path("/home/data/t120631/project/ipf_stage6_pilot")
if root.exists():
    raise SystemExit("Pilot root already exists; inspect rather than recreate it.")
os.umask(0o027)
root.mkdir()
for name in ["code", "config", "data/raw", "data/acquisition_records", "runs", "logs", "tmp", "cache"]:
    (root / name).mkdir(parents=True, exist_ok=True)
versions = {name: importlib.metadata.version(name) for name in ["scanpy", "anndata", "scipy", "numpy", "pandas", "h5py", "Pillow", "matplotlib"]}
(root / "environment.json").write_text(json.dumps({"created_utc": datetime.now(timezone.utc).isoformat(), "interpreter": sys.executable, "mode": "Reuse the existing validated smoke-test virtual environment without modifying packages.", "versions": versions, "computational_threads": 4, "planned_process_virtual_memory_ceiling_GiB": 48}, indent=2) + "\n")
with (root / "environment.freeze.txt").open("x") as f:
    subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=f, check=True)
print(json.dumps({"root": str(root), "versions": versions}, indent=2))

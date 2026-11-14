"""Run one declared component with immutable code snapshots and resource controls."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
p.add_argument("--component", choices=["GSE135893", "GSE159354", "Franzen", "GSE159585"], required=True)
p.add_argument("--run-id", required=True)
args = p.parse_args()
if not args.run_id.replace("_", "").replace("-", "").isalnum():
    raise SystemExit("Unsafe run identifier")
root = args.root.resolve()
run = root / "runs" / args.run_id
output = run / args.component
if output.exists():
    raise SystemExit("Output already exists; retain it and use a new run ID")
run.mkdir(parents=True, exist_ok=True)
source = run / "source" / args.component
source.mkdir(parents=True, exist_ok=False)
scripts = {"GSE135893": ["qc_habermann.py"], "GSE159354": ["qc_depianto.py", "qc_spatial.py"], "Franzen": ["qc_spatial.py"], "GSE159585": ["qc_snrna.py"]}[args.component]
hashes = {}
for name in scripts:
    shutil.copyfile(root / "code" / name, source / name)
    hashes[name] = hashlib.sha256((source / name).read_bytes()).hexdigest()
os.environ.update({"TMPDIR": str(root / "tmp"), "MPLCONFIGDIR": str(root / "cache/matplotlib"), "MPLBACKEND": "Agg", "PYTHONDONTWRITEBYTECODE": "1"})
for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[key] = "4"
resource.setrlimit(resource.RLIMIT_AS, (48 * 1024**3, 48 * 1024**3))
required = [r for r in json.loads((root / "config/acquisition_manifest.json").read_text())["files"] if r["cohort"] == args.component]
for row in required:
    path = root / "data/raw" / args.component / row["filename"]
    rec = json.loads((root / "data/acquisition_records" / (args.component + "__" + row["filename"] + ".json")).read_text())
    if rec.get("status") != "complete" or path.stat().st_size != rec["bytes"]:
        raise SystemExit("Incomplete input: " + row["filename"])
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            checksum.update(block)
    if checksum.hexdigest() != rec["sha256"]:
        raise SystemExit("Input checksum differs from acquisition record: " + row["filename"])
command = [sys.executable, "-u", str(source / scripts[0]), "--raw-dir", str(root / "data/raw" / args.component), "--manifest", str(root / "config/pilot_manifest.json"), "--output-dir", str(output)]
start = time.monotonic()
record = {"started_utc": datetime.now(timezone.utc).isoformat(), "command": command, "code_sha256": hashes, "threads": 4, "virtual_memory_ceiling_GiB": 48}
result = subprocess.run(command, check=False)
record.update({"exit_code": result.returncode, "completed_utc": datetime.now(timezone.utc).isoformat(), "elapsed_seconds": time.monotonic()-start, "peak_child_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024})
(run / (args.component + "_execution.json")).write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2), flush=True)
raise SystemExit(result.returncode)

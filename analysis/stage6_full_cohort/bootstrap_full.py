"""Create an isolated full-cohort run, preserving pilot outputs and source records."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.umask(0o027)
root = Path("/home/data/t120631/project/ipf_stage6_pilot")
run = root / "runs/preparation_001"
run.mkdir(exist_ok=False)
for name in ("code", "config", "logs", "rna_raw", "rna_review", "spatial_review", "genes", "review_pack", "dependencies", "acquisition/data"):
    (run / name).mkdir(parents=True, exist_ok=True)
for name, target in (("raw", root / "data/raw"), ("acquisition_records", root / "data/acquisition_records")):
    (run / "acquisition/data" / name).symlink_to(target, target_is_directory=True)
(run / "RUN_STATUS.json").write_text(json.dumps({
    "started_utc": datetime.now(timezone.utc).isoformat(),
    "authorized_task": "Full native cohort preparation and human QC/identity/histology review package",
    "status": "preparing", "raw_counts_policy": "Original counts and pilot checkpoints preserved",
    "review_boundary": "Provisional QC/identity evidence, no expert anatomical signoff or final eligibility",
    "runtime": "/home/data/t120631/project/ipf_stage6_smoke_test_20261113/.venv/bin/python",
    "per_process_threads": 4, "planned_peak_parallel_compute_processes": 2,
    "external_candidate_outcomes": "reserved",
}, indent=2) + "\n")
print(run)

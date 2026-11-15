"""Create a completion manifest only after all declared technical components pass."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile
import h5py

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
p.add_argument("--run-id", required=True)
args = p.parse_args()
root = args.root.resolve()
run = root / "runs" / args.run_id
success = run / "SUCCESS.json"
if success.exists():
    raise SystemExit("Completion manifest already exists; do not replace it")
checks = []


def check(label, passed):
    checks.append({"check": label, "passed": bool(passed)})
    if not passed:
        raise ValueError(label)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


declared = json.loads((root / "config/acquisition_manifest.json").read_text())["files"]
acquired = json.loads((root / "data/acquisition.json").read_text())
check("All 24 complete source files", acquired["status"] == "complete" and len(acquired["files"]) == len(declared) == 24 and not acquired["errors"])
for item in acquired["files"]:
    path = Path(item["path"])
    check("Source integrity: " + item["filename"], item["status"] == "complete" and path.stat().st_size == item["bytes"] and sha256(path) == item["sha256"])
    initial = next((r for r in item["attempts"] if r.get("offset") == 0 and "http_status" in r), None)
    check("Complete initial HTTP response: " + item["filename"], initial is not None and initial["http_status"] == 200 and not initial["headers"].get("Content-Range"))
archive = root / "data/raw/Franzen/hs_visium_spaceranger_output.zip"
with zipfile.ZipFile(archive) as z:
    check("Complete spatial ZIP CRC", z.testzip() is None)
reports = {c: json.loads((run / c / "summary.json").read_text()) for c in ["GSE135893", "GSE159354", "GSE159585", "Franzen"]}
check("Habermann complete read and source join", reports["GSE135893"]["technical_execution_complete"] and reports["GSE135893"]["selected_shape_cells_genes"] == [8577, 33694] and reports["GSE135893"]["source_metric_disagreements"] == {"nCount_RNA": 0, "nFeature_RNA": 0})
check("DePianto sixteen library roundtrips", reports["GSE159354"]["status"] == "technical_qc_passed" and len(reports["GSE159354"]["libraries"]) == 16 and all(r["roundtrip_passed"] for r in reports["GSE159354"]["libraries"]))
check("snRNA full EOF and selected-cell mapping", reports["GSE159585"]["status"] == "technical_qc_passed" and reports["GSE159585"]["selected_cells_retained"] == 28970 and reports["GSE159585"]["full_count_file_read_to_eof"] and reports["GSE159585"]["gzip_trailer_crc_passed"])
fr = reports["Franzen"]
check("Six original spatial donor sections", fr["status"] == "passed_technical_checks" and fr["sections_completed"] == fr["biological_donors"] == 6 and all(s["capture_positions"] == 4992 and s["reload_verified"] for s in fr["sections"]))
for c in reports:
    execution = json.loads((run / (c + "_execution.json")).read_text())
    check("Successful component process: " + c, execution["exit_code"] == 0)
    for name, expected in execution["code_sha256"].items():
        check("Executed code snapshot: " + c + "/" + name, sha256(run / "source" / c / name) == expected)
checkpoints = sorted(run.glob("*/*.h5ad")) + sorted(run.glob("Franzen/*/*.h5ad"))
check("Twenty-four raw-count checkpoints", len(checkpoints) == 24)
for path in checkpoints:
    with h5py.File(path, "r") as f:
        check("Sparse checkpoint header: " + str(path.relative_to(run)), f["X"].attrs["encoding-type"] in ("csr_matrix", "csc_matrix") and all(int(x) > 0 for x in f["X"].attrs["shape"]))
figures = sorted((run / "figures").glob("*.png"))
overlays = sorted((run / "Franzen").glob("*/coordinate_overlay.png"))
check("Technical plots and coordinate overlays", len(figures) == 4 and len(overlays) == 6)
artifacts = []
for path in sorted(run.rglob("*")):
    if path.is_file():
        artifacts.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path), "large_checkpoint_retained_on_server": path.suffix == ".h5ad"})
result = {"status": "technical_pilot_complete", "completed_utc": datetime.now(timezone.utc).isoformat(), "client_date": "2026-11-15", "root": str(root), "run_id": args.run_id, "checks": checks, "source_files": len(declared), "source_bytes": sum(r["bytes"] for r in acquired["files"]), "raw_checkpoints": len(checkpoints), "artifacts": artifacts, "biological_analysis_ready": False, "limits": ["Pilot-subset technical QC only; full-cohort cell/state coverage and final inclusion criteria remain pending.", "No cell/gene filtering, normalization, disease effects, state discovery, niche tests or target ranking.", "snRNA disease/procurement confounding and spatial anatomical/spot-mixture limitations remain.", "Doublet/ambient-RNA strategy, gene harmonization and cross-source donor independence are not resolved by format validation."]}
success.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": result["status"], "source_files": result["source_files"], "source_bytes": result["source_bytes"], "checkpoints": len(checkpoints), "checks_passed": len(checks)}, indent=2), flush=True)

"""Bundle only declared small pilot review artifacts; exclude raw data/H5AD."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
p.add_argument("--run-id", required=True)
p.add_argument("--output", type=Path, required=True)
args = p.parse_args()
root = args.root.resolve()
run = root / "runs" / args.run_id
if not (run / "SUCCESS.json").exists():
    raise SystemExit("Require the verified completion manifest before export")
paths = [root / "data/acquisition.json", root / "environment.json", root / "environment.freeze.txt", run / "SUCCESS.json"]
for pattern in ["*_execution.json", "*/summary.json", "GSE135893/donor_library_QC.csv", "GSE135893/source_broad_annotation_QC.csv", "GSE159354/library_qc.csv", "GSE159585/donor_qc.tsv", "Franzen/section_donor_qc.tsv", "Franzen/*/section_summary.json", "Franzen/*/coordinate_overlay.png", "figures/*.png", "capture_markers/marker_capture.tsv"]:
    paths.extend(sorted(run.glob(pattern)))
paths = sorted(set(paths))
if any(not path.is_file() or path.suffix == ".h5ad" for path in paths):
    raise SystemExit("Invalid small-artifact selection")
with zipfile.ZipFile(args.output, "x", compression=zipfile.ZIP_DEFLATED) as bundle:
    for path in paths:
        bundle.write(path, str(path.relative_to(root)))
h = hashlib.sha256(args.output.read_bytes()).hexdigest()
record = {"archive": str(args.output), "bytes": args.output.stat().st_size, "sha256": h, "files": len(paths), "large_data_included": False}
args.output.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2), flush=True)

"""Verify a complete public-source download relayed through the user's Mac."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--root", type=Path, required=True)
p.add_argument("--download-record", type=Path, required=True)
args = p.parse_args()
root = args.root.resolve()
client = json.loads(args.download_record.read_text())
rows = json.loads((root / "config/acquisition_manifest.json").read_text())["files"]
matches = [r for r in rows if r["cohort"] == client["cohort"] and r["filename"] == client["filename"] and r["url"] == client["url"]]
if len(matches) != 1 or client.get("status") != "complete":
    raise SystemExit("Relay does not identify one complete declared source file")
row = matches[0]
folder = root / "data/raw" / row["cohort"]
received = folder / (row["filename"] + ".relay")
destination = folder / row["filename"]
if destination.exists():
    raise SystemExit("Final destination already exists; do not replace it")
h = hashlib.sha256()
with received.open("rb") as f:
    for block in iter(lambda: f.read(8 * 1024**2), b""):
        h.update(block)
if received.stat().st_size != client["bytes"] or h.hexdigest() != client["sha256"]:
    raise SystemExit("SSH-received file size/hash differs from the full source download")
if row.get("expected_bytes") and row["expected_bytes"] != received.stat().st_size:
    raise SystemExit("Relay size differs from readiness manifest")
record_path = root / "data/acquisition_records" / (row["cohort"] + "__" + row["filename"] + ".json")
previous = json.loads(record_path.read_text()) if record_path.exists() else None
received.rename(destination)
record = {**client, "path": str(destination), "transport": "Complete source download to temporary Mac storage, then SSH copy with equal server SHA256.", "relay_source_record": str(args.download_record), "server_verified_utc": datetime.now(timezone.utc).isoformat(), "prior_direct_server_record": previous, "prior_direct_partial_preserved": str(folder / (row["filename"] + ".part")), "server_verified_sha256": h.hexdigest()}
record_path.write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps({"status": "complete", "bytes": record["bytes"], "sha256": record["sha256"], "path": str(destination)}, indent=2))

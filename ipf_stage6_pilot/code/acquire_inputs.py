"""Resumable bounded-concurrency acquisition; preserve partial data on failures."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import threading
import time
import urllib.error
import urllib.request


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    rows = json.loads(args.manifest.read_text())["files"]
    record_dir = root / "data/acquisition_records"
    record_dir.mkdir(parents=True, exist_ok=True)
    lock = threading.Lock()

    def log(message):
        with lock:
            print(now(), message, flush=True)

    def acquire(row):
        folder = root / "data/raw" / row["cohort"]
        folder.mkdir(parents=True, exist_ok=True)
        dest = folder / row["filename"]
        partial = dest.with_name(dest.name + ".part")
        record_path = record_dir / (row["cohort"] + "__" + row["filename"] + ".json")
        record = {**row, "started_utc": now(), "status": "in_progress", "attempts": []}
        if dest.exists():
            if not record_path.exists():
                raise RuntimeError("Existing file has no acquisition record: " + str(dest))
            old = json.loads(record_path.read_text())
            if old.get("status") != "complete" or old.get("sha256") != digest(dest):
                raise RuntimeError("Existing file integrity/status mismatch: " + str(dest))
            log("VERIFIED existing " + row["filename"])
            return old
        if partial.exists() and record_path.exists():
            record = json.loads(record_path.read_text())
            if not record.get("entity_validator"):
                for earlier in record.get("attempts", []):
                    saved = earlier.get("headers", {})
                    etag = saved.get("ETag") or ""
                    validator = etag if etag and not etag.startswith("W/") else saved.get("Last-Modified")
                    if validator:
                        record["entity_validator"] = validator
                        break
        expected_total = row.get("expected_bytes") or record.get("expected_bytes_observed")
        for attempt in range(1, 4):
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"User-Agent": "IPF-Stage6-technical-pilot/1.0", "Accept-Encoding": "identity"}
            if offset:
                headers["Range"] = f"bytes={offset}-"
                if record.get("entity_validator"):
                    headers["If-Range"] = record["entity_validator"]
            event = {"attempt": attempt, "offset": offset, "started_utc": now()}
            record["attempts"].append(event)
            try:
                log(f"GET {row['cohort']}/{row['filename']} offset={offset}")
                request = urllib.request.Request(row["url"], headers=headers)
                with urllib.request.urlopen(request, timeout=120) as response:
                    event.update({"http_status": response.status, "actual_url": response.url, "headers": {k: response.headers.get(k) for k in ("Content-Length", "Content-Range", "ETag", "Last-Modified", "Content-Type")}})
                    if "text/html" in response.headers.get("Content-Type", "").lower():
                        raise RuntimeError("HTML response is not a source data file")
                    if offset:
                        m = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", ""))
                        if response.status != 206 or not m or int(m[1]) != offset:
                            raise RuntimeError("Server did not honor exact resume offset")
                        total = int(m[3])
                    else:
                        if response.status != 200 or response.headers.get("Content-Range"):
                            raise RuntimeError("Initial transfer must be a complete HTTP200 response")
                        total = int(response.headers["Content-Length"]) if response.headers.get("Content-Length") else None
                    etag = response.headers.get("ETag", "")
                    validator = etag if etag and not etag.startswith("W/") else response.headers.get("Last-Modified")
                    if record.get("entity_validator") and validator != record["entity_validator"]:
                        raise RuntimeError("Source entity changed during transfer; partial preserved")
                    if validator:
                        record["entity_validator"] = validator
                    if expected_total and total and expected_total != total:
                        raise RuntimeError("Source size changed from declared metadata")
                    expected_total = expected_total or total
                    record["expected_bytes_observed"] = expected_total
                    record_path.write_text(json.dumps(record, indent=2) + "\n")
                    last_report = time.monotonic()
                    with partial.open("ab") as out:
                        while True:
                            block = response.read(4 * 1024 * 1024)
                            if not block:
                                break
                            out.write(block)
                            if time.monotonic() - last_report > 30:
                                out.flush()
                                log(f"PROGRESS {row['filename']} bytes={out.tell()} expected={expected_total}")
                                last_report = time.monotonic()
                if expected_total is not None and partial.stat().st_size != expected_total:
                    raise RuntimeError("Incomplete response length")
                sha = digest(partial)
                if row.get("expected_sha256") and sha != row["expected_sha256"]:
                    raise RuntimeError("Metadata differs from readiness snapshot; keep partial for review")
                partial.rename(dest)
                record.update({"status": "complete", "completed_utc": now(), "bytes": dest.stat().st_size, "sha256": sha, "path": str(dest)})
                record_path.write_text(json.dumps(record, indent=2) + "\n")
                log(f"COMPLETE {row['filename']} bytes={record['bytes']} sha256={sha}")
                return record
            except Exception as exc:
                event["error"] = repr(exc)
                record.update({"status": "partial_or_failed", "last_error": repr(exc)})
                record_path.write_text(json.dumps(record, indent=2) + "\n")
                log(f"RETRYABLE {row['filename']}: {exc!r}")
                if attempt < 3:
                    time.sleep(5)
        raise RuntimeError("Acquisition failed after three attempts: " + row["filename"])

    results, errors = [], []
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = {pool.submit(acquire, row): row for row in rows}
        for future in as_completed(jobs):
            try:
                results.append(future.result())
            except Exception as exc:
                errors.append({"file": jobs[future]["filename"], "error": repr(exc)})
    manifest = {"completed_utc": now(), "status": "complete" if not errors else "failed", "files": results, "errors": errors}
    (root / "data/acquisition.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if errors:
        raise SystemExit(1)
    log("ALL ACQUISITIONS COMPLETE")


if __name__ == "__main__":
    main()

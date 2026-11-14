#!/usr/bin/env python3
"""Exercise acquisition against four tiny localhost HTTP cases using only stdlib.

The work directory is retained. No public downloads or existing project data are
used. Production main() runs unchanged; only its retry sleep is replaced locally.
"""

import argparse
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import json
from pathlib import Path
import sys
import threading
import time
from types import SimpleNamespace


GOLDEN = b"ACGT\nIPF_pilot_counts\nend\n"
CHANGED = b"TGCA\nIPF_pilot_counts\nend\n"
PREFIX = GOLDEN[:7]
OLD_ETAG = '"fixture-v1"'
NEW_ETAG = '"fixture-v2"'


class FixtureHandler(BaseHTTPRequestHandler):
    def log_message(self, *_args):
        pass

    def do_GET(self):
        range_header = self.headers.get("Range")
        if_range = self.headers.get("If-Range")
        status, body, etag, content_range = 200, GOLDEN, OLD_ETAG, None
        if self.path == "/unexpected206":
            status, body = 206, PREFIX
            content_range = f"bytes 0-{len(PREFIX) - 1}/{len(GOLDEN)}"
        elif self.path == "/resume":
            if range_header == "bytes=7-" and if_range == OLD_ETAG:
                status, body = 206, GOLDEN[7:]
                content_range = f"bytes 7-{len(GOLDEN) - 1}/{len(GOLDEN)}"
        elif self.path == "/changed":
            # A same-size new entity ignores the stale If-Range and sends200.
            body, etag = CHANGED, NEW_ETAG
        elif self.path != "/complete":
            self.send_error(404)
            return
        self.server.requests.append({"path": self.path, "range": range_header,
                                     "if_range": if_range, "status": status})
        self.send_response(status)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("ETag", etag)
        if content_range:
            self.send_header("Content-Range", content_range)
        self.end_headers()
        self.wfile.write(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--acquisition-script", type=Path,
                        default=Path(__file__).with_name("acquire_inputs.py"))
    args = parser.parse_args()
    work_dir = args.work_dir.resolve()
    if work_dir.exists() and any(work_dir.iterdir()):
        parser.error("--work-dir must be absent or empty; earlier fixture files are retained")
    work_dir.mkdir(parents=True, exist_ok=True)
    source_path = args.acquisition_script.resolve(strict=True)
    spec = importlib.util.spec_from_file_location("fixture_acquisition", source_path)
    module = importlib.util.module_from_spec(spec)
    original_bytecode_setting = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = original_bytecode_setting
    # Keep production logic untouched and do not replace stdlib time globally.
    module.time = SimpleNamespace(monotonic=time.monotonic, sleep=lambda _seconds: None)
    assert len(CHANGED) == len(GOLDEN) and CHANGED != GOLDEN
    server = ThreadingHTTPServer(("127.0.0.1", 0), FixtureHandler)
    server.requests = []
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    results = []
    try:
        for case, endpoint, resume, succeeds in [
            ("complete200", "/complete", False, True),
            ("resume_same_entity", "/resume", True, True),
            ("initial206_rejected", "/unexpected206", False, False),
            ("changed_entity200_preserves_partial", "/changed", True, False),
        ]:
            root = work_dir / case
            root.mkdir()
            row = {"cohort": "fixture", "filename": "golden.bin",
                   "url": f"http://127.0.0.1:{server.server_port}{endpoint}"}
            # Deliberately leave unexpected206 without a size/hash: HTTP status
            # validation, rather than a coincidental length mismatch, must reject.
            if succeeds:
                row.update(expected_bytes=len(GOLDEN),
                           expected_sha256=hashlib.sha256(GOLDEN).hexdigest())
            destination = root / "data/raw/fixture/golden.bin"
            partial = destination.with_name("golden.bin.part")
            record_path = root / "data/acquisition_records/fixture__golden.bin.json"
            if resume:
                partial.parent.mkdir(parents=True)
                partial.write_bytes(PREFIX)
                record_path.parent.mkdir(parents=True)
                record_path.write_text(json.dumps({**row, "status": "partial_or_failed",
                    "entity_validator": OLD_ETAG, "expected_bytes_observed": len(GOLDEN),
                    "attempts": []}))
            manifest_path = root / "manifest.json"
            manifest_path.write_text(json.dumps({"files": [row]}))
            before_requests = len(server.requests)
            original_argv = sys.argv
            exit_code = 0
            try:
                sys.argv = [str(source_path), "--root", str(root), "--manifest", str(manifest_path)]
                module.main()
            except SystemExit as exc:
                exit_code = exc.code
            finally:
                sys.argv = original_argv
            record = json.loads(record_path.read_text())
            acquisition = json.loads((root / "data/acquisition.json").read_text())
            requests = server.requests[before_requests:]
            if succeeds:
                assert exit_code == 0 and acquisition["status"] == "complete", case
                assert destination.read_bytes() == GOLDEN, case
                assert not partial.exists(), case
                assert record["status"] == "complete" and record["bytes"] == len(GOLDEN), case
                assert record["sha256"] == hashlib.sha256(GOLDEN).hexdigest(), case
                assert record["entity_validator"] == OLD_ETAG, case
                assert len(requests) == 1, case
                assert requests[0]["status"] == (206 if resume else 200), case
            else:
                assert exit_code == 1 and acquisition["status"] == "failed", case
                assert record["status"] == "partial_or_failed" and not destination.exists(), case
                assert len(requests) == 3, case
                if resume:
                    assert partial.read_bytes() == PREFIX, case
                    assert record["entity_validator"] == OLD_ETAG, case
                    assert all(r["status"] == 200 for r in requests), case
                else:
                    assert not partial.exists(), case
                    assert all(r["status"] == 206 for r in requests), case
                    assert all("complete HTTP200" in a["error"] for a in record["attempts"]), case
            if resume:
                assert all(r["range"] == "bytes=7-" and r["if_range"] == OLD_ETAG for r in requests), case
            else:
                assert all(r["range"] is None and r["if_range"] is None for r in requests), case
            results.append({"case": case, "result": "PASS", "exit_code": exit_code,
                            "requests": requests, "original_partial_preserved": resume and not succeeds})
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)
    verdict = {"fixture": "PASS", "source_script": str(source_path),
               "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
               "golden_hex": GOLDEN.hex(), "golden_sha256": hashlib.sha256(GOLDEN).hexdigest(),
               "cases": results, "retry_sleep_skipped_only_in_fixture_module": True,
               "work_directory_retained": str(work_dir)}
    (work_dir / "fixture_validation.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()

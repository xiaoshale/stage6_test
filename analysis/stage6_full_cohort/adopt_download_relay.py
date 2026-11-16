"""Adopt SHA-verified complete public downloads; preserve partial/older records."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024**2), b''): h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--relay', type=Path, required=True)
    a = p.parse_args()
    manifest = json.loads((a.run / 'config/acquisition_manifest.json').read_text())['files']
    declared = {(r['cohort'], r['filename'], r['url']): r for r in manifest}
    for record_path in sorted((a.relay / 'data/acquisition_records').glob('*.json')):
        client = json.loads(record_path.read_text())
        key = (client['cohort'], client['filename'], client['url'])
        if key not in declared or client.get('status') != 'complete': raise ValueError('Undeclared/incomplete relay')
        received = a.relay / 'data/raw' / client['cohort'] / client['filename']
        if sha(received) != client['sha256'] or received.stat().st_size != client['bytes']: raise ValueError('Relay transport mismatch')
        dest = a.root / 'data/raw' / client['cohort'] / client['filename']
        server_record = a.root / 'data/acquisition_records' / record_path.name
        previous = json.loads(server_record.read_text()) if server_record.exists() else None
        if dest.exists():
            if sha(dest) != client['sha256'] or not previous or previous.get('status') != 'complete': raise ValueError('Existing destination mismatch')
            print('Verified matching existing', dest.name, flush=True)
            continue
        received.rename(dest)
        server_record.write_text(json.dumps({**client, 'path': str(dest), 'transport': 'Complete public-source Mac download relayed over SSH; equal server size and SHA256', 'server_verified_utc': datetime.now(timezone.utc).isoformat(), 'prior_server_record': previous, 'partial_preserved': str(dest.with_name(dest.name + '.part'))}, indent=2) + '\n')
        print('Adopted verified relay', dest.name, flush=True)


if __name__ == '__main__': main()

#!/usr/bin/env python3
"""Acquire bounded official identity assets, never expression or sequence data."""
import argparse
import csv
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
import urllib.parse
import urllib.request

MAX_TOTAL = 30 * 1024 * 1024
QUERY = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
         '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">'
         '<Dataset name="hsapiens_gene_ensembl" interface="default">'
         '<Attribute name="ensembl_gene_id"/><Attribute name="ensembl_gene_id_version"/>'
         '<Attribute name="external_gene_name"/><Attribute name="chromosome_name"/>'
         '<Attribute name="gene_biotype"/></Dataset></Query>')
HGNC_BASE = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(url, path, remaining):
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite reference asset: {path}')
    request = urllib.request.Request(url, headers={'User-Agent': 'Stage6-reference-preparation/1.0'})
    with urllib.request.urlopen(request, timeout=45) as response:
        data = response.read(remaining + 1)
        if len(data) > remaining:
            raise ValueError(f'Official download exceeds remaining {remaining}-byte budget')
        if response.status != 200 or not data:
            raise ValueError('Empty or unsuccessful reference download')
        headers = {key: response.headers.get(key) for key in
                   ('Last-Modified', 'ETag', 'Content-Type', 'Content-Length')}
        final_url = response.geturl()
    path.write_bytes(data)
    return {'path': path.name, 'url': url, 'final_url': final_url,
            'bytes': len(data), 'sha256': sha256(path), 'http_headers': headers,
            'acquired_utc': datetime.now(timezone.utc).isoformat()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--ensembl-host', default='https://jun2026.archive.ensembl.org')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / 'reference_manifest.json'
    if manifest_path.exists():
        raise FileExistsError(f'Refusing to overwrite {manifest_path}')
    assets = {}
    total = 0
    service = args.ensembl_host.rstrip('/') + '/biomart/martservice'
    for key, name, url in (
        ('ensembl_registry', 'ensembl116_registry.xml', service + '?type=registry'),
        ('ensembl_genes', 'ensembl116_genes.tsv', service + '?query=' + urllib.parse.quote(QUERY)),
        ('hgnc_complete', 'hgnc_complete_set.tsv', HGNC_BASE + 'hgnc_complete_set.txt'),
        ('hgnc_withdrawn', 'hgnc_withdrawn.tsv', HGNC_BASE + 'withdrawn.txt'),
    ):
        asset = fetch(url, args.output_dir / name, MAX_TOTAL - total)
        assets[key] = asset
        total += asset['bytes']
        print(json.dumps({'asset': key, 'bytes': asset['bytes'], 'sha256': asset['sha256']}), flush=True)
    registry = (args.output_dir / assets['ensembl_registry']['path']).read_text()
    if 'ensembl_mart_116' not in registry:
        raise ValueError('BioMart registry does not establish Ensembl release 116')
    with (args.output_dir / assets['ensembl_genes']['path']).open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        expected = ['Gene stable ID', 'Gene stable ID version', 'Gene name',
                    'Chromosome/scaffold name', 'Gene type']
        if header != expected:
            raise ValueError(f'Unexpected Ensembl header: {header}')
        rows = list(reader)
    if len(rows) < 30000 or any(len(row) != 5 or not row[0].startswith('ENSG') for row in rows):
        raise ValueError('Ensembl gene axis is incomplete or malformed')
    manifest = {
        'schema_version': 1, 'common_identity_target': 'Homo sapiens GRCh38 Ensembl 116 (June 2026)',
        'created_utc': datetime.now(timezone.utc).isoformat(), 'download_budget_bytes': MAX_TOTAL,
        'acquired_bytes': total, 'assets': assets, 'ensembl_gene_rows': len(rows),
        'source_counting_reference_status': 'Not recovered by this identity-asset acquisition',
        'history_status': 'HGNC withdrawal records, including replacement reports where supplied, acquired; complete HGNC/Ensembl split/merge history not acquired. Missing target-release source IDs remain held.',
        'boundary': 'Frozen downloaded identity assets; no counts, matrices, historical read reassignment or expression harmonization.',
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(str(manifest_path), flush=True)


if __name__ == '__main__':
    main()

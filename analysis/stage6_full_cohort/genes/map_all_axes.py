"""Map every native intake axis once by content, preserving each unit crosswalk."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from map_features import REQUIRED, read_tsv, sha256, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve(strict=True)
    if args.output_dir.exists():
        raise FileExistsError('Mapping batch output must be a new directory; previous results are preserved.')
    intake_path = run / 'config/cohort_intake.json'
    intake = json.loads(intake_path.read_text())
    expected = {unit['unit_uid']: (resource, unit)
                for resource, cohort in intake['resources'].items() for unit in cohort['units']}
    sources = []
    manifests = [intake_path]
    for resource in ('GSE135893', 'GSE159354', 'GSE159585'):
        manifest = run / 'rna_raw' / resource / 'inputs.json'
        manifests.append(manifest)
        records = json.loads(manifest.read_text())['records']
        for record in records:
            if record['resource'] != resource:
                raise ValueError('RNA input manifest resource mismatch')
            sources.append({'resource': resource, 'unit_id': record['library_id'],
                            'donor_id': record['donor_id'], 'axis_path': record['feature_axis'],
                            'declared_sha256': record['feature_axis_sha256'],
                            'expected_features': record['features']})
    spatial_manifest = run / 'spatial_review/summary.json'
    manifests.append(spatial_manifest)
    for record in json.loads(spatial_manifest.read_text())['sections']:
        axis = record['feature_axis']
        sources.append({'resource': 'Franzen', 'unit_id': record['section_id'],
                        'donor_id': record['donor_id'], 'axis_path': axis['path'],
                        'declared_sha256': axis['sha256'], 'expected_features': record['features']})
    identifiers = [f"{row['resource']}::{row['unit_id']}" for row in sources]
    if len(identifiers) != len(set(identifiers)) or set(identifiers) != set(expected):
        raise ValueError('Native axis inputs do not exactly cover the locked intake units')

    groups = {}
    axis_cache = {}
    unit_rows = []
    for record in sources:
        path = Path(record['axis_path']).resolve(strict=True)
        path.relative_to(run)
        if path not in axis_cache:
            raw_hash = sha256(path)
            rows, fields = read_tsv(path)
            if not set(REQUIRED).issubset(fields) or len(set(fields)) != len(fields):
                raise ValueError(f'Unexpected native feature schema: {path}')
            rows = [{key: row[key] for key in REQUIRED} for row in rows]
            normalized = json.dumps(rows, ensure_ascii=False, separators=(',', ':'), sort_keys=True).encode()
            content_hash = hashlib.sha256(normalized).hexdigest()
            axis_cache[path] = {'raw_sha256': raw_hash, 'content_sha256': content_hash,
                                'features': len(rows)}
            if content_hash not in groups:
                groups[content_hash] = {'representative': str(path), 'resources': set(),
                                        'units': [], 'features': len(rows)}
        axis_info = axis_cache[path]
        if axis_info['raw_sha256'] != record['declared_sha256'] or axis_info['features'] != record['expected_features']:
            raise ValueError(f'Native feature hash/length mismatch: {path}')
        unit_uid = f"{record['resource']}::{record['unit_id']}"
        if expected[unit_uid][1]['donor_id'] != record['donor_id']:
            raise ValueError(f'Feature unit/donor crosswalk mismatch: {unit_uid}')
        content_hash = axis_info['content_sha256']
        groups[content_hash]['resources'].add(record['resource'])
        groups[content_hash]['units'].append(unit_uid)
        unit_rows.append({'unit_uid': unit_uid, 'resource': record['resource'],
                          'donor_id': record['donor_id'], 'library_or_section_id': record['unit_id'],
                          'source_axis_path': str(path), 'source_axis_sha256': axis_info['raw_sha256'],
                          'source_features': axis_info['features'], 'content_sha256': content_hash})

    args.output_dir.mkdir(parents=True)
    summaries = {}
    mappings = {}
    for content_hash, group in sorted(groups.items()):
        label = '+'.join(sorted(group['resources'])) + '__' + content_hash[:12]
        folder = args.output_dir / 'axes' / label
        subprocess.run([sys.executable, str(Path(__file__).with_name('map_features.py')),
                        '--features', group['representative'], '--resource', '+'.join(sorted(group['resources'])),
                        '--references', str(args.references), '--output-dir', str(folder)], check=True)
        summary = json.loads((folder / 'mapping_summary.json').read_text())
        if summary['source_rows'] != group['features'] or summary['counts_read_or_changed'] is not False:
            raise ValueError('Mapping changed row accounting or count policy')
        mapping = folder / 'gene_mapping.tsv'
        mappings[content_hash] = {'label': label, 'path': str(mapping), 'sha256': sha256(mapping)}
        summaries[label] = {**summary, 'content_sha256': content_hash,
                            'native_units': len(group['units']), 'source_resources': sorted(group['resources'])}
    for row in unit_rows:
        row['mapping_path'] = mappings[row['content_sha256']]['path']
        row['mapping_sha256'] = mappings[row['content_sha256']]['sha256']
        row['mapping_status'] = 'metadata_mapped_human_review_pending'
    write_tsv(args.output_dir / 'unit_feature_mapping_crosswalk.tsv', unit_rows, list(unit_rows[0]))
    command = [sys.executable, str(Path(__file__).with_name('coverage.py'))]
    for info in mappings.values():
        command.extend(['--mapping', info['label'] + '=' + info['path']])
    command.extend(['--output-dir', str(args.output_dir / 'coverage')])
    subprocess.run(command, check=True)
    coverage = json.loads((args.output_dir / 'coverage/coverage_summary.json').read_text())
    result = {
        'status': 'complete_identity_metadata_human_review_pending',
        'native_units_mapped': len(unit_rows), 'units_by_resource': dict(Counter(r['resource'] for r in unit_rows)),
        'unique_native_axis_files': len(axis_cache), 'unique_axes_by_content': len(groups),
        'content_fingerprint': 'SHA256 of ordered five-field source metadata rows in canonical JSON; no gene renaming or row sorting',
        'unit_crosswalk_sha256': sha256(args.output_dir / 'unit_feature_mapping_crosswalk.tsv'),
        'input_manifest_sha256': {str(p): sha256(p) for p in manifests},
        'reference_manifest_sha256': sha256(args.references),
        'mapping_summaries': summaries, 'coverage_summary': coverage,
        'counts_read_or_changed': False,
        'boundary': 'Complete native axis identity/coverage metadata. Symbol-only matches remain provisional, held genes remain visible, no matrix projection or final biological eligibility claimed.',
    }
    (args.output_dir / 'batch_mapping_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('status', 'native_units_mapped', 'units_by_resource',
                                                'unique_native_axis_files', 'unique_axes_by_content')}, indent=2))


if __name__ == '__main__':
    main()

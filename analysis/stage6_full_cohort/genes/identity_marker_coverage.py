"""Audit fixed RNA identity seeds against native axes; never infer expression."""
import argparse
import ast
import json
from pathlib import Path

from map_features import read_tsv, sha256, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--seed-code', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    programs = None
    for node in ast.parse(args.seed_code.read_text()).body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == 'PROGRAMS' for target in node.targets):
            programs = ast.literal_eval(node.value)
    if not programs:
        raise ValueError('Fixed PROGRAMS dictionary not found in the actual RNA review code')
    markers = sorted({gene for _, genes, _ in programs.values() for gene in genes.split()})
    batch = json.loads(args.batch.read_text())
    crosswalk, _ = read_tsv(args.batch.parent / 'unit_feature_mapping_crosswalk.tsv')
    resource_axes = {}
    for row in crosswalk:
        resource_axes[(row['resource'], row['mapping_path'])] = row['content_sha256']
    rows = []
    totals = {}
    for (resource, mapping_path), axis_hash in sorted(resource_axes.items()):
        mappings, _ = read_tsv(mapping_path)
        marker_records = {}
        for marker in markers:
            found = [row for row in mappings if row['source_symbol'] == marker]
            status = ('not_in_source_export_coverage_unknown' if not found else
                      'multiple_source_rows_review' if len(found) > 1 else
                      'source_id_unique' if found[0]['primary_eligible'] == 'true' else
                      'symbol_only_provisional' if found[0]['provisional_unique'] == 'true' else
                      'source_label_present_mapping_held')
            marker_records[marker] = status
            rows.append({'resource': resource, 'axis_content_sha256': axis_hash, 'marker': marker,
                         'identity_programs': '|'.join(name for name, (_, genes, _) in programs.items() if marker in genes.split()),
                         'source_symbol_rows': len(found), 'identity_coverage_status': status,
                         'source_gene_ids': '|'.join(sorted({r['source_gene_id'] for r in found if r['source_gene_id']})),
                         'canonical_gene_ids': '|'.join(sorted({r['canonical_gene_id'] for r in found if r['canonical_gene_id']})),
                         'mapping_classes': '|'.join(sorted({r['mapping_class'] for r in found})),
                         'hold_reasons': '|'.join(sorted({r['hold_reason'] for r in found if r['hold_reason']}))})
        totals[resource + '__' + axis_hash[:12]] = {
            'markers': len(markers), 'source_label_present': sum(v != 'not_in_source_export_coverage_unknown' for v in marker_records.values()),
            'statuses': {state: sum(v == state for v in marker_records.values()) for state in sorted(set(marker_records.values()))},
            'mapping_path': mapping_path, 'mapping_sha256': sha256(mapping_path),
        }
    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(args.output_dir / 'identity_marker_coverage.tsv', rows, list(rows[0]))
    summary = {'markers': markers, 'fixed_marker_dictionary_sha256': sha256(args.seed_code),
               'batch_mapping_summary_sha256': sha256(args.batch), 'resources': totals,
               'boundary': 'Native feature-label and identity-mapping availability only; no counts read, detection, marker specificity, cell identity or expert review inferred. Missing exported labels remain coverage unknown.'}
    (args.output_dir / 'identity_marker_coverage_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'fixed_markers': len(markers), 'resource_axes': len(totals), 'resources': totals}, indent=2))


if __name__ == '__main__':
    main()

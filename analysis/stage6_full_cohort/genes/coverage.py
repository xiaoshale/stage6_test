#!/usr/bin/env python3
"""Describe feature coverage; absence is not converted to an expression zero."""
import argparse
import json
from pathlib import Path
from map_features import read_tsv, sha256, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mapping', action='append', required=True, metavar='LABEL=PATH')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    axes = {}
    for value in args.mapping:
        label, path = value.split('=', 1)
        if not label or label in axes:
            raise ValueError('Mapping labels must be nonempty and unique')
        rows, _ = read_tsv(path)
        axes[label] = {
            'strict': {row['canonical_gene_id'] for row in rows if row['primary_eligible'] == 'true'},
            'provisional': {row['canonical_gene_id'] for row in rows if row['provisional_unique'] == 'true'},
            'candidates': {gene for row in rows for gene in row['candidate_gene_ids'].split('|') if gene},
            'path': path, 'sha256': sha256(path),
        }
    universe = sorted(set().union(*(axis['candidates'] for axis in axes.values())))
    rows = []
    for gene in universe:
        row = {'canonical_gene_id': gene}
        for label, axis in axes.items():
            row[label] = ('source_id_unique' if gene in axis['strict'] else
                          'symbol_only_provisional' if gene in axis['provisional'] else
                          'mapping_held' if gene in axis['candidates'] else
                          'not_in_mapped_export_coverage_unknown')
        rows.append(row)
    strict = set.intersection(*(axis['strict'] for axis in axes.values()))
    candidate = set.intersection(*(axis['strict'] | axis['provisional'] for axis in axes.values()))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in ('gene_coverage.tsv', 'shared_strict_features.tsv', 'shared_unique_candidates.tsv', 'coverage_summary.json'):
        if (args.output_dir / name).exists():
            raise FileExistsError(f'Refusing to overwrite {name}')
    write_tsv(args.output_dir / 'gene_coverage.tsv', rows, ['canonical_gene_id', *axes])
    for name, genes in (('shared_strict_features.tsv', strict), ('shared_unique_candidates.tsv', candidate)):
        write_tsv(args.output_dir / name, [{'canonical_gene_id': gene} for gene in sorted(genes)], ['canonical_gene_id'])
    summary = {
        'inputs': {label: {'path': axis['path'], 'sha256': axis['sha256'], 'strict_features': len(axis['strict']),
                           'symbol_only_provisional_features': len(axis['provisional'])} for label, axis in axes.items()},
        'union_candidate_genes': len(universe), 'shared_strict_features': len(strict),
        'shared_unique_candidates_including_provisional': len(candidate),
        'boundary': 'Metadata coverage only. Export absence is not off-panel status or zero expression. No counts, gene projection, imputation, expression harmonization or final cohort lock.',
    }
    (args.output_dir / 'coverage_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Conservatively map exported native features without reading or changing counts."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

RULE_VERSION = 'stage6-gene-identity-v1'
GENE_ID = re.compile(r'^(ENSG[0-9]{11})(?:\.([0-9]+))?$')
REQUIRED = ('source_row_1based', 'source_gene_id', 'source_symbol', 'genome', 'feature_type')
OUTPUT = ('resource_id', 'source_file_sha256', *REQUIRED, 'canonical_gene_id',
          'canonical_gene_version', 'target_gene_symbol', 'target_chromosome',
          'hgnc_ids', 'approved_symbols', 'mapping_class', 'mapping_method',
          'candidate_gene_ids', 'primary_eligible', 'provisional_unique',
          'hold_reason', 'source_reference_status', 'source_mt_prefix',
          'common_reference_mt', 'rule_version')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'Missing TSV header: {path}')
        rows = []
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'Malformed TSV row {number}: {path}')
            rows.append(row)
        return rows, reader.fieldnames


class References:
    def __init__(self, genes, hgnc, withdrawn):
        self.genes = {}
        self.hgnc_by_id = {}
        self.by_ensembl = defaultdict(list)
        self.symbol_matches = defaultdict(list)
        self.withdrawn_symbols = defaultdict(list)
        for row in genes:
            gene = row['Gene stable ID']
            if not GENE_ID.fullmatch(gene) or '.' in gene:
                raise ValueError(f'Invalid target Ensembl gene ID: {gene}')
            if gene in self.genes:
                raise ValueError(f'Duplicate target gene ID: {gene}')
            full = GENE_ID.fullmatch(row['Gene stable ID version'])
            if not full or full.group(1) != gene:
                raise ValueError(f'Invalid target Ensembl version for {gene}')
            self.genes[gene] = row
        if not self.genes:
            raise ValueError('Empty Ensembl target axis')
        for row in hgnc:
            if row['status'] != 'Approved':
                continue
            if row['hgnc_id'] in self.hgnc_by_id:
                raise ValueError('Duplicate HGNC record ID')
            self.hgnc_by_id[row['hgnc_id']] = row
            ids = set(filter(None, row.get('ensembl_gene_id', '').split('|')))
            if any(not GENE_ID.fullmatch(gene) for gene in ids):
                raise ValueError(f'Unexpected HGNC Ensembl ID format: {ids}')
            row['_gene_ids'] = ids
            for gene in ids:
                self.by_ensembl[gene].append(row)
            for field, method in (('symbol', 'approved'), ('prev_symbol', 'previous'), ('alias_symbol', 'alias')):
                for symbol in filter(None, row.get(field, '').split('|')):
                    self.symbol_matches[symbol].append((row, method))
        for row in withdrawn:
            replacement_ids = set()
            for hgnc_id in re.findall(r'HGNC:[0-9]+', row.get('MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)', '')):
                if hgnc_id in self.hgnc_by_id:
                    replacement_ids.update(self.hgnc_by_id[hgnc_id]['_gene_ids'])
            row['_gene_ids'] = replacement_ids
            self.withdrawn_symbols[row['WITHDRAWN_SYMBOL']].append(row)

    @classmethod
    def load(cls, manifest_path):
        manifest_path = Path(manifest_path)
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('common_identity_target') != 'Homo sapiens GRCh38 Ensembl 116 (June 2026)':
            raise ValueError('Reference manifest is not the specified Ensembl-116 identity target')
        rows = {}
        for key in ('ensembl_registry', 'ensembl_genes', 'hgnc_complete', 'hgnc_withdrawn'):
            entry = manifest['assets'][key]
            path = manifest_path.parent / entry['path']
            if sha256(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
                raise ValueError(f'Frozen reference hash/size mismatch: {key}')
            if key == 'ensembl_registry':
                if 'ensembl_mart_116' not in path.read_text():
                    raise ValueError('Ensembl registry does not establish release 116')
            else:
                rows[key], _ = read_tsv(path)
        return cls(rows['ensembl_genes'], rows['hgnc_complete'], rows['hgnc_withdrawn']), manifest

    def symbol_candidates(self, symbol):
        matches = self.symbol_matches.get(symbol, [])
        withdrawn = self.withdrawn_symbols.get(symbol, [])
        ids = set()
        for record, _ in matches:
            ids.update(record['_gene_ids'])
        for record in withdrawn:
            ids.update(record['_gene_ids'])
        return ids, matches, withdrawn


def map_one(source, refs, resource, source_hash):
    result = {key: '' for key in OUTPUT}
    result.update(source)
    result.update(resource_id=resource, source_file_sha256=source_hash,
                  mapping_class='unmapped', primary_eligible='false', provisional_unique='false',
                  source_reference_status='unknown_not_reconstructed', rule_version=RULE_VERSION)
    symbol = source['source_symbol']
    source_id = source['source_gene_id']
    result['source_mt_prefix'] = str(symbol.upper().startswith('MT-')).lower()
    feature_type = source['feature_type']
    if feature_type and feature_type != 'Gene Expression':
        result.update(mapping_class='non_gene_feature', hold_reason='not_gene_expression')
        return result
    candidates, matches, withdrawn = refs.symbol_candidates(symbol)
    gene = None
    id_match = GENE_ID.fullmatch(source_id) if source_id else None
    if source_id:
        if not id_match:
            result['hold_reason'] = 'unsupported_source_identifier_no_symbol_fallback'
            return result
        gene = id_match.group(1)
        if gene not in refs.genes:
            result.update(mapping_class='retired_unresolved', mapping_method='source_id_absent_target_release',
                          candidate_gene_ids='|'.join(sorted(candidates)),
                          hold_reason='source_id_not_in_target_release_history_unavailable')
            return result
        if candidates and gene not in candidates:
            result.update(mapping_class='ambiguous', mapping_method='source_id_symbol_conflict',
                          candidate_gene_ids='|'.join(sorted(candidates | {gene})),
                          hold_reason='source_id_conflicts_with_frozen_symbol_mapping')
            return result
        result.update(mapping_class='resolved_unique', mapping_method='source_ensembl_in_target_release',
                      primary_eligible='true')
    else:
        hgnc_records = {row['hgnc_id'] for row, _ in matches}
        if len(candidates) > 1 or len(hgnc_records) > 1:
            result.update(mapping_class='ambiguous', mapping_method='symbol_multiple_records',
                          candidate_gene_ids='|'.join(sorted(candidates)), hold_reason='nonunique_symbol_or_history')
            return result
        if len(candidates) != 1:
            result.update(mapping_class='retired_unresolved' if withdrawn else 'unmapped',
                          hold_reason='no_unique_current_gene_id')
            return result
        gene = next(iter(candidates))
        if gene not in refs.genes:
            result.update(mapping_class='retired_unresolved', candidate_gene_ids=gene,
                          hold_reason='symbol_target_absent_ensembl116')
            return result
        methods = sorted({kind for _, kind in matches})
        if withdrawn:
            methods.append('withdrawn_history')
        result.update(mapping_class='resolved_unique', mapping_method='hgnc_' + '_'.join(methods),
                      provisional_unique='true', hold_reason='symbol_only_source_reference_not_reconstructed')
    target = refs.genes[gene]
    hgnc_rows = refs.by_ensembl.get(gene, [])
    result.update(canonical_gene_id=gene, candidate_gene_ids=gene,
                  canonical_gene_version=GENE_ID.fullmatch(target['Gene stable ID version']).group(2),
                  target_gene_symbol=target['Gene name'], target_chromosome=target['Chromosome/scaffold name'],
                  common_reference_mt=str(target['Chromosome/scaffold name'] in ('MT', 'chrM')).lower(),
                  hgnc_ids='|'.join(sorted(row['hgnc_id'] for row in hgnc_rows)),
                  approved_symbols='|'.join(sorted({row['symbol'] for row in hgnc_rows})))
    if source['genome'] and not (source['genome'].startswith('GRCh38') or source['genome'] == 'hg38'):
        result['primary_eligible'] = 'false'
        result['hold_reason'] = (result['hold_reason'] + ';' if result['hold_reason'] else '') + 'source_genome_requires_review'
    return result


def map_rows(rows, refs, resource, source_hash):
    expected = [str(i) for i in range(1, len(rows) + 1)]
    if [row['source_row_1based'] for row in rows] != expected:
        raise ValueError('Source rows must preserve complete contiguous original feature order starting at 1')
    if not rows or any(not row['source_gene_id'] and not row['source_symbol'] for row in rows):
        raise ValueError('Empty axis or feature lacking both ID and symbol')
    mapped = [map_one(row, refs, resource, source_hash) for row in rows]
    original_ids = Counter(row['source_gene_id'] for row in rows if row['source_gene_id'])
    canonical = Counter(row['canonical_gene_id'] for row in mapped if row['canonical_gene_id'])
    for row in mapped:
        reasons = []
        if row['source_gene_id'] and original_ids[row['source_gene_id']] > 1:
            reasons.append('duplicate_source_id')
        if row['canonical_gene_id'] and canonical[row['canonical_gene_id']] > 1:
            reasons.append('multiple_rows_same_canonical_no_aggregation')
        if reasons:
            row['primary_eligible'] = 'false'
            row['provisional_unique'] = 'false'
            row['hold_reason'] = ';'.join(filter(None, [row['hold_reason'], *reasons]))
    return mapped


def write_tsv(path, rows, fields):
    with Path(path).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--features', type=Path, required=True)
    parser.add_argument('--resource', required=True)
    parser.add_argument('--references', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    refs, manifest = References.load(args.references)
    rows, fields = read_tsv(args.features)
    if not set(REQUIRED).issubset(fields):
        raise ValueError(f'Feature TSV requires columns {REQUIRED}; found {fields}')
    if len(set(fields)) != len(fields):
        raise ValueError('Duplicate feature TSV column names')
    # Project only declared feature metadata; do not carry expression or arbitrary columns.
    rows = [{key: row[key] for key in REQUIRED} for row in rows]
    source_hash = sha256(args.features)
    mapped = map_rows(rows, refs, args.resource, source_hash)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / 'gene_mapping.tsv'
    summary_path = args.output_dir / 'mapping_summary.json'
    if output_path.exists() or summary_path.exists():
        raise FileExistsError('Refusing to overwrite existing mapping output')
    write_tsv(output_path, mapped, OUTPUT)
    summary = {
        'resource_id': args.resource, 'rule_version': RULE_VERSION, 'source_features': str(args.features),
        'source_file_sha256': source_hash, 'source_rows': len(rows),
        'mapping_class_counts': dict(Counter(row['mapping_class'] for row in mapped)),
        'primary_eligible_rows': sum(row['primary_eligible'] == 'true' for row in mapped),
        'provisional_unique_rows': sum(row['provisional_unique'] == 'true' for row in mapped),
        'duplicate_source_id_rows': sum('duplicate_source_id' in row['hold_reason'] for row in mapped),
        'duplicate_source_symbol_rows': sum(count - 1 for count in Counter(row['source_symbol'] for row in rows if row['source_symbol']).values()),
        'blank_source_symbol_rows': sum(not row['source_symbol'] for row in rows),
        'canonical_colliding_rows': sum('multiple_rows_same_canonical' in row['hold_reason'] for row in mapped),
        'reference_manifest_sha256': sha256(args.references),
        'reference_asset_sha256': {key: value['sha256'] for key, value in manifest['assets'].items()},
        'mapping_sha256': sha256(output_path), 'counts_read_or_changed': False,
        'count_conservation_check': 'Not applicable: mapping metadata only; no matrix projection performed.',
        'source_reference_status': 'Unknown: this mapper does not reconstruct original counting references.',
        'primary_eligible_boundary': 'Unique source ENSG in target release, without symbol/genome/row conflicts. Identity eligibility does not establish identical quantification or analytical cohort readiness.',
        'symbol_only_boundary': 'Unique nomenclature candidates remain provisional until source-reference crosswalk review; they are not strict-primary features.',
        'mitochondrial_boundary': 'Common-reference chromosome is an identity annotation; source MT-prefix is provisional QC and never overwritten.',
    }
    summary_path.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

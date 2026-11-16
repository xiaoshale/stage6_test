"""Narrow identity safeguards: ambiguity, dot preservation, IDs and collisions."""
import json
from pathlib import Path
from map_features import References, map_rows


ONE = 'ENSG00000000001'
TWO = 'ENSG00000000002'
THREE = 'ENSG00000000003'


def gene(identifier, symbol, version='4', chromosome='1'):
    return {'Gene stable ID': identifier, 'Gene stable ID version': identifier + '.' + version,
            'Gene name': symbol, 'Chromosome/scaffold name': chromosome, 'Gene type': 'protein_coding'}


def hgnc(identifier, symbol, gene_id, alias='', previous=''):
    return {'hgnc_id': identifier, 'symbol': symbol, 'status': 'Approved',
            'ensembl_gene_id': gene_id, 'alias_symbol': alias, 'prev_symbol': previous}


def source(identifier='', symbol='REAL1', feature_type='Gene Expression', genome='GRCh38'):
    return {'source_row_1based': '1', 'source_gene_id': identifier,
            'source_symbol': symbol, 'genome': genome, 'feature_type': feature_type}


def main():
    refs = References([gene(ONE, 'REAL1'), gene(TWO, 'RP11-34P13.3'), gene(THREE, 'SAME')],
                      [hgnc('HGNC:1', 'REAL1', ONE, 'SHARED', 'OLD1'),
                       hgnc('HGNC:2', 'RP11-34P13.3', TWO, 'SHARED'),
                       hgnc('HGNC:3', 'SAME', THREE)], [])
    passed = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        passed.append(name)

    def one(row):
        return map_rows([row], refs, 'fixture', 'fixture_hash')[0]

    row = one(source(ONE + '.2'))
    check('validated ENSG version split preserves full source and records target version',
          row['canonical_gene_id'] == ONE and row['source_gene_id'] == ONE + '.2'
          and row['canonical_gene_version'] == '4' and row['primary_eligible'] == 'true')
    row = one(source(symbol='RP11-34P13.3'))
    check('dotted source symbols are never version-stripped', row['source_symbol'] == 'RP11-34P13.3'
          and row['canonical_gene_id'] == TWO and row['provisional_unique'] == 'true')
    row = one(source(symbol='SHARED'))
    check('ambiguous alias is held without copying to two genes', row['mapping_class'] == 'ambiguous'
          and row['canonical_gene_id'] == '' and set(row['candidate_gene_ids'].split('|')) == {ONE, TWO})
    row = one(source('ENSG00000000009'))
    check('retired source ID is not silently remapped through current symbol', row['mapping_class'] == 'retired_unresolved'
          and row['canonical_gene_id'] == '' and row['primary_eligible'] == 'false')
    row = one(source(ONE, symbol='RP11-34P13.3'))
    check('source ID versus symbol conflict is held', row['mapping_class'] == 'ambiguous'
          and row['primary_eligible'] == 'false')
    row = one(source(symbol='REAL1'))
    check('unique symbol-only identity remains provisional', row['canonical_gene_id'] == ONE
          and row['provisional_unique'] == 'true' and row['primary_eligible'] == 'false')
    row = one(source(symbol='OLD1'))
    check('previous name is recorded as a different mapping method', row['mapping_method'] == 'hgnc_previous'
          and row['provisional_unique'] == 'true')
    row = one(source(ONE, feature_type='Antibody Capture'))
    check('non-gene features are held out of gene identity', row['mapping_class'] == 'non_gene_feature'
          and row['canonical_gene_id'] == '')
    row = one(source(ONE, genome='GRCh37'))
    check('different genome requires review without losing source ID', row['primary_eligible'] == 'false'
          and row['source_gene_id'] == ONE and 'source_genome_requires_review' in row['hold_reason'])
    a, b = source(ONE), source(ONE)
    b['source_row_1based'] = '2'
    rows = map_rows([a, b], refs, 'fixture', 'fixture_hash')
    check('duplicate source or canonical rows never automatically aggregate',
          len(rows) == 2 and all(r['primary_eligible'] == 'false' for r in rows)
          and all('multiple_rows_same_canonical' in r['hold_reason'] for r in rows))
    a, b = source(ONE), source(symbol='REAL1')
    b['source_row_1based'] = '2'
    rows = map_rows([a, b], refs, 'fixture', 'fixture_hash')
    check('distinct source rows mapping to one gene are both quarantined',
          all(r['primary_eligible'] == 'false' and r['provisional_unique'] == 'false' for r in rows))
    try:
        row = source(ONE)
        row['source_row_1based'] = '2'
        map_rows([row], refs, 'fixture', 'fixture_hash')
    except ValueError:
        check('missing or reordered original source rows fail', True)
    else:
        raise AssertionError('non-contiguous source rows accepted')
    result = {'checks_passed': len(passed), 'checks': passed,
              'boundary': 'Synthetic gene-identity safeguards, not real feature counts or matrix validation.'}
    (Path(__file__).resolve().parent / 'mapping_fixture_result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

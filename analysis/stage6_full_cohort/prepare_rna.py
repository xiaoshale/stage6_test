"""Make per-library raw integer-count checkpoints for the locked native intake."""
import argparse
import csv
import gc
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--component', choices=['GSE135893', 'GSE159354', 'GSE159585'], required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--pilot-code', type=Path, required=True)
    a = p.parse_args()
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        os.environ[name] = '4'
    sys.path.insert(0, str(a.pilot_code))
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import io
    from threadpoolctl import threadpool_limits
    from qc_spatial import read_counts

    resource = a.component
    units = json.loads((a.run / 'config/cohort_intake.json').read_text())['resources'][resource]['units']
    raw = a.root / 'data/raw' / resource
    out = a.run / 'rna_raw' / resource
    out.mkdir(exist_ok=False)
    axis_dir = a.run / 'genes/feature_axes'
    axis_dir.mkdir(parents=True, exist_ok=True)
    records = []

    def write(counts, obs, var, unit, provenance):
        counts = counts.tocsr()
        counts.sum_duplicates(); counts.eliminate_zeros(); counts.sort_indices()
        if counts.shape != (len(obs), len(var)) or not obs.index.is_unique or not var.index.is_unique:
            raise ValueError('Non-unique/inconsistent source axes')
        if not np.issubdtype(counts.data.dtype, np.integer) or (counts.data < 0).any():
            raise ValueError('Nonnegative integer counts required')
        lib = unit['library_or_section_id']
        obs = obs.copy()
        for key, value in {'donor_id': unit['donor_id'], 'diagnosis': unit['source_diagnosis'], 'analysis_group': unit['analysis_group'], 'library_id': lib, 'resource': resource, 'facs_fraction': unit.get('facs_fraction') or 'unknown', 'project_final_eligibility': 'not_adjudicated'}.items():
            obs[key] = value
        obs.index = pd.Index([f'{resource}::{lib}::{v}' for v in obs.index], name='cell_id')
        totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
        mt = var['source_symbol'].str.upper().str.startswith('MT-').to_numpy()
        mcounts = np.asarray(counts[:, mt].sum(axis=1)).ravel().astype(np.int64)
        obs['total_counts'] = totals
        obs['n_genes_by_counts'] = np.diff(counts.indptr)
        obs['pct_counts_mt'] = np.divide(100. * mcounts, totals, out=np.full(len(totals), np.nan), where=totals > 0) if mt.any() else np.nan
        var = var.copy(); var['mt'] = mt
        axis = var[['source_row_1based', 'source_gene_id', 'source_symbol', 'genome', 'feature_type']].to_csv(sep='\t', index=False)
        digest = hashlib.sha256(axis.encode()).hexdigest()
        apath = axis_dir / f'{resource}__{digest[:16]}.tsv'
        if apath.exists() and sha(apath) != digest:
            raise ValueError('Feature axis hash collision')
        if not apath.exists(): apath.write_text(axis)
        obj = ad.AnnData(X=counts, obs=obs, var=var)
        obj.uns['preparation'] = {'source_json': json.dumps(provenance), 'unit_json': json.dumps(unit), 'counts': 'Original counts, all source-called profiles/features retained; final review pending'}
        path = out / f'{lib}.h5ad'
        obj.write_h5ad(path, compression='gzip')
        restored = ad.read_h5ad(path)
        if restored.shape != obj.shape or not restored.obs_names.equals(obj.obs_names) or not restored.var_names.equals(obj.var_names) or (restored.X != counts).nnz:
            raise ValueError('Raw count/axis roundtrip failed')
        pd.testing.assert_frame_equal(restored.obs, obj.obs, check_dtype=False, check_categorical=False)
        pd.testing.assert_frame_equal(restored.var, obj.var, check_dtype=False, check_categorical=False)
        records.append({'path': str(path), 'resource': resource, 'library_id': lib, 'donor_id': unit['donor_id'], 'diagnosis': unit['source_diagnosis'], 'analysis_group': unit['analysis_group'], 'assay': 'snRNA' if resource == 'GSE159585' else 'scRNA', 'fraction': unit.get('facs_fraction'), 'source_label_columns': [c for c in obs if c.startswith('source_')], 'capture_identity_status': 'verified' if resource == 'GSE159354' else 'proxy_unresolved', 'cells': obj.n_obs, 'features': obj.n_vars, 'nnz': int(counts.nnz), 'original_count_sum': int(totals.sum()), 'checkpoint_sha256': sha(path), 'count_roundtrip': True, 'feature_axis': str(apath), 'feature_axis_sha256': digest})
        (out / 'inputs.json').write_text(json.dumps({'records': records}, indent=2) + '\n')
        print(resource, lib, obj.shape, 'raw count roundtrip PASS', flush=True)
        del restored, obj
        gc.collect()

    def symbols_var(symbols, genome):
        return pd.DataFrame({'source_row_1based': np.arange(1, len(symbols) + 1), 'source_gene_id': '', 'source_symbol': list(symbols), 'genome': genome, 'feature_type': 'Gene Expression'}, index=list(symbols))

    if resource == 'GSE135893':
        meta = pd.read_csv(raw / 'GSE135893_IPF_metadata.csv.gz', index_col=0, keep_default_na=False)
        symbols = pd.read_csv(raw / 'GSE135893_genes.tsv.gz', header=None, sep='\t', dtype=str)[0]
        barcodes = pd.Index(pd.read_csv(raw / 'GSE135893_barcodes.tsv.gz', header=None, sep='\t', dtype=str)[0])
        if not meta.index.is_unique or not barcodes.is_unique or not symbols.is_unique or set(meta['orig.ident']) != {u['library_or_section_id'] for u in units}:
            raise ValueError('Source/intake axes mismatch')
        with threadpool_limits(limits=4), gzip.open(raw / 'GSE135893_matrix.mtx.gz', 'rb') as f:
            matrix = io.mmread(f)
            if f.read(1): raise ValueError('Trailing MTX payload')
        if matrix.shape != (len(symbols), len(barcodes)) or not np.issubdtype(matrix.data.dtype, np.integer) or (matrix.data < 0).any():
            raise ValueError('Invalid full sparse MTX')
        if matrix.data.max() <= np.iinfo(np.int32).max: matrix.data = matrix.data.astype(np.int32)
        matrix = matrix.tocsc()
        var = symbols_var(symbols, 'hg38_source_release_unknown')
        provenance = {'source_sha256': sha(raw / 'GSE135893_matrix.mtx.gz'), 'calling': 'Exact author metadata only; unannotated columns excluded'}
        for unit in units:
            subset = meta.loc[meta['orig.ident'] == unit['library_or_section_id']]
            if subset.empty or set(subset['Sample_Name']) != {unit['donor_id']} or set(subset['Diagnosis']) != {unit['source_diagnosis']}:
                raise ValueError('Donor/diagnosis crosswalk mismatch')
            positions = barcodes.get_indexer(subset.index)
            if (positions < 0).any(): raise ValueError('Missing author barcode')
            counts = matrix[:, positions].T.tocsr()
            if not np.array_equal(np.asarray(counts.sum(axis=1)).ravel(), pd.to_numeric(subset['nCount_RNA'])) or not np.array_equal(np.diff(counts.indptr), pd.to_numeric(subset['nFeature_RNA'])):
                raise ValueError('Author count/feature metric disagreement')
            write(counts, subset.add_prefix('source_'), var, unit, provenance)

    elif resource == 'GSE159354':
        ledger = json.loads((a.run / 'config/donor_library_ledger.json').read_text())[resource]
        by_id = {r['library_id']: r for r in ledger['library_rows']}
        for unit in units:
            row = by_id[unit['library_or_section_id']]; url, = row['count_or_component_urls']
            path = raw / url.rsplit('/', 1)[1]
            acq = json.loads((a.root / 'data/acquisition_records' / f'{resource}__{path.name}.json').read_text())
            if acq['status'] != 'complete' or sha(path) != acq['sha256']: raise ValueError('Unverified source download')
            counts, barcodes, v = read_counts(path.read_bytes())
            var = pd.DataFrame({'source_row_1based': np.arange(1, len(v) + 1), 'source_gene_id': v['gene_id'].to_numpy(), 'source_symbol': v['gene_symbol'].to_numpy(), 'genome': v['genome'].to_numpy() if 'genome' in v else 'unknown', 'feature_type': v['feature_type'].to_numpy() if 'feature_type' in v else 'unknown'}, index=v.index)
            write(counts, pd.DataFrame({'source_barcode': barcodes}, index=barcodes), var, unit, {'source_sha256': acq['sha256'], 'url': url})
    else:
        from qc_snrna import stream_counts
        all_barcodes, meta = set(), {}
        donors = {u['donor_id'] for u in units}
        with gzip.open(raw / 'GSE159585_snRNAseq_metadata.csv.gz', 'rt') as f:
            for row in csv.DictReader(f, delimiter=';'):
                barcode = row['']
                if barcode in all_barcodes: raise ValueError('Duplicate nucleus barcode')
                all_barcodes.add(barcode)
                if row['patient_ID'] in donors: meta[barcode] = row
        matrix, symbols, barcodes, totals, detected, mt, gene_cells = stream_counts(raw / 'GSE159585_snRNAseq_raw_counts.csv.gz', all_barcodes, meta)
        matrix = matrix.tocsr()
        if len(set(symbols)) != len(symbols): raise ValueError('Duplicate nucleus gene labels')
        obs = pd.DataFrame([meta[b] for b in barcodes], index=barcodes).drop(columns=[''])
        var = symbols_var(symbols, 'source_reference_unknown')
        provenance = {'source_sha256': sha(raw / 'GSE159585_snRNAseq_raw_counts.csv.gz'), 'calling': 'Native deposited nuclei; COVID and bulk excluded'}
        for unit in units:
            mask = obs['sample_ID'].eq(unit['library_or_section_id']).to_numpy()
            if not mask.any() or set(obs.loc[mask, 'patient_ID']) != {unit['donor_id']}: raise ValueError('Nucleus crosswalk mismatch')
            write(matrix[mask], obs.loc[mask].add_prefix('source_'), var, unit, provenance)
    if len(records) != len(units): raise ValueError('Incomplete intake')
    (out / 'summary.json').write_text(json.dumps({'resource': resource, 'libraries': len(records), 'profiles': sum(r['cells'] for r in records), 'donors': len({r['donor_id'] for r in records}), 'all_count_roundtrips': True, 'final_eligibility': 'not_adjudicated'}, indent=2) + '\n')


if __name__ == '__main__':
    main()

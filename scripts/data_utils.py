"""Bounded-memory 10x HDF5 validation and donor aggregation helpers."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / 'data/raw/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5'
META = ROOT / 'data/raw/GSE174367_snRNA-seq_cell_meta.csv.gz'

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def decode(values):
    return np.array([v.decode('utf-8') for v in values])

def load_metadata(path=META):
    # Preserve the literal pathology category "None"; only NA/empty mean missing.
    frame = pd.read_csv(path, keep_default_na=False, na_values=['NA', ''])
    required = ['Barcode','SampleID','Diagnosis','Batch','Cell.Type','Age','Sex','PMI','RIN']
    if any(c not in frame for c in required):
        raise ValueError('Required metadata fields absent')
    if frame.Barcode.isna().any() or not frame.Barcode.is_unique:
        raise ValueError('Metadata barcodes missing or duplicated')
    # groupby drops missing keys; reject them before donor-level consistency checks.
    if frame.SampleID.isna().any() or frame.SampleID.astype(str).str.strip().eq('').any():
        raise ValueError('Donor identifiers missing or blank')
    if not frame.Diagnosis.isin(['AD', 'Control']).all():
        raise ValueError('Ambiguous diagnosis')
    for field in ['Diagnosis', 'Batch', 'Age', 'Sex', 'PMI', 'RIN']:
        if frame.groupby('SampleID')[field].nunique(dropna=False).max() > 1:
            raise ValueError(f'Inconsistent donor metadata: {field}')
    return frame

def align_barcodes(barcodes, meta):
    if len(set(barcodes)) != len(barcodes):
        raise ValueError('Duplicated matrix barcodes')
    locations = pd.Index(barcodes).get_indexer(meta.Barcode)
    if np.any(locations < 0):
        raise ValueError('Metadata nuclei absent from expression matrix; no silent drops allowed')
    return locations

def blocks(group, block_size=2048):
    """Yield raw gene-by-nucleus CSC blocks; never densify the full matrix."""
    shape = group['shape'][:]
    if shape.shape != (2,) or not np.issubdtype(shape.dtype, np.integer) or (shape < 0).any():
        raise ValueError('Invalid integer CSC shape')
    n_genes, n_cells = map(int, shape)
    pointers = group['indptr'][:]
    if (pointers.ndim != 1 or not np.issubdtype(pointers.dtype, np.integer)
            or len(pointers) != n_cells + 1 or pointers[0] != 0
            or np.any(pointers[1:] < pointers[:-1])):
        raise ValueError('Invalid CSC pointers')
    if pointers[-1] != len(group['data']) or len(group['indices']) != len(group['data']):
        raise ValueError('CSC storage length mismatch')
    for start in range(0, n_cells, block_size):
        stop = min(start + block_size, n_cells)
        lo, hi = int(pointers[start]), int(pointers[stop])
        data, indices = group['data'][lo:hi], group['indices'][lo:hi]
        if not np.issubdtype(data.dtype, np.integer) or (data < 0).any():
            raise ValueError('Counts are not nonnegative raw integers')
        if not np.issubdtype(indices.dtype, np.integer):
            raise ValueError('Gene indices must use integer storage')
        if (indices < 0).any() or (indices >= n_genes).any():
            raise ValueError('Gene index out of bounds')
        matrix = sparse.csc_matrix((data.astype(np.int64), indices, pointers[start:stop+1] - lo), shape=(n_genes, stop-start))
        # 10x files need not store row indices in ascending order. Sorting
        # changes representation only; duplicate coordinates remain a failure.
        matrix.sort_indices()
        if not matrix.has_canonical_format:
            raise ValueError('Duplicated entries within a nucleus')
        yield start, stop, matrix

def aggregate_block(matrix, labels, donor_ids):
    lookup = {d: i for i, d in enumerate(donor_ids)}
    selected = [(i, lookup[v]) for i, v in enumerate(labels) if v in lookup]
    if not selected:
        return np.zeros((matrix.shape[0], len(donor_ids)), dtype=np.int64)
    cells, donors = zip(*selected)
    indicator = sparse.csr_matrix((np.ones(len(cells), dtype=np.int64), (cells, donors)), shape=(matrix.shape[1], len(donor_ids)))
    return (matrix @ indicator).toarray()

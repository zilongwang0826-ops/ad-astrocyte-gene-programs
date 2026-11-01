#!/usr/bin/env python3
"""Sum raw astrocyte counts by donor and verify conservation."""
import json
import h5py
import numpy as np
import pandas as pd
from data_utils import ROOT, COUNTS, blocks, aggregate_block, write_json

audit=json.loads((ROOT/'results/audit/dataset_audit.json').read_text())
if audit['status']!='PASS': raise RuntimeError('Data gate is not PASS')
meta=pd.read_csv(ROOT/'results/audit/aligned_cell_metadata.tsv.gz',sep='\t',keep_default_na=False,na_values=['NA',''])
donors=pd.read_csv(ROOT/'results/audit/donor_audit.tsv',sep='\t',keep_default_na=False,na_values=['NA',''])
donors=donors[donors.included].sort_values(['Diagnosis','SampleID']).copy()
donor_ids=donors.SampleID.tolist()
assert len(donor_ids)==len(set(donor_ids))
astro=meta[(meta['Cell.Type']=='ASC') & meta.SampleID.isin(donor_ids)]
labels=np.full(audit['matrix_nuclei'],'',dtype=object)
labels[astro.matrix_column]=astro.SampleID
counts=np.zeros((audit['matrix_genes'],len(donor_ids)),dtype=np.int64)
with h5py.File(COUNTS) as h:
    for start,stop,matrix in blocks(h['matrix']):
        counts+=aggregate_block(matrix,labels[start:stop],donor_ids)
expected=astro.groupby('SampleID').raw_umi.sum().reindex(donor_ids).to_numpy()
assert np.array_equal(counts.sum(axis=0),expected)
assert int(counts.sum())==int(astro.raw_umi.sum())==audit['selected_astrocyte_umi']
assert (counts>=0).all()
genes=pd.read_csv(ROOT/'results/tables/gene_annotation.tsv',sep='\t')
frame=pd.DataFrame(counts,index=genes.gene_id,columns=donor_ids);frame.index.name='gene_id'
frame.to_csv(ROOT/'results/tables/pseudobulk_counts.tsv.gz',sep='\t')
donors=donors.rename(columns={'SampleID':'donor_id','Diagnosis':'diagnosis'})
donors['library_size']=counts.sum(axis=0)
donors.to_csv(ROOT/'results/tables/pseudobulk_samples.tsv',sep='\t',index=False)
write_json(ROOT/'results/audit/pseudobulk_checks.json',{'integer_nonnegative':True,'donor_names_aligned':True,'no_duplicated_donors':True,'donor_sum_conservation':True,'global_sum_conservation':True,'total_umi':int(counts.sum()),'dimensions':list(counts.shape)})
print('Pseudobulk checks passed:',counts.shape,'total UMIs',int(counts.sum()))

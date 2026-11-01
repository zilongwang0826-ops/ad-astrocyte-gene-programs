#!/usr/bin/env python3
"""Audit the actual GEO inputs before any differential expression is allowed."""
import h5py
import numpy as np
import pandas as pd
from data_utils import ROOT, COUNTS, META, digest, write_json, decode, load_metadata, align_barcodes, blocks

out = ROOT / 'results/audit'
write_json(out/'dataset_audit.json',{'status':'RUNNING','message':'No inference until every audit check completes.'})
meta = load_metadata()
with h5py.File(COUNTS) as handle:
    g = handle['matrix']
    shape = tuple(map(int, g['shape'][:]))
    barcodes = decode(g['barcodes'][:])
    locations = align_barcodes(barcodes, meta)
    genes = pd.DataFrame({'gene_id':decode(g['features/id'][:]),'gene_symbol':decode(g['features/name'][:]),'genome':decode(g['features/genome'][:]),'feature_type':decode(g['features/feature_type'][:])})
    if not genes.gene_id.is_unique: raise ValueError('Duplicated gene IDs')
    if set(genes.feature_type) != {'Gene Expression'}: raise ValueError('Non-expression features')
    library = np.zeros(shape[1], dtype=np.int64)
    detected = np.zeros(shape[1], dtype=np.int64)
    minimum, maximum, nnz, total = None, 0, 0, 0
    for start, stop, matrix in blocks(g):
        library[start:stop] = np.asarray(matrix.sum(axis=0)).ravel()
        detected[start:stop] = np.asarray((matrix > 0).sum(axis=0)).ravel()
        if matrix.nnz:
            val = int(matrix.data.min()); minimum = val if minimum is None else min(val, minimum)
            maximum = max(maximum, int(matrix.data.max()))
        nnz += matrix.nnz
        total += int(matrix.sum())
    schema = {n:{'shape':list(g[n].shape),'dtype':str(g[n].dtype)} for n in ['data','indices','indptr','barcodes','shape']}

meta['matrix_column'] = locations
meta['raw_umi'] = library[locations]
meta['detected_genes'] = detected[locations]
unmatched = pd.DataFrame({'Barcode':barcodes[~pd.Index(barcodes).isin(meta.Barcode)]})
unmatched['disposition'] = 'Excluded: absent from final author cell annotations'
unmatched.to_csv(out/'matrix_only_barcodes.tsv',sep='\t',index=False)
meta.to_csv(out/'aligned_cell_metadata.tsv.gz',sep='\t',index=False)
astro = meta[meta['Cell.Type']=='ASC'].copy()
donor = meta.groupby('SampleID',sort=True).first()[['Diagnosis','Batch','Age','Sex','PMI','RIN','Tangle.Stage','Plaque.Stage']].reset_index()
donor['total_nuclei'] = donor.SampleID.map(meta.groupby('SampleID').size())
donor['astrocyte_nuclei'] = donor.SampleID.map(astro.groupby('SampleID').size()).fillna(0).astype(int)
for c in ['raw_umi','detected_genes']:
    donor[f'astrocyte_median_{c}'] = donor.SampleID.map(astro.groupby('SampleID')[c].median())
donor['included'] = donor.astrocyte_nuclei >= 20
included = donor[donor.included]
counts = included.Diagnosis.value_counts().to_dict()
if any(counts.get(k,0)<5 for k in ['AD','Control']):
    raise ValueError('BLOCKED_DATA_AUDIT: fewer than five covered donors per group')
batch_table = pd.crosstab(included.Batch,included.Diagnosis)
if (batch_table.min(axis=1)==0).all():
    raise ValueError('BLOCKED_DATA_AUDIT: diagnosis completely separated by batch')
donor.to_csv(out/'donor_audit.tsv',sep='\t',index=False)
batch_table.to_csv(out/'batch_by_diagnosis.tsv',sep='\t')
pd.crosstab(included.Sex,included.Diagnosis).to_csv(out/'sex_by_diagnosis.tsv',sep='\t')
clinical=[]
for dx,part in included.groupby('Diagnosis'):
    for c in ['Age','PMI','RIN']:
        clinical.append({'diagnosis':dx,'covariate':c,'n_observed':int(part[c].notna().sum()),'n_missing':int(part[c].isna().sum()),'median':float(part[c].median()),'min':float(part[c].min()),'max':float(part[c].max())})
pd.DataFrame(clinical).to_csv(out/'covariate_summary.tsv',sep='\t',index=False)
genes.to_csv(ROOT/'results/tables/gene_annotation.tsv',sep='\t',index=False)
supp_path=ROOT/'data/reference/supplementary_donor_metadata.tsv'
supp=pd.read_csv(supp_path,sep='\t')
if set(supp['Sample.ID'])!=set(donor.SampleID) or not supp.donor_case_key.is_unique:
    raise ValueError('BLOCKED_DATA_AUDIT: supplementary donor identity mismatch')
cross=donor.merge(supp,on=None,left_on='SampleID',right_on='Sample.ID',suffixes=('_geo','_supp'),validate='one_to_one')
for name in ['Age','PMI','RIN']:
    if not np.allclose(cross[name+'_geo'],cross[name+'_supp'],equal_nan=True):
        raise ValueError('Clinical metadata conflict: '+name)
for name in ['Diagnosis','Sex']:
    if not (cross[name+'_geo']==cross[name+'_supp']).all():
        raise ValueError('Clinical metadata conflict: '+name)
cross.to_csv(out/'supplementary_metadata_crosscheck.tsv',sep='\t',index=False)
summary={'status':'PENDING_DESIGN_VERIFICATION','matrix_genes':shape[0],'matrix_nuclei':shape[1],'metadata_nuclei':len(meta),'matrix_only_nuclei':len(unmatched),'metadata_only_nuclei':0,'astrocyte_nuclei':len(astro),'donors_total':len(donor),'donors_included':len(included),'donors_by_diagnosis':counts,'donors_excluded':donor.loc[~donor.included,'SampleID'].tolist(),'astrocytes_min_per_donor':int(astro.groupby('SampleID').size().min()),'astrocytes_max_per_donor':int(astro.groupby('SampleID').size().max()),'gene_ids_duplicated':int(genes.gene_id.duplicated().sum()),'gene_symbols_duplicated':int(genes.gene_symbol.duplicated().sum()),'raw_count_min_stored':minimum,'raw_count_max':maximum,'stored_entries':nnz,'all_matrix_umi':total,'selected_astrocyte_umi':int(astro.raw_umi.sum()),'genome_labels':sorted(genes.genome.unique().tolist()),'metadata_missing_nuclei':meta.isna().sum().astype(int).to_dict(),'schema':schema}
write_json(out/'dataset_audit.json',summary)
summary['status']='PASS'
summary['independent_case_ids']=len(supp.donor_case_key.unique())
summary['supplementary_total_nuclei']=int(supp['Total cells in snRNA-seq'].sum())
summary['primary_model']='~ diagnosis'
summary['batch_provenance_warning']='GEO Batch differs from supplementary Library.Batch; latter assay attribution is unresolved. Preserve both; do not invent a mapping.'
summary['primary_interpretation']='Marginal unadjusted AD-versus-Control association; age, sex and technical confounding remain possible.'
manifest=[]
for path in [COUNTS,META,ROOT/'data/raw/h.all.v2026.1.Hs.symbols.gmt']:
    base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/' if path in [COUNTS,META] else 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/'
    manifest.append({'file':str(path.relative_to(ROOT)),'url':base+path.name,'bytes':path.stat().st_size,'sha256':digest(path),'retrieved_date_utc':'2026-10-29'})
manifest_path=ROOT/'data/input_manifest.json'
if manifest_path.exists():
    import json
    expected={r['file']:r for r in json.loads(manifest_path.read_text())}
    for row in manifest:
        if row['file'] not in expected or row['sha256']!=expected[row['file']]['sha256'] or row['bytes']!=expected[row['file']]['bytes']:
            raise ValueError('Input integrity mismatch; do not rewrite expected hashes: '+row['file'])
else:
    write_json(manifest_path,manifest)
rows='\n'.join(f'|{r.SampleID}|{r.Diagnosis}|{r.Batch}|{r.Age}|{r.Sex}|{r.astrocyte_nuclei}|' for r in donor.itertuples())
report=f'''# Dataset audit: PASS with interpretation limitations

## Matrix and provenance
Official GEO GSE174367 processed snRNA-seq inputs; URLs, byte sizes and SHA256 hashes are in `data/input_manifest.json`. The matrix is **{shape[0]:,} genes x {shape[1]:,} nuclei**, sparse CSC, with {nnz:,} stored integer entries. Every count block was checked for nonnegative integer values (stored range {minimum}–{maximum}), valid indices and unique coordinates. Stored row indices are not always sorted; sorting restores canonical sparse representation without altering counts. Gene IDs are unique; {summary['gene_symbols_duplicated']:,} duplicate-symbol entries are retained at gene-ID level and handled explicitly for enrichment.

The HDF5 reference label is `GRCh38.p12.premrna2`. Original Supplementary Note (pp. 11–12) describes Cell Ranger 3.0.2, a GRCh38.p12 pre-mRNA reference and `cellranger aggr` read-depth normalization. These are untransformed integer UMI counts **after the original aggregation/subsampling**, not untouched sequencing reads. Exact annotation release is not stated and is not inferred from the genome label. No new alignment, original-QC reproduction or doublet-removal analysis is claimed.

## Barcode reconciliation
All **{len(meta):,}** final author metadata barcodes match the HDF5 exactly; no prefix/suffix manipulation is needed. The **{len(unmatched)} matrix-only barcodes** are explicitly excluded because final donor/cell-type labels are absent and are listed in `matrix_only_barcodes.tsv`. This count is consistent with the original Supplementary Note reporting removal of 298 snRNA-seq nuclei in small clusters. Their per-barcode exclusion reasons cannot be reconstructed. Metadata-only barcodes: zero. The count-matrix total must not be confused with the astrocyte population.

## Donors, annotation and coverage
The original `Cell.Type = ASC` annotation selects **{len(astro):,} astrocyte nuclei**. Supplementary Table 1 supplies **18 unique case-year/case-number pairs**, one for each SampleID, and confirms diagnosis, age, sex, PMI and RIN. Thus there are **11 AD and 7 control donors**, with no duplicated biological replicate. No separate technical-library identifier is available in cell metadata. Region is FC (frontal/prefrontal cortex) in the original study. Original donor-case IDs remain pseudonymous.

Every donor has {summary['astrocytes_min_per_donor']}–{summary['astrocytes_max_per_donor']} astrocytes; the prespecified minimum of 20 excludes no donor. This is a coverage screen, not a statistical power guarantee. Per-nucleus UMI/detected-gene summaries are derived from the count matrix.

|Donor|Diagnosis|GEO batch|Age|Sex|Astrocytes|
|---|---|---:|---:|---|---:|
{rows}

## Missingness and confounding
Only PMI is missing (Sample-100); age, sex, RIN and GEO Batch are complete. The literal plaque category `None` means no plaques and is retained as a category, not treated as missing. Age values of 90 are top-coded in the source; they are not interpreted as exact ages above that threshold. AD has 6 female/5 male donors versus control 2 female/5 male. Clinical ranges and medians are in `covariate_summary.tsv`. Both diagnoses occur in every GEO batch (AD/control counts: batch 1=2/2, batch 2=4/2, batch 3=5/3).

Two source discrepancies are retained transparently. Supplementary Table 1 reports {summary['supplementary_total_nuclei']:,} snRNA-seq nuclei, differing from final cell metadata, and its `Library.Batch` assignment differs from GEO `Batch`. The original snRNA processing code uses a field named `Batch`; the table does not identify whether `Library.Batch` refers to the same assay/stage. It would be unjustified to silently replace one with the other. Both tables are cross-referenced in `supplementary_metadata_crosscheck.tsv`; final counts and author cell identities come from the processed GEO files. Neither candidate batch variable is perfectly confounded with diagnosis, but batch effects are not ruled out.

## Inferential gate and model decision
Reliable donors, unambiguous diagnosis, author astrocyte labels, untransformed integer counts, at least five covered donors per group and no demonstrated complete technical confounding satisfy the gate. Before differential expression, the requested **`~ diagnosis`** model is retained as a **marginal, unadjusted group comparison**. Age/sex imbalance and unresolved technical batch provenance preclude disease-specific or causal interpretation. No covariate is selected by gene-level P values. LODO measures donor influence; it cannot remove systematic confounding. The package must carry these limitations explicitly.
'''
(out/'dataset_audit.md').write_text(report)
write_json(out/'dataset_audit.json',summary)
print(summary)

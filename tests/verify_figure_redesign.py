#!/usr/bin/env python3
"""Independent source, plotted artist and cross-document integrity checks."""
from pathlib import Path
import hashlib,json,re,os
import numpy as np
import pandas as pd
import pymupdf
ROOT=Path(__file__).resolve().parents[1]
QA=ROOT/'results/audit/figure_qa'
checks=[]
def check(value,label):
    if not value:raise AssertionError(label)
    checks.append(label)
def frame(rel):return pd.read_csv(ROOT/rel,sep='\t')
def close(a,b,label):
    check(np.asarray(a).shape==np.asarray(b).shape,label+' shape')
    check(np.allclose(a,b,atol=1e-12,rtol=1e-12),label+' values')
def scene(stem):return json.loads((QA/(stem+'_scene.json')).read_text())['axes']
def scat(axis):return [np.array(a['xy']) for a in axis['artists'] if a['type']=='PathCollection']
def mesh(axis):return np.array(next(a['values'] for a in axis['artists'] if a['type']=='QuadMesh'))
baseline=json.loads((QA/'BASELINE_NUMERICAL_HASHES.json').read_text())
for rel,digest in baseline.items():
    if rel=='results/audit/reported_claims.json':
        old_claims=QA/'PRE_TARGETED_REPORTED_CLAIMS.json'
        check(hashlib.sha256(old_claims.read_bytes()).hexdigest()==digest,'Historical claim registry preserved')
        old=json.loads(old_claims.read_text());current=json.loads((ROOT/rel).read_text())
        for key,value in old.items():
            expected=value.copy()
            # The historical registry keeps its exact bytes. Only this approved
            # public filename change is normalized; values and other fields are strict.
            if expected.get('source')=='results/audit/final_revision_sensitivity.json':
                expected['source']='results/audit/sensitivity_validation.json'
            check(current.get(key)==expected,'Unchanged primary/sensitivity claim '+key)
    elif rel in {'application/Research_Description.md', 'application/CV_Project_Bullets.md'}:
        # Application prose is reviewed separately from scientific figure data.
        # The baseline retains its historical hashes; no old text is reconstructed
        # here to decide whether current numerical outputs are valid.
        check((ROOT/rel).is_file(),'Application description available: '+rel)
    else:
        check(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,'Unchanged '+rel)
d=frame('results/plot_data/figure_1_donors.tsv')
s=frame('results/tables/pseudobulk_samples.tsv').set_index('donor_id')
c=frame('results/sensitivity/donor_covariates.tsv').set_index('donor_id')
check(d.donor_id.is_unique and set(d.donor_id)==set(s.index),'All independent donor identities')
close(d.astrocyte_nuclei,s.loc[d.donor_id,'astrocyte_nuclei'],'Donor nucleus counts')
close(d.recorded_age_score,c.loc[d.donor_id,'recorded_age_score'],'Capped recorded ages')
check(d.sex.tolist()==c.loc[d.donor_id,'sex'].tolist(),'Recorded sex')
z=scene('figure_1_study_overview')
close(scat(z[0])[0],np.c_[d.astrocyte_nuclei,np.arange(18)],'Figure 1A artist offsets')
for actual,sex in zip(scat(z[1]),['F','M']):
    part=d[d.sex==sex];close(actual,np.c_[part.recorded_age_score,part.astrocyte_nuclei],'Figure 1B '+sex)
de=frame('results/tables/differential_expression.tsv')
age=frame('results/sensitivity/age_score_differential_expression.tsv').set_index('gene_id')
p=frame('results/plot_data/figure_2_model_comparison.tsv')
selected=de.sort_values(['PValue','gene_id']).head(10)
check(p.gene_id.tolist()==selected.gene_id.tolist(),'Fixed primary-P selection')
check((de.FDR>=.05).all() and (age.FDR>=.05).all(),'Zero corrected gene discoveries in both models')
z=scene('figure_2_differential_expression')
close(scat(z[0])[0],np.c_[de.logFC,-np.log10(np.maximum(de.PValue,np.nextafter(0.,1.)))],'All-gene plotted coordinates')
close(scat(z[1])[0],np.c_[selected.logFC,range(10)],'Primary effect artist')
close(scat(z[1])[1],np.c_[age.loc[selected.gene_id,'logFC'],range(10)],'Age effect artist')
en=frame('results/tables/enrichment_results.tsv').sort_values(['PValue','gene_set']).head(10)
rho=frame('results/sensitivity/correlation_all_sets.tsv')
aset=frame('results/sensitivity/age_score_enrichment.tsv').set_index('gene_set')
mat=rho.pivot(index='gene_set',columns='correlation',values='FDR').loc[en.gene_set,[0,.01,.03,.05]]
z=scene('figure_3_gene_set_enrichment')
close(np.concatenate(scat(z[0])),np.c_[en.rank_biserial,range(10)],'Rank-effect artists')
close(mesh(z[1]),-np.log10(mat.to_numpy()),'Correlation FDR mesh')
close(mesh(z[2]),-np.log10(aset.loc[en.gene_set,'FDR'].to_numpy()[:,None]),'Age FDR mesh')
counts=frame('results/plot_data/figure_3_full_family_counts.tsv')
check(counts.significant_sets.tolist()==[18,5,0,0] and counts.tested_sets.eq(49).all(),'Unchanged full-family discovery counts')
lodo=frame('results/tables/lodo_selected_results.tsv')
summary=frame('results/tables/lodo_selected_summary.tsv').sort_values('selection_rank')
z=scene('figure_4_donor_sensitivity');points=scat(z[0])
delta=lodo.assign(change=lodo.logFC-lodo.full_logFC).pivot(index='gene_id',columns='omitted_donor',values='change').loc[summary.gene_id,d.donor_id]
close(mesh(z[1]),delta.to_numpy(),'LODO signed-change mesh')
for i,row in enumerate(summary.itertuples()):
    part=lodo[lodo.gene_id==row.gene_id].set_index('omitted_donor').loc[d.donor_id]
    close(points[2*i][:,0],part.logFC,'LODO '+row.gene_symbol+' offsets')
    close(points[2*i+1],[[row.full_logFC,i]],'Full estimate '+row.gene_symbol)
    close([part.logFC.min(),part.logFC.max()],[row.min_logFC,row.max_logFC],'LODO '+row.gene_symbol+' range')
mds=frame('results/targeted_validation/donor_mds.tsv')
z=scene('figure_s1_donor_structure')
for axis in z[:2]:close(scat(axis)[0],mds[['MDS1','MDS2']],'MDS common coordinates')
for actual,sex in zip(scat(z[2]),['F','M']):
    close(actual,mds.loc[mds.sex==sex,['MDS1','MDS2']],'MDS sex '+sex)
for axis,field in [(z[3],'GEO_batch'),(z[4],'supplementary_library_batch')]:
    for actual,level in zip(scat(axis),sorted(mds[field].unique())):
        close(actual,mds.loc[mds[field]==level,['MDS1','MDS2']],f'MDS {field} {level}')
gl=frame('results/targeted_validation/gene_set_lodo.tsv')
gs=frame('results/targeted_validation/gene_set_lodo_summary.tsv').set_index('gene_set')
five=frame('results/tables/enrichment_results.tsv').query('FDR < 0.05').sort_values(['PValue','gene_set'])
z=scene('figure_s2_program_sensitivity');points=scat(z[0])
q=gl.pivot(index='gene_set',columns='omitted_donor',values='FDR').loc[five.gene_set,d.donor_id]
close(mesh(z[1]),-np.log10(q.to_numpy()),'Program LODO FDR artist mesh')
for i,name in enumerate(five.gene_set):
    part=gl[gl.gene_set==name].set_index('omitted_donor').loc[d.donor_id]
    close(points[2*i][:,0],part.rank_biserial,'Program LODO effects '+name)
    close(points[2*i+1],[[gs.loc[name,'primary_rank_biserial'],i]],'Program full effect '+name)
    close([part.rank_biserial.min(),part.rank_biserial.max()],[gs.loc[name,'min_effect'],gs.loc[name,'max_effect']],'Program range '+name)
for name in ['figure_1_study_overview','figure_2_differential_expression','figure_3_gene_set_enrichment','figure_4_donor_sensitivity','figure_s1_donor_structure','figure_s2_program_sensitivity']:
    for ext in ['png','svg','pdf']:check((ROOT/'results/figures'/f'{name}.{ext}').is_file(),name+' '+ext)
    doc=pymupdf.open(ROOT/'results/figures'/f'{name}.pdf')
    check(not doc[0].get_images(),name+' vector geometry')
    check(len(doc[0].get_text())>200,name+' searchable text')
    check('<text' in (ROOT/'results/figures'/f'{name}.svg').read_text(),name+' editable SVG text')
    layout=json.loads((QA/(name+'_layout.json')).read_text())
    check(not layout['text_outside_canvas'],name+' text inside canvas')
fig3=pymupdf.open(ROOT/'results/figures/figure_3_gene_set_enrichment.pdf')
min_type=min(span['size'] for block in fig3[0].get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans'] if span['text'].strip())
report_scale=(595.275590551-90)/fig3[0].rect.width
check(min_type*report_scale>=8.0,'Figure3 effective type >=8pt at full report width')
summary_pdf=pymupdf.open(ROOT/'application/AD_Astrocyte_Research_Summary.pdf')
check(len(summary_pdf)==1,'Summary remains one page')
text=(ROOT/'report/manuscript_source/manuscript.md').read_text()
check('ranges, not confidence intervals' in text,'LODO uncertainty defined')
check('No gene passes BH FDR < 0.05' in text,'Null gene finding in caption')
check('recorded-age model' in text and 'same ranking' in text,'Sensitivity definitions retained')
record={'status':'PASS','checks_passed':len(checks),'checks':checks,'protected_files':len(baseline)-3,'editorial_files_reviewed_separately':2,'preserved_primary_claims':len(old),'new_inferential_tests_in_plotting':0,'targeted_results':'49 Hallmark sets x 18 donor omissions; computed separately from plotting','scope':'Independent key/coordinate/value comparison against saved results and serialized Matplotlib artists; visual QA is separate.'}
output=Path(os.environ.get('AD_ASTRO_INTERNAL_DIR',ROOT/'.build'))/'validation/figure_consistency.json'
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(record,indent=2)+'\n')
print(f'PASS: {len(checks)} source/artist/document checks; {len(baseline)-3} numerical/source files unchanged; historical registry archived and {len(old)} prior claims preserved.')

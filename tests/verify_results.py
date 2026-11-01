#!/usr/bin/env python3
"""Independent output, claim, multiplicity and document consistency checks."""
import hashlib
import json
from pathlib import Path
import re
import sys
import os
import gzip
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
import pymupdf as fitz
ROOT=Path(__file__).resolve().parents[1]
checks=[]
def check(condition,label):
    if not bool(condition):raise AssertionError(label)
    checks.append(label)
def load(path):return json.loads((ROOT/path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
audit=load('results/audit/dataset_audit.json')
d=load('results/audit/analysis_summary.json');e=load('results/audit/enrichment_summary.json');l=load('results/audit/lodo_summary.json')
de=pd.read_csv(ROOT/'results/tables/differential_expression.tsv',sep='\t')
en=pd.read_csv(ROOT/'results/tables/enrichment_results.tsv',sep='\t')
samples=pd.read_csv(ROOT/'results/tables/pseudobulk_samples.tsv',sep='\t')
counts=pd.read_csv(ROOT/'results/tables/pseudobulk_counts.tsv.gz',sep='\t',index_col=0)
check(audit['status']=='PASS','Data gate PASS')
check(counts.columns.tolist()==samples.donor_id.tolist(),'Pseudobulk donor alignment')
check(samples.donor_id.is_unique and len(samples)==18,'No duplicated independent donor')
check(np.array_equal(counts.sum().values,samples.library_size.values),'Pseudobulk library-size conservation')
check(int(counts.to_numpy().sum())==audit['selected_astrocyte_umi'],'Astrocyte global sum conservation')
check(len(de)==d['tested_genes'] and de.gene_id.is_unique,'Complete unique tested-gene results')
check(np.allclose(de.FDR,false_discovery_control(de.PValue.to_numpy()),atol=1e-12,rtol=1e-12),'Independent SciPy BH check for genes')
check(np.allclose(en.FDR,false_discovery_control(en.PValue.to_numpy()),atol=1e-12,rtol=1e-12),'Independent SciPy BH check for sets')
check(int((de.FDR<.05).sum())==d['significant_genes'],'Gene significance count')
check(int((en.FDR<.05).sum())==e['significant_gene_sets'],'Set significance count')
check((de.contrast=='AD-Control').all() and d['reference_level']=='Control','Contrast orientation')
panel=pd.read_csv(ROOT/'results/tables/lodo_selected_results.tsv',sep='\t')
selected=pd.read_csv(ROOT/'results/tables/lodo_selected_genes.tsv',sep='\t')
check(selected.gene_id.tolist()==de.sort_values(['PValue','gene_id']).head(10).gene_id.tolist(),'Fixed full-data P-based display selection')
check(len(panel)==180 and panel.groupby('gene_id').omitted_donor.nunique().eq(18).all(),'All LODO panel fits present')
check(int(panel.groupby('gene_id').direction_matches_full.all().sum())==l['selected_all_directions_consistent'],'LODO reported sign consistency')
figures=list((ROOT/'results/figures').glob('*.pdf'))
main_figures=[p for p in figures if re.match(r'figure_[1-4]_',p.name)]
supplementary_figures=[p for p in figures if re.match(r'figure_s[12]_',p.name)]
check(len(main_figures)==4 and len(supplementary_figures)==2 and len(figures)==6,'Four main and two supplementary figures')
for f in figures:
    check(f.with_suffix('.png').exists() and f.with_suffix('.svg').exists(),f.stem+' exports')
    provenance=load('results/audit/figure_qa/'+f.stem+'_source.json')
    for name,value in provenance['sources'].items():
        source_path=ROOT/name if name.startswith('results/') else ROOT/'results/tables'/name
        check(sha(source_path)==value,f.stem+' source hash: '+name)
claims=load('results/audit/reported_claims.json')
targeted=load('results/targeted_validation/validation_summary.json')
check(targeted['status']=='PASS','Completed targeted validation')
programs=pd.read_csv(ROOT/'results/targeted_validation/gene_set_lodo.tsv',sep='\t')
program_summary=pd.read_csv(ROOT/'results/targeted_validation/gene_set_lodo_summary.tsv',sep='\t')
primary_names=set(en.loc[en.FDR<.05,'gene_set'])
five=program_summary[program_summary.gene_set.isin(primary_names)]
check(len(five)==5 and len(programs)==18*49,'Complete targeted validation output')
for donor,part in programs.groupby('omitted_donor'):
    check(set(part.gene_set)==set(en.gene_set),donor+' complete program family')
    check(np.allclose(part.FDR,false_discovery_control(part.PValue.to_numpy()),atol=1e-14,rtol=1e-12),donor+' independent program BH')
for key,value in {
    'program_direction_stable':int(five.direction_consistent_n.eq(18).sum()),
    'program_min_fdr_retained':int(five.FDR_below_n.min()),
    'program_max_fdr_retained':int(five.FDR_below_n.max()),
    'program_max_effect_change':float(five.max_absolute_effect_change.max()),
}.items():
    check(np.isclose(claims[key]['value'],value,atol=1e-14,rtol=1e-12),key+' complete-table source')
mds=load('results/targeted_validation/mds_summary.json')
for i in range(2):
    check(np.isclose(claims[f'mds{i+1}_variance']['value'],mds['variance_percent'][i]),f'MDS axis {i+1} variance source')
check(claims['mds_top_genes']['value']==mds['top_genes']==500,'Fixed MDS gene count')
template=(ROOT/'report/manuscript_source/manuscript_template.md').read_text()
text=(ROOT/'report/manuscript_source/manuscript.md').read_text()
cohort=pd.read_csv(ROOT/'results/audit/cohort_table_display.tsv',sep='\t',dtype=str)
cohort_md='| Characteristic | AD | Control |\n|---|---|---|\n'+'\n'.join('| '+r.Characteristic+' | '+r.AD+' | '+r.Control+' |' for r in cohort.itertuples(index=False))
expected=re.sub(r'\[\[([a-z0-9_]+)\]\]',lambda m:claims[m[1]]['rendered'],template.replace('[[COHORT_TABLE]]',cohort_md))
check(text==expected,'Every manuscript result placeholder resolves exactly from the claim registry')
check('[[' not in text and 'AUTHOR_INPUT_NEEDED' not in text,'No unfinished manuscript finding placeholders')
figure_captions=[int(n) for n in re.findall(r'^\*\*Figure (\d+)\s*\|',text,re.M)]
figure_references={int(n) for n in re.findall(r'\bFig\.\s*(\d+)\b',text)}
figure_images=[int(n) for n in re.findall(r'!\[Figure (\d+)\]',text)]
check(figure_captions==[1,2,3,4],'Figure captions numbered continuously from 1 to 4')
check(figure_references==set(figure_captions),'Every figure has a valid main-text cross-reference')
check(figure_images==figure_captions,'Figure image labels match caption order')
for number,target in re.findall(r'!\[Figure (\d+)\]\(([^)]+)\)',text):
    check(f'figure_{number}_' in target and (ROOT/'report/manuscript_source'/target).resolve().is_file(),f'Figure {number} resolves to the matching image file')
table_captions=[int(n) for n in re.findall(r'^\*\*Table (\d+)\s*\|',text,re.M)]
table_references={int(n) for n in re.findall(r'(?<!Supplementary )\bTable\s+(\d+)\b',text)}
check(table_captions==list(range(1,len(table_captions)+1)) and table_references==set(table_captions),'Main-table numbering and references agree; original-study Supplementary Table 1 is distinct')
check(claims['tested_genes']['value']==len(de) and claims['astrocytes']['value']==int(samples.astrocyte_nuclei.sum()),'Core reported numbers independently reconciled')
for key,field in [('emx2_logfc','logFC'),('emx2_f','F'),('emx2_p','PValue'),('emx2_fdr','FDR')]:
    check(np.isclose(claims[key]['value'],de.set_index('gene_symbol').loc['EMX2',field],rtol=1e-14),key+' source statistic')
for key,name,field in [('tnfa_fdr','TNFA_SIGNALING_VIA_NFKB','FDR'),('tnfa_effect','TNFA_SIGNALING_VIA_NFKB','rank_biserial'),('inflammatory_fdr','INFLAMMATORY_RESPONSE','FDR'),('upr_fdr','UNFOLDED_PROTEIN_RESPONSE','FDR'),('il6_fdr','IL6_JAK_STAT3_SIGNALING','FDR'),('bile_fdr','BILE_ACID_METABOLISM','FDR'),('bile_effect','BILE_ACID_METABOLISM','rank_biserial')]:
    check(np.isclose(claims[key]['value'],en.set_index('gene_set').loc['HALLMARK_'+name,field],rtol=1e-14),key+' source statistic')
check({p.name for p in (ROOT/'application').glob('*.pdf')}=={'AD_Astrocyte_Research_Report.pdf','AD_Astrocyte_Research_Summary.pdf'},'Only the two final application PDFs are distributed')
for file in ['AD_Astrocyte_Research_Report.pdf','AD_Astrocyte_Research_Summary.pdf']:
    doc=fitz.open(ROOT/'application'/file)
    check(len(doc)==1 if 'Summary' in file else 5<=len(doc)<=7,file+' requested page-count range')
    source_name='summary.md' if 'Summary' in file else 'manuscript.md'
    source_text=(ROOT/'report/manuscript_source'/source_name).read_text()
    body_title=next(line[2:].strip() for line in source_text.splitlines() if line.startswith('# '))
    check(doc.metadata['title']==body_title,file+' metadata title matches body heading')
    check(doc.metadata['author']=='Zilong Wang',file+' author metadata')
    check(doc.metadata['subject']=='Reanalysis of public snRNA-seq data (GSE174367)',file+' subject metadata')
    check(all('Zilong Wang · Project period: August 2026' in page.get_text() for page in doc),file+' requested footer on every page')
    contents='\n'.join(p.get_text() for p in doc)
    check(len(contents)>1500 and '12,331' in contents and '4,756' in contents,file+' searchable actual text')
    check('\ufffd' not in contents,file+' no replacement glyph')
    for index,page in enumerate(doc):
        for word in page.get_text('words'):
            check(word[0]>=-1 and word[1]>=-1 and word[2]<=page.rect.width+1 and word[3]<=page.rect.height+1,f'{file} page {index+1} word within page: {word[4]}')
required=['application/CV_Project_Bullets.md','application/Research_Description.md','README.md','LICENSE','run.sh','data/input_manifest.json','environment/requirements.lock.txt','environment/renv.lock','results/audit/dataset_audit.md']
for path in required:check((ROOT/path).is_file(),'Required deliverable '+path)

# New sensitivity families are independently reconstructed with SciPy BH.
sen=load('results/audit/sensitivity_validation.json')
scenarios=pd.read_csv(ROOT/'results/sensitivity/correlation_all_sets.tsv',sep='\t')
age_de=pd.read_csv(ROOT/'results/sensitivity/age_score_differential_expression.tsv',sep='\t')
age_sets=pd.read_csv(ROOT/'results/sensitivity/age_score_enrichment.tsv',sep='\t')
check(sorted(scenarios.correlation.unique())==[0,.01,.03,.05],'Fixed correlation scenarios')
for rho,part in scenarios.groupby('correlation'):
    check(len(part)==len(en) and set(part.gene_set)==set(en.gene_set),f'Fixed eligible family at correlation {rho}')
    check(np.allclose(part.FDR,false_discovery_control(part.PValue.to_numpy()),rtol=1e-12,atol=1e-14),f'Independent scenario BH at correlation {rho}')
    reported=next(r for r in sen['correlation_results'] if r['correlation']==rho)
    check(int((part.FDR<.05).sum())==reported['significant_sets'],f'Reported set count at correlation {rho}')
check(scenarios.groupby('gene_set').Direction.nunique().eq(1).all(),'Fixed-rank directions across correlation settings')
for name,frame in [('age genes',age_de),('age sets',age_sets)]:
    check(np.allclose(frame.FDR,false_discovery_control(frame.PValue.to_numpy()),rtol=1e-12,atol=1e-14),'Independent BH '+name)
check(set(age_de.gene_id)==set(de.gene_id),'Adjusted model preserves tested gene universe')
original_five=en.loc[en.FDR<.05,'gene_set']
retained=int((age_sets.set_index('gene_set').loc[original_five,'FDR']<.05).sum())
check(retained==sen['original_five_significant_after_age_score_adjustment']==claims['age_retained_original']['value'],'Original-five retention distinct from total adjusted count')
check(int((age_sets.FDR<.05).sum())==sen['adjusted_significant_sets']==claims['age_all_set_discoveries']['value'],'Total adjusted set count')
check(int((age_de.FDR<.05).sum())==claims['age_gene_discoveries']['value'],'Adjusted null gene-level claim')
emx=age_de.set_index('gene_symbol').loc['EMX2']
check(np.isclose(emx.logFC,claims['age_emx2_logfc']['value']) and np.isclose(emx.FDR,claims['age_emx2_fdr']['value']),'Adjusted EMX2 claim source')
check(sen['adjusted_model']['rank']==claims['age_design_rank']['value'] and sen['adjusted_model']['residual_df']==claims['age_residual_df']['value'],'Design-rank and nominal degrees-of-freedom claims')
for row in load('results/audit/reproducibility_comparison.json')['tables']:
    p=ROOT/row['table'];raw=gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes()
    check(hashlib.sha256(raw).hexdigest()==row['content_sha256'],'Preserved primary output '+row['table'])
artifact_paths=['report/manuscript_source/manuscript.md','report/manuscript_source/summary.md','report/manuscript_source/Supplementary_Methods_and_Tables.md','README.md','application/CV_Project_Bullets.md','application/Research_Description.md']
for name in artifact_paths:
    content=(ROOT/name).read_text()
    check(not re.search(r'AUTHOR[-_]TODO|\bTODO\b|\bFIXME\b|\bTBD\b|AUTHOR_INPUT_NEEDED',content),'No unfinished editorial marker in '+name)
    check('[[' not in content,'No unresolved numeric macro in '+name)
for block in ['1','2','3']:
    check('Supplementary Table S'+block in text and 'Supplementary Table S'+block in (ROOT/'report/manuscript_source/Supplementary_Methods_and_Tables.md').read_text(),'Supplementary table S'+block+' citation resolves')
refs=[int(x) for x in re.findall(r'^(\d+)\. ',text.split('## References')[1],re.M)]
cited={int(x) for group in re.findall(r'\[([\d,]+)\]',text.split('## References')[0]) for x in group.split(',')}
check(refs==list(range(1,12)) and cited==set(refs),'Eleven consecutive verified references all cited')
manual=[]
without=re.sub(r'\[\[[a-z0-9_]+\]\]','',template)
for m in re.finditer(r'(?<![A-Za-z])\d+(?:\.\d+)*(?![A-Za-z])',without):
    manual.append({'literal':m.group(),'context':without[max(0,m.start()-55):m.end()+55].replace('\n',' ')})
result={'status':'PASS','checks_passed':len(checks),'checks':checks,'result_claim_fields':len(claims),'manual_numeric_contexts_for_review':manual,'scope':'Programmatic consistency and page-boundary checks; visual inspection and scientific judgement are recorded separately.'}
check_output=Path(os.environ.get('AD_ASTRO_INTERNAL_DIR',ROOT/'.build'))/'validation/consistency_checks.json'
check_output.parent.mkdir(parents=True,exist_ok=True)
check_output.write_text(json.dumps(result,indent=2)+'\n')
print(f'PASS: {len(checks)} checks including document word boundaries; {len(claims)} evidence-linked numerical/parameter fields.')

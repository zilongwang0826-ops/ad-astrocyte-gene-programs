#!/usr/bin/env python3
"""Resolve claims from results and build searchable reports with vector figures."""
from pathlib import Path
import html
import json
import re
import tempfile
import os
import pandas as pd
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Flowable,KeepTogether,Table,TableStyle,CondPageBreak
from reportlab.lib.pagesizes import A4
from pypdf import PdfReader,PdfWriter,Transformation
from data_utils import ROOT,write_json,digest
from typography import register_report_fonts

SOURCE=ROOT/'report/manuscript_source'
APP=ROOT/'application'
def read_json(name): return json.loads((ROOT/'results/audit'/name).read_text())
audit=read_json('dataset_audit.json');de=read_json('analysis_summary.json')
enr=read_json('enrichment_summary.json');lodo=read_json('lodo_summary.json')
dt=pd.read_csv(ROOT/'results/tables/differential_expression.tsv',sep='\t').set_index('gene_symbol')
et=pd.read_csv(ROOT/'results/tables/enrichment_results.tsv',sep='\t').set_index('gene_set')
lt=pd.read_csv(ROOT/'results/tables/lodo_selected_summary.tsv',sep='\t').set_index('gene_symbol')
cov=pd.read_csv(ROOT/'results/audit/covariate_summary.tsv',sep='\t')
claims={}
def add(key,value,source,field,fmt=None):
    rendered=format(value,fmt) if fmt else str(value)
    claims[key]={'value':value,'rendered':rendered,'source':source,'field':field}
    return rendered
for key,field in [('astrocytes','astrocyte_nuclei'),('ad_donors','n_ad'),('control_donors','n_control'),('donors','donors'),('tested_genes','tested_genes'),('input_genes','input_genes')]:
    add(key,de[field],'results/audit/analysis_summary.json',field,',d')
for key,field in [('matrix_nuclei','matrix_nuclei'),('metadata_nuclei','metadata_nuclei'),('matrix_only','matrix_only_nuclei'),('min_astrocytes','astrocytes_min_per_donor'),('max_astrocytes','astrocytes_max_per_donor'),('astrocyte_umis','selected_astrocyte_umi')]:
    add(key,audit[field],'results/audit/dataset_audit.json',field,',d')
for key,field in [('tested_sets','tested_gene_sets'),('source_sets','source_gene_sets'),('significant_sets','significant_gene_sets'),('ranked_symbols','unique_symbols'),('min_set_size','min_mapped_size'),('max_set_size','max_mapped_size')]:
    add(key,enr[field],'results/audit/enrichment_summary.json',field,',d')
for key,field in [('displayed_genes','selected_genes'),('lodo_iterations','iterations')]:
    add(key,lodo[field],'results/audit/lodo_summary.json',field,'d')
add('lodo_donors',de['donors']-1,'results/audit/analysis_summary.json','donors minus one','d')
add('min_gene_fdr',de['min_fdr'],'results/audit/analysis_summary.json','min_fdr','.3f')
add('fdr_threshold',de['inference']['fdr_threshold'],'results/audit/analysis_summary.json','inference.fdr_threshold','.2f')
for label,field in [('emx2_logfc','logFC'),('emx2_f','F'),('emx2_p','PValue'),('emx2_fdr','FDR')]:
    add(label,float(dt.loc['EMX2',field]),'results/tables/differential_expression.tsv','EMX2.'+field,'.2e' if field=='PValue' else '.3f')
for key,gene,field in [('tnfa_effect','TNFA_SIGNALING_VIA_NFKB','rank_biserial'),('tnfa_fdr','TNFA_SIGNALING_VIA_NFKB','FDR'),('inflammatory_fdr','INFLAMMATORY_RESPONSE','FDR'),('upr_fdr','UNFOLDED_PROTEIN_RESPONSE','FDR'),('il6_fdr','IL6_JAK_STAT3_SIGNALING','FDR'),('bile_effect','BILE_ACID_METABOLISM','rank_biserial'),('bile_fdr','BILE_ACID_METABOLISM','FDR')]:
    add(key,float(et.loc['HALLMARK_'+gene,field]),'results/tables/enrichment_results.tsv',gene+'.'+field,'.4f' if field=='FDR' else '.3f')
for key,field in [('emx2_lodo_min','min_logFC'),('emx2_lodo_max','max_logFC')]:
    add(key,float(lt.loc['EMX2',field]),'results/tables/lodo_selected_summary.tsv','EMX2.'+field,'.3f')
add('max_lodo_change',lodo['selected_max_absolute_change'],'results/audit/lodo_summary.json','selected_max_absolute_change','.3f')
for dx,key in [('AD','ad_age_median'),('Control','control_age_median')]:
    add(key,float(cov[(cov.diagnosis==dx)&(cov.covariate=='Age')]['median'].iloc[0]),'results/audit/covariate_summary.tsv',dx+'.Age.median','.0f')
for key,field in [('r_version','R'),('edger_version','edgeR'),('limma_version','limma')]:
    add(key,de['versions'][field],'results/audit/analysis_summary.json','versions.'+field)
add('set_correlation',enr['inter_gene_correlation'],'results/audit/enrichment_summary.json','inter_gene_correlation','.2f')
add('hallmark_version','2026.1.Hs','results/audit/enrichment_summary.json','collection')
coverage=pd.read_csv(ROOT/'results/tables/enrichment_set_coverage.tsv',sep='\t').set_index('gene_set')
add('excluded_set_size',int(coverage.loc['HALLMARK_PANCREAS_BETA_CELLS','mapped_genes']),'results/tables/enrichment_set_coverage.tsv','HALLMARK_PANCREAS_BETA_CELLS.mapped_genes','d')
add('min_nuclei_rule',20,'results/audit/preanalysis_plan.md','coverage threshold','d')
add('min_donors_rule',5,'results/audit/preanalysis_plan.md','minimum group size','d')

sensitivity=read_json('sensitivity_validation.json')
correlation=pd.read_csv(ROOT/'results/sensitivity/correlation_original_five.tsv',sep='\t').set_index('gene_set')
age_sets=pd.read_csv(ROOT/'results/sensitivity/age_score_original_five.tsv',sep='\t').set_index('gene_set')
donors=pd.read_csv(ROOT/'results/sensitivity/donor_covariates.tsv',sep='\t')
for key,value in [('rho_zero',0.),('rho_mid',.03),('rho_high',.05)]:
    add(key,value,'results/sensitivity/sensitivity_plan.md','fixed correlation scenario','.0f' if value==0 else '.2f')
for key,rho in [('rho0_discoveries',0),('rho03_discoveries',.03),('rho05_discoveries',.05)]:
    row=next(r for r in sensitivity['correlation_results'] if r['correlation']==rho)
    add(key,row['significant_sets'],'results/sensitivity/correlation_summary.tsv',f'correlation={rho}.significant_sets','d')
for key,field in [('age_retained_original','original_five_significant_after_age_score_adjustment'),('age_gene_discoveries','adjusted_significant_genes'),('age_all_set_discoveries','adjusted_significant_sets')]:
    add(key,sensitivity[field],'results/audit/sensitivity_validation.json',field,'d')
for key,field,fmt in [('age_design_rank','rank','d'),('age_residual_df','residual_df','d'),('age_vif','predictor_VIF','.3f')]:
    add(key,sensitivity['adjusted_model'][field],'results/sensitivity/design_preflight.json',field,fmt)
for key,field in [('age_emx2_logfc','logFC'),('age_emx2_fdr','FDR')]:
    add(key,sensitivity['emx2_model_comparison'][1][field],'results/sensitivity/emx2_model_comparison.tsv','recorded-age-score.'+field,'.3f')
for key,column in [('rho03_original_fdr','FDR_rho_0.03'),('rho05_original_fdr','FDR_rho_0.05')]:
    values=correlation[column].to_numpy()
    assert max(values)-min(values)<1e-12
    add(key,float(values[0]),'results/sensitivity/correlation_original_five.tsv','all original five.'+column,'.3f')
for key,name in [('age_retained_fdr','TNFA_SIGNALING_VIA_NFKB'),('age_il6_fdr','IL6_JAK_STAT3_SIGNALING'),('age_bile_fdr','BILE_ACID_METABOLISM')]:
    add(key,float(age_sets.loc['HALLMARK_'+name,'adjusted_FDR']),'results/sensitivity/age_score_original_five.tsv',name+'.adjusted_FDR','.4f')
for key,dx,sex in [('ad_female','AD','F'),('ad_male','AD','M'),('control_female','Control','F'),('control_male','Control','M')]:
    add(key,int(((donors.diagnosis==dx)&(donors.sex==sex)).sum()),'results/sensitivity/donor_covariates.tsv',dx+'.'+sex+' count','d')

targeted=ROOT/'results/targeted_validation'
program_summary=pd.read_csv(targeted/'gene_set_lodo_summary.tsv',sep='\t').set_index('gene_set')
program_five=program_summary.loc[et.index[et.FDR<.05]].sort_values('primary_rank')
mds_summary=json.loads((targeted/'mds_summary.json').read_text())
assert len(program_five)==5 and program_five.n_iterations.eq(18).all()
for key,value,field,fmt in [
    ('program_direction_stable',int((program_five.direction_consistent_n==program_five.n_iterations).sum()),'original five.direction_consistent_n == n_iterations','d'),
    ('program_min_fdr_retained',int(program_five.FDR_below_n.min()),'original five.min(FDR_below_n)','d'),
    ('program_max_fdr_retained',int(program_five.FDR_below_n.max()),'original five.max(FDR_below_n)','d'),
    ('program_max_effect_change',float(program_five.max_absolute_effect_change.max()),'original five.max(max_absolute_effect_change)','.3f')]:
    add(key,value,'results/targeted_validation/gene_set_lodo_summary.tsv',field,fmt)
for i in [0,1]:
    add('mds'+str(i+1)+'_variance',float(mds_summary['variance_percent'][i]),'results/targeted_validation/mds_summary.json',f'variance_percent[{i}]','.1f')
add('mds_top_genes',int(mds_summary['top_genes']),'results/targeted_validation/mds_summary.json','top_genes','d')

cohort=[]
for dx in ['AD','Control']:
    part=donors[donors.diagnosis==dx]
    age=part.recorded_age_score
    values={'Donors':str(len(part)),'Recorded age, years':f'{age.median():.0f} ({age.min():.0f}–90+)',
       'Female / male':f'{int((part.sex=="F").sum())} / {int((part.sex=="M").sum())}',
       'PMI (deposited)':f'{part.PMI.median():.3f} ({part.PMI.min():.2f}–{part.PMI.max():.2f})',
       'RIN':f'{part.RIN.median():.1f} ({part.RIN.min():.1f}–{part.RIN.max():.1f})',
       'GEO batch 1 / 2 / 3':' / '.join(str(int((part.GEO_batch==i).sum())) for i in [1,2,3]),
       'Supplementary batch 1 / 2 / 3':' / '.join(str(int((part.supplementary_library_batch==i).sum())) for i in [1,2,3])}
    cohort.append(values)
cohort_rows=[{'Characteristic':k,'AD':cohort[0][k],'Control':cohort[1][k]} for k in cohort[0]]
pd.DataFrame(cohort_rows).to_csv(ROOT/'results/audit/cohort_table_display.tsv',sep='\t',index=False)
cohort_markdown='| Characteristic | AD | Control |\n|---|---|---|\n'+'\n'.join('| '+r['Characteristic']+' | '+r['AD']+' | '+r['Control']+' |' for r in cohort_rows)

def resolve(text):
    def sub(m):
        if m[1] not in claims: raise ValueError('Unresolved claim '+m[1])
        return claims[m[1]]['rendered']
    return re.sub(r'\[\[([a-z0-9_]+)\]\]',sub,text.replace('[[COHORT_TABLE]]',cohort_markdown))

manuscript=resolve((SOURCE/'manuscript_template.md').read_text())
(SOURCE/'manuscript.md').write_text(manuscript)
write_json(ROOT/'results/audit/reported_claims.json',claims)

def markdown_table(frame):
    def cell(value):
        if pd.isna(value): return 'NA'
        if isinstance(value,float): return f'{value:.6g}'
        return str(value)
    rows=['| '+' | '.join(map(str,frame.columns))+' |',
          '| '+' | '.join('---' for _ in frame.columns)+' |']
    rows.extend('| '+' | '.join(cell(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None))
    return '\n'.join(rows)

supplement=(SOURCE/'supplementary_template.md').read_text()
donor_display=donors[['donor_id','diagnosis','age_display','sex','PMI','RIN','astrocyte_nuclei','GEO_batch','supplementary_library_batch']].copy()
donor_display.columns=['Donor','Diagnosis','Age','Sex','PMI','RIN','Astrocyte nuclei','GEO batch','Supplementary batch']
cor_display=correlation.reset_index().drop(columns=['rank_biserial'])
cor_display['gene_set']=cor_display.gene_set.str.removeprefix('HALLMARK_')
age_display=age_sets.reset_index()[['gene_set','primary_direction','primary_FDR','adjusted_direction','adjusted_FDR']]
age_display['gene_set']=age_display.gene_set.str.removeprefix('HALLMARK_')
emx_display=pd.read_csv(ROOT/'results/sensitivity/emx2_model_comparison.tsv',sep='\t').drop(columns=['gene_id'])
for token,frame in [('[[DONOR_TABLE]]',donor_display),('[[CORRELATION_COUNTS]]',pd.read_csv(ROOT/'results/sensitivity/correlation_summary.tsv',sep='\t')),
                    ('[[CORRELATION_TABLE]]',cor_display),('[[AGE_TABLE]]',age_display),('[[EMX2_TABLE]]',emx_display)]:
    supplement=supplement.replace(token,markdown_table(frame))
program_display=program_five.reset_index()[['gene_set','primary_rank_biserial','min_effect','max_effect','direction_consistent_n','FDR_below_n','min_rank','max_rank']].copy()
program_display['gene_set']=program_display.gene_set.str.removeprefix('HALLMARK_')
program_display.columns=['Hallmark set','Full effect','Min effect','Max effect','Same direction / 18','FDR < 0.05 / 18','Best rank','Worst rank']
supplement=supplement.replace('[[PROGRAM_LODO_TABLE]]',markdown_table(program_display))
(SOURCE/'Supplementary_Methods_and_Tables.md').write_text(resolve(supplement))

summary_template=(SOURCE/'summary_template.md').read_text()
summary=resolve(summary_template)
(SOURCE/'summary.md').write_text(summary)

# Application prose and CV bullets require their own review; rebuilding the
# research report must not silently rewrite those separate deliverables.

# Resolve text without producing a PDF when an editorial pass is in progress.
import sys
if '--sources-only' in sys.argv:
    print('Manuscript, summary and supplementary text resolved from saved results.')
    raise SystemExit(0)

# Fonts are installed locally; the pipeline never downloads or redistributes them.
font_dir=register_report_fonts()
styles={
 'body':ParagraphStyle('body',allowWidows=0,allowOrphans=0,fontName='Body',fontSize=10,leading=13.2,spaceAfter=5,textColor=colors.HexColor('#202020')),
 'title':ParagraphStyle('title',fontName='BodyBold',fontSize=18,leading=21,spaceAfter=10,textColor=colors.black),
 'h2':ParagraphStyle('h2',fontName='BodyBold',fontSize=12,leading=15,spaceBefore=7,spaceAfter=5,keepWithNext=True,textColor=colors.black),
 'legend':ParagraphStyle('legend',fontName='Sans',fontSize=8.5,leading=11,spaceAfter=9),
 'ref':ParagraphStyle('ref',fontName='Body',fontSize=9.1,leading=11.5,spaceAfter=5),
 'h3':ParagraphStyle('h3',fontName='SansBold',fontSize=10.5,leading=13,spaceBefore=8,spaceAfter=6,keepWithNext=True),
 'table':ParagraphStyle('table',fontName='Body',fontSize=9.2,leading=11.5),
 'meta':ParagraphStyle('meta',fontName='Sans',fontSize=9,leading=13,spaceAfter=8),
}

def markup(s):
    s=html.escape(s)
    s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
    s=re.sub(r'\*(.+?)\*',r'<i>\1</i>',s)
    s=re.sub(r'(https://[^\s<]+)',r'<link href="\1" color="#303030">\1</link>',s)
    return s

placements=[]
class VectorFigure(Flowable):
    def __init__(self,path,width):
        super().__init__();self.path=path
        page=PdfReader(path).pages[0];self.origw=float(page.mediabox.width);self.origh=float(page.mediabox.height)
        self.width=width;self.height=width*self.origh/self.origw;self.hAlign='CENTER'
    def draw(self):
        m=self.canv._currentMatrix
        placements.append({'page':self.canv.getPageNumber()-1,'path':self.path,'x':m[4],'y':m[5],'scale':self.width/self.origw})

def footer(canvas,doc):
    canvas.setFont('Sans',7.4);canvas.setFillColor(colors.HexColor('#666666'))
    canvas.drawString(45,27,'Zilong Wang · Project period: August 2026')
    canvas.drawRightString(A4[0]-45,27,str(doc.page))

def build(text,target,is_summary=False):
    placements.clear();story=[];section='';width=A4[0]-90
    body_title=next(line[2:].strip() for line in text.splitlines() if line.startswith('# '))
    paragraphs=[p.strip() for p in re.split(r'\n\s*\n',text.strip())]
    i=0
    while i<len(paragraphs):
        p=paragraphs[i];i+=1
        if p=='<!-- PAGEBREAK -->':story.append(PageBreak());continue
        if p=='## References':
            # Keep the short bibliography together instead of stranding a few
            # citations on an otherwise empty continuation page.
            bibliography=ParagraphStyle('bibliography',parent=styles['ref'],fontSize=8.3,leading=10.5,spaceAfter=4)
            items=[Paragraph(markup(item.replace('\n',' ')),bibliography) for item in paragraphs[i:]]
            midpoint=(len(items)+1)//2
            table=Table([[items[:midpoint],'',items[midpoint:]]],colWidths=[(width-18)/2,18,(width-18)/2],hAlign='LEFT')
            table.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),0),('RIGHTPADDING',(0,0),(-1,-1),0),('TOPPADDING',(0,0),(-1,-1),0),('BOTTOMPADDING',(0,0),(-1,-1),0)]))
            story.append(KeepTogether([Paragraph('References',styles['h2']),table]));break
        if p.startswith('!['):
            stem=re.search(r'figure_[^)]+',p).group(0).replace('.png','.pdf')
            path=ROOT/'results/figures'/stem
            # Use the full content width in both documents so small figure type
            # is not reduced a second time in the one-page summary.
            group=[VectorFigure(path,width),Spacer(1,3)]
            if i<len(paragraphs) and paragraphs[i].startswith('**Figure'):
                group.append(Paragraph(markup(paragraphs[i].replace('\n',' ')),styles['legend']));i+=1
            story.append(KeepTogether(group));continue
        if p.startswith('| '):
            rows=[[cell.strip() for cell in row.strip().strip('|').split('|')] for row in p.splitlines() if not re.match(r'^\|[-: |]+\|$',row)]
            cells=[[Paragraph(markup(cell),styles['table']) for cell in row] for row in rows]
            table=Table(cells,colWidths=[width*.48,width*.26,width*.26],repeatRows=1,hAlign='LEFT')
            table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.white),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#555555')),('LINEBELOW',(0,-1),(-1,-1),.6,colors.HexColor('#555555')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]))
            group=[table,Spacer(1,5)]
            if i<len(paragraphs) and paragraphs[i].startswith('**Table'):
                group.append(Paragraph(markup(paragraphs[i]),styles['legend']));i+=1
            story.append(KeepTogether(group));continue
        if p.startswith('# '):style=styles['title'];p=p[2:]
        elif p.startswith('## '):style=styles['h2'];p=p[3:];section=p
        elif p.startswith('### '):style=styles['h3'];p=p[4:]
        elif p.startswith('**Figure ') or p.startswith('**Reading the figure.'):style=styles['legend']
        elif section=='References':style=styles['ref']
        elif p.startswith('Université') or p.startswith('Research project ') or p=='**Zilong Wang**':style=styles['meta']
        else:style=styles['body']
        if is_summary:
            if style is styles['body']:style=ParagraphStyle('summary_body',parent=style,fontSize=9.4,leading=12,spaceAfter=6)
            if style is styles['title']:style=ParagraphStyle('summary_title',parent=style,fontSize=15,leading=18,spaceAfter=8)
            if style is styles['h2']:style=ParagraphStyle('summary_h2',parent=style,fontSize=10.3,leading=13,spaceBefore=4,spaceAfter=5)
        story.append(Paragraph(markup(p.replace('\n',' ')),style))
    tmp=target.with_suffix('.text-only.pdf')
    doc=SimpleDocTemplate(str(tmp),pagesize=A4,rightMargin=45,leftMargin=45,topMargin=39,bottomMargin=42,
        title=body_title,author='Zilong Wang')
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    reader=PdfReader(tmp);writer=PdfWriter()
    for page in reader.pages:writer.add_page(page)
    for item in placements:
        source=PdfReader(item['path']).pages[0]
        transform=Transformation().scale(item['scale']).translate(item['x'],item['y'])
        writer.pages[item['page']].merge_transformed_page(source,transform,over=True)
    writer.add_metadata({'/Title':body_title,'/Author':'Zilong Wang','/Subject':'Reanalysis of public snRNA-seq data (GSE174367)'})
    with target.open('wb') as handle:writer.write(handle)
    # Intermediate is task-created and retained internally for audit, not in public package.
    internal=Path(os.environ.get('AD_ASTRO_INTERNAL_DIR',ROOT/'.build'))/'pdf-preview'
    internal.mkdir(parents=True,exist_ok=True)
    tmp.replace(internal/tmp.name)
    return len(writer.pages)

counts={}
final_report=APP/'AD_Astrocyte_Research_Report.pdf'
final_summary=APP/'AD_Astrocyte_Research_Summary.pdf'
counts['report_pages']=build(manuscript,final_report)
counts['summary_pages']=build(summary,final_summary,True)
assert counts['summary_pages']==1,'The research summary must remain one page.'
counts['final_report']=str(final_report.relative_to(ROOT))
counts['final_summary']=str(final_summary.relative_to(ROOT))
counts['fonts']={'body':'Times New Roman','headings':'Times New Roman (title and sections); Arial (subsections)','source':'Installed local TrueType files; no remote fonts'}
write_json(ROOT/'results/audit/document_build.json',counts)
print(counts)

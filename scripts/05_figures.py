#!/usr/bin/env python3
"""Four main and two supplementary plates from saved donor-level results."""
from pathlib import Path
import json, os, platform
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'.build/matplotlib'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from data_utils import ROOT,digest,write_json
from typography import register_plot_fonts
register_plot_fonts()
CFG=json.loads((ROOT/'scripts/figure_style.json').read_text())
AD,CONTROL,DARK,GRAY=[CFG['colors'][k] for k in ['AD','Control','ink','neutral']]
OUT=ROOT/'results/figures';DATA=ROOT/'results/plot_data';QA=ROOT/'results/audit/figure_qa'
for p in [OUT,DATA,QA]:p.mkdir(exist_ok=True,parents=True)
FDR=json.loads((ROOT/'results/audit/analysis_summary.json').read_text())['inference']['fdr_threshold']
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial'],'font.size':9,'axes.labelsize':9,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.65,'text.color':DARK,'axes.labelcolor':DARK,'pdf.fonttype':42,'ps.fonttype':42,
 'svg.fonttype':'none','svg.hashsalt':'astrocyte-v2','savefig.dpi':400})
sources={}
def load(rel):
    sources[rel]=digest(ROOT/rel)
    return pd.read_csv(ROOT/rel,sep='\t')
def data(name,frame):
    p=DATA/(name+'.tsv');frame.to_csv(p,sep='\t',index=False,float_format='%.17g');return str(p.relative_to(ROOT))
def panel(fig,x,text,letter):
    fig.text(x,.955,letter,weight='bold',fontsize=12,va='top')
    fig.text(x+.035,.952,text,fontsize=10,va='top')
def handles():
    return [Line2D([],[],marker='o',ls='',color=AD,label='AD'),Line2D([],[],marker='o',ls='',color=CONTROL,label='Control')]
def save(fig,name,inputs,derived,meaning):
    for axis in fig.axes:
        for artist in axis.collections:artist.set_rasterized(False)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[]
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1>fig.bbox.width+1 or b.y1>fig.bbox.height+1:bad.append(t.get_text())
    if bad:raise ValueError((name,bad))
    for ext in ['pdf','svg','png']:
        meta={'Title':meaning,'Author':'Zilong Wang','CreationDate':None,'ModDate':None} if ext=='pdf' else {'Date':None} if ext=='svg' else None
        fig.savefig(OUT/f'{name}.{ext}',dpi=400,metadata=meta)
    scene=[]
    for axis in fig.axes:
        artists=[]
        for item in axis.collections:
            kind=type(item).__name__
            if kind=='PathCollection':
                artists.append({'type':kind,'xy':np.asarray(item.get_offsets(),dtype=float).tolist()})
            elif kind=='QuadMesh':
                artists.append({'type':kind,'values':np.asarray(item.get_array(),dtype=float).tolist()})
        scene.append({'xlim':list(axis.get_xlim()),'ylim':list(axis.get_ylim()),'xlabel':axis.get_xlabel(),'ylabel':axis.get_ylabel(),'artists':artists})
    write_json(QA/f'{name}_scene.json',{'axes':scene})
    write_json(QA/f'{name}_layout.json',{'status':'PASS','text_outside_canvas':bad,'size_inches':list(fig.get_size_inches()),'axes_bounds':[list(a.get_position().bounds) for a in fig.axes]})
    write_json(QA/f'{name}_source.json',{'figure':name,'meaning':meaning,'sources':{s:sources[s] for s in inputs},
      'plotting_data':{s:digest(ROOT/s) for s in derived},'points_are_empirical':True,'command':'python scripts/05_figures.py',
      'configuration':'scripts/figure_style.json','configuration_sha256':digest(ROOT/'scripts/figure_style.json'),
      'versions':{'Python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__}})
    plt.close(fig)

# Individual donors, with capped recorded age shown separately from coverage.
sp='results/tables/pseudobulk_samples.tsv';cp='results/sensitivity/donor_covariates.tsv'
s=load(sp);c=load(cp)
d=s[['donor_id','diagnosis','astrocyte_nuclei']].merge(c[['donor_id','recorded_age_score','age_display','age_top_coded','sex']],on='donor_id',validate='one_to_one')
d['group_order']=d.diagnosis.map({'AD':0,'Control':1});d['id_order']=d.donor_id.str.removeprefix('Sample-').astype(int)
d=d.sort_values(['group_order','id_order']).reset_index(drop=True)
assert len(d)==18 and d.astrocyte_nuclei.sum()==4756
fig=plt.figure(figsize=(8.7,4.9));a=fig.add_axes([.12,.14,.39,.70]);b=fig.add_axes([.67,.22,.29,.59])
panel(fig,.03,'Nuclei retained by donor','A');panel(fig,.59,'Recorded age and coverage','B')
col=d.diagnosis.map({'AD':AD,'Control':CONTROL})
for i,r in d.iterrows():
    a.hlines(i,0,r.astrocyte_nuclei,color=col.iloc[i],lw=1.25,alpha=.6)
    a.text(r.astrocyte_nuclei+12,i,str(r.astrocyte_nuclei),va='center',fontsize=7.5)
a.scatter(d.astrocyte_nuclei,np.arange(18),c=col,s=25)
a.set_yticks(range(18),d.donor_id.str.removeprefix('Sample-'));a.invert_yaxis()
a.set_ylabel('Donor ID (Sample-)');a.set_xlim(0,620);a.set_xticks([0,200,400,600]);a.set_xlabel('Astrocyte nuclei')
a.axvline(20,color=DARK,lw=.7,ls='--');a.axhline(10.5,color='#D7D5D6',lw=.6)
fig.legend(handles=handles(),loc='upper left',bbox_to_anchor=(.12,.895),frameon=False,ncol=2)
for sex,marker in [('F','o'),('M','^')]:
    z=d[d.sex==sex];b.scatter(z.recorded_age_score,z.astrocyte_nuclei,s=42,c=z.diagnosis.map({'AD':AD,'Control':CONTROL}),marker=marker,edgecolor='white',lw=.5)
b.set_xlim(72,92);b.set_xticks([75,80,85,90],['75','80','85','90+']);b.set_ylim(0,620)
b.set_xlabel('Recorded age (years)');b.set_ylabel('Astrocyte nuclei')
b.legend(handles=[Line2D([],[],marker='o',ls='',color=DARK,label='Female'),Line2D([],[],marker='^',ls='',color=DARK,label='Male')],loc='upper left',frameon=False,fontsize=7.6)
fig.text(.12,.035,'18 independent donors: 11 AD, 7 control  |  Dashed line: 20-nucleus inclusion screen',fontsize=8)
fig.text(.67,.115,'90+ is a capped recorded score.',fontsize=7.8)
save(fig,'figure_1_study_overview',[sp,cp],[data('figure_1_donors',d)],'Donor coverage and recorded-age structure')

# All tested genes remain nonsignificant; use the existing age-score estimates.
dp='results/tables/differential_expression.tsv';ap='results/sensitivity/age_score_differential_expression.tsv'
de=load(dp);age=load(ap);assert not (de.FDR<FDR).any() and not (age.FDR<FDR).any()
de['minus_log10_raw_p']=-np.log10(np.maximum(de.PValue,np.nextafter(0.,1.)))
selected=de.sort_values(['PValue','gene_id']).head(10).copy()
paired=selected[['gene_id','gene_symbol','logFC','PValue','FDR']].merge(age[['gene_id','logFC','FDR']],on='gene_id',validate='one_to_one',suffixes=('_primary','_age'))
fig=plt.figure(figsize=(8.7,4.5));a=fig.add_axes([.09,.18,.44,.64]);b=fig.add_axes([.71,.18,.25,.64])
panel(fig,.025,'Gene-level evidence','A');panel(fig,.61,'Recorded-age sensitivity','B')
a.scatter(de.logFC,de.minus_log10_raw_p,s=4,c=GRAY,alpha=.5,edgecolors='none')
a.axvline(0,color='#A8A3AA',lw=.7,ls=':');lim=max(4.,float(np.ceil(de.logFC.abs().max())))
a.set_xlim(-lim,lim);a.set_ylim(0,5.9);a.set_xlabel('log2 fold change (AD / control)');a.set_ylabel('-log10(raw P value)')
positions=[(-2.6,5.1),(.75,5.45),(3.05,4.8),(1.8,4.35),(3.15,5.55),(-2.6,3.8)]
for (_,r),pos in zip(de.sort_values(['PValue','gene_id']).head(6).iterrows(),positions):
    a.annotate(r.gene_symbol,xy=(r.logFC,r.minus_log10_raw_p),xytext=pos,ha='center',fontsize=7.5,arrowprops={'arrowstyle':'-','lw':.55,'color':'#928994','shrinkA':3,'shrinkB':3})
for i,r in paired.iterrows():b.plot([r.logFC_primary,r.logFC_age],[i,i],color='#B7AFB9',lw=1)
b.scatter(paired.logFC_primary,range(10),marker='D',s=26,c=AD,label='Primary')
b.scatter(paired.logFC_age,range(10),s=27,facecolor='white',edgecolor=AD,label='Age-score adjusted',lw=.9)
b.set_yticks(range(10),paired.gene_symbol,fontsize=7.7);b.invert_yaxis();b.axvline(0,color='#A8A3AA',ls=':',lw=.7)
b.set_xlim(-1.3,3.3);b.set_xticks([-1,0,1,2,3]);b.set_xlabel('log2 fold change')
b.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.05,1.16),fontsize=7.6)
fig.text(.09,.87,f'{len(de):,} genes; no BH FDR < {FDR:.2f}',fontsize=9,weight='bold')
fig.text(.09,.055,'Labels: six smallest primary P values. Right: ten primary-P-selected genes; no gene passes FDR in either model.',fontsize=7.8)
save(fig,'figure_2_differential_expression',[dp,ap],[data('figure_2_all_genes',de),data('figure_2_model_comparison',paired)],'Gene-level null discoveries and recorded-age effect sensitivity')

# Align effects and FDR by gene-set key; no new hypothesis tests.
ep='results/tables/enrichment_results.tsv';rp='results/sensitivity/correlation_all_sets.tsv';asp='results/sensitivity/age_score_enrichment.tsv'
en=load(ep);rho=load(rp);aset=load(asp)
top=en.sort_values(['PValue','gene_set']).head(10).copy();rhos=[0.,.01,.03,.05]
matrix=rho.pivot(index='gene_set',columns='correlation',values='FDR').reindex(index=top.gene_set,columns=rhos)
age_q=aset.set_index('gene_set').loc[top.gene_set,'FDR'].to_numpy()
assert np.allclose(matrix[.01].values,top.FDR,rtol=1e-12,atol=1e-14)
display={'TNFA_SIGNALING_VIA_NFKB':'TNFα signaling via NF-κB','INFLAMMATORY_RESPONSE':'Inflammatory response',
 'UNFOLDED_PROTEIN_RESPONSE':'Unfolded protein response','BILE_ACID_METABOLISM':'Bile acid metabolism',
 'IL6_JAK_STAT3_SIGNALING':'IL6-JAK-STAT3 signaling','MTORC1_SIGNALING':'mTORC1 signaling',
 'PI3K_AKT_MTOR_SIGNALING':'PI3K-AKT-mTOR signaling','APOPTOSIS':'Apoptosis',
 'ESTROGEN_RESPONSE_EARLY':'Estrogen response (early)',
 'EPITHELIAL_MESENCHYMAL_TRANSITION':'Epithelial-mesenchymal transition'}
names=[display.get(g.removeprefix('HALLMARK_'),g.removeprefix('HALLMARK_').replace('_',' ').lower()) for g in top.gene_set]
fig=plt.figure(figsize=(8.7,5.6));a=fig.add_axes([.32,.32,.22,.52]);b=fig.add_axes([.61,.32,.25,.52]);c=fig.add_axes([.91,.32,.065,.52])
panel(fig,.025,'Rank effects','A');panel(fig,.60,'Correlation scenarios','B');panel(fig,.88,'Age model','C')
for i,r in enumerate(top.itertuples()):
    color=AD if r.rank_biserial>0 else CONTROL
    a.hlines(i,0,r.rank_biserial,lw=1.4,color=color)
    a.scatter([r.rank_biserial],[i],s=43,facecolors=color if r.FDR<FDR else 'white',edgecolors=color,lw=1.1,zorder=3)
a.axvline(0,color='#8B858E',lw=.65);a.set_ylim(9.5,-.5);a.set_yticks(range(10),names,fontsize=10.5)
a.set_xlim(-.34,.38);a.set_xticks([-.3,0,.3]);a.tick_params(labelsize=10.5)
a.set_xlabel('Rank-biserial effect\nControl ← 0 → AD',fontsize=10.5)
cmap=LinearSegmentedColormap.from_list('evidence',['#FAF7FA','#D3B8CD',AD]);vmax=float(np.ceil(-np.log10(min(matrix.to_numpy().min(),age_q.min()))));norm=Normalize(0,vmax)
for axis,values in [(b,matrix.to_numpy()),(c,age_q[:,None])]:
    axis.pcolormesh(np.arange(values.shape[1]+1)-.5,np.arange(11)-.5,-np.log10(values),cmap=cmap,norm=norm,edgecolors='white',lw=.8)
    axis.set_ylim(9.5,-.5);axis.set_yticks([]);axis.spines[['top','right','bottom','left']].set_visible(False)
    for i in range(10):
        for j in range(values.shape[1]):
            q=values[i,j];label=f'{q:.3f}' if q>=.001 else f'{q:.1e}'
            axis.text(j,i,label,ha='center',va='center',fontsize=10.5,color='white' if norm(-np.log10(q))>.65 else DARK)
            if q<FDR:axis.add_patch(Rectangle((j-.46,i-.46),.92,.92,fill=False,ec=DARK,lw=.85))
b.set_xticks(range(4),['0','0.01','0.03','0.05']);b.xaxis.tick_top();b.tick_params(length=0,labelsize=10.5)
c.set_xticks([0],['ρ = 0.01']);c.xaxis.tick_top();c.tick_params(length=0,labelsize=10.5)
counts=[int(((rho.correlation==r)&(rho.FDR<FDR)).sum()) for r in rhos]
fig.text(.61,.25,'Discoveries / 49 sets:',fontsize=10.5)
for j,count in enumerate(counts):fig.text(.61+.25*(j+.5)/4,.205,str(count),ha='center',fontsize=11,weight='bold')
bar=fig.add_axes([.72,.12,.22,.018]);legend_cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=bar,orientation='horizontal',ticks=np.arange(vmax+1))
legend_cb.solids.set_edgecolor('face')
bar.set_title('-log10(FDR)',fontsize=10.5,pad=5);bar.tick_params(labelsize=10.5,length=2)
fig.text(.025,.055,'Top 10 primary-P-ranked sets. Filled points / outlines: FDR < 0.05.',fontsize=10.5)
fig.text(.025,.018,'B: fixed ranks. C: recorded-age model with a refitted ranking.',fontsize=10.5)
plot=top.copy()
for r in rhos:plot[f'FDR_rho_{r:g}']=matrix[r].values
plot['FDR_age_score']=age_q
save(fig,'figure_3_gene_set_enrichment',[ep,rp,asp],[data('figure_3_program_evidence',plot),data('figure_3_full_family_counts',pd.DataFrame({'correlation':rhos,'significant_sets':counts,'tested_sets':len(en)}))],'Hallmark rank effects and separate correlation and age-model evidence')

# Ranges and donor-resolved changes are descriptive sensitivities, not CIs.
lp='results/tables/lodo_selected_results.tsv';su='results/tables/lodo_selected_summary.tsv'
lodo=load(lp);summary=load(su).sort_values('selection_rank').reset_index(drop=True)
delta=lodo.copy();delta['delta_logFC']=delta.logFC-delta.full_logFC;order=d.donor_id.tolist()
heat=delta.pivot(index='gene_id',columns='omitted_donor',values='delta_logFC').reindex(index=summary.gene_id,columns=order)
assert not heat.isna().any().any() and heat.shape==(10,18)
fig=plt.figure(figsize=(8.7,4.4));a=fig.add_axes([.13,.23,.34,.59]);b=fig.add_axes([.57,.23,.39,.59])
panel(fig,.025,'Full effects and omission ranges','A');panel(fig,.55,'Change after omitting each donor','B')
for i,r in enumerate(summary.itertuples()):
    z=lodo[lodo.gene_id==r.gene_id].set_index('omitted_donor').loc[order]
    a.hlines(i,r.min_logFC,r.max_logFC,color='#AAA0AE',lw=1.3)
    a.scatter(z.logFC,i+np.linspace(-.17,.17,18),s=12,c=z.omitted_diagnosis.map({'AD':AD,'Control':CONTROL}),edgecolors='none',alpha=.8)
    a.scatter([r.full_logFC],[i],marker='D',s=25,c=DARK,edgecolor='white',lw=.4,zorder=4)
a.set_yticks(range(10),summary.gene_symbol,fontsize=8);a.set_ylim(9.5,-.5);a.set_xlim(-1.3,3.2)
a.set_xticks([-1,0,1,2,3]);a.axvline(0,color='#AAA0AE',lw=.7,ls=':');a.set_xlabel('log2 fold change (AD / control)')
cmap_delta=LinearSegmentedColormap.from_list('delta',[CONTROL,'#FCFAFC',AD]);bound=float(np.ceil(abs(heat.to_numpy()).max()*10)/10)
mesh=b.pcolormesh(np.arange(19)-.5,np.arange(11)-.5,heat.to_numpy(),cmap=cmap_delta,norm=TwoSlopeNorm(0,-bound,bound),edgecolors='white',lw=.3)
b.set_ylim(9.5,-.5);b.set_yticks([]);b.set_xticks(range(18),[x.removeprefix('Sample-') for x in order],rotation=90,fontsize=7)
b.tick_params(length=0);b.set_xlabel('Omitted donor ID (Sample-)',fontsize=8);b.spines[['top','right','bottom','left']].set_visible(False)
b.axvline(10.5,color=DARK,lw=.8)
for start,end,label,color in [(0,11,'AD omitted',AD),(11,18,'Control omitted',CONTROL)]:
    b.plot([start/18,end/18],[1.035,1.035],transform=b.transAxes,color=color,lw=3,clip_on=False)
    b.text((start+end)/36,1.075,label,transform=b.transAxes,ha='center',fontsize=7.5,color=color)
bar=fig.add_axes([.64,.075,.25,.025]);legend_cb=fig.colorbar(mesh,cax=bar,orientation='horizontal',ticks=[-bound,0,bound])
legend_cb.solids.set_edgecolor('face')
bar.set_xlabel('Omission log2FC − full-data log2FC',fontsize=7.7,labelpad=2);bar.tick_params(labelsize=7,length=2)
fig.legend(handles=[Line2D([],[],marker='D',ls='',color=DARK,label='Full data')]+handles(),loc='upper left',bbox_to_anchor=(.13,.906),frameon=False,ncol=3,fontsize=7.7)
fig.text(.025,.025,'Ten genes selected by full-data P; 18 omissions each. Spans are LODO ranges, not confidence intervals.',fontsize=7.8)
save(fig,'figure_4_donor_sensitivity',[lp,su,cp,sp],[data('figure_4_omission_effects',delta),data('figure_4_ranges',summary)],'Selected gene ranges and donor-specific effect changes')

# One common donor embedding, shown with separately retained source covariates.
mp='results/targeted_validation/donor_mds.tsv'
mjp='results/targeted_validation/mds_summary.json'
mds=load(mp)
sources[mjp]=digest(ROOT/mjp)
mds_info=json.loads((ROOT/mjp).read_text())
assert mds.donor_id.is_unique and set(mds.donor_id)==set(d.donor_id)
assert np.isfinite(mds[['MDS1','MDS2']].to_numpy()).all()
variance=np.asarray(mds_info['variance_percent'],dtype=float)
fig=plt.figure(figsize=(8.7,6.9))
axs=[fig.add_axes([x,y,.245,.285]) for y in [.59,.15] for x in [.08,.405,.73]]
xspan=np.ptp(mds.MDS1);yspan=np.ptp(mds.MDS2)
for axis in axs[:5]:
    axis.set_xlim(mds.MDS1.min()-.12*xspan,mds.MDS1.max()+.12*xspan)
    axis.set_ylim(mds.MDS2.min()-.12*yspan,mds.MDS2.max()+.12*yspan)
    axis.set_xlabel(f'Dimension 1 ({variance[0]:.1f}%)',fontsize=10)
    axis.set_ylabel(f'Dimension 2 ({variance[1]:.1f}%)',fontsize=10)
    axis.tick_params(labelsize=9)
for axis,title,letter in zip(axs,['Diagnosis','Recorded age','Recorded sex','GEO Batch','Supplementary Library.Batch','Reading this diagnostic'],'ABCDEF'):
    axis.set_title(f'{letter}  {title}',loc='left',fontsize=10.5,pad=31)
axs[0].scatter(mds.MDS1,mds.MDS2,c=mds.diagnosis.map({'AD':AD,'Control':CONTROL}),s=38,edgecolors='white',linewidths=.5)
axs[0].legend(handles=handles(),loc='lower left',bbox_to_anchor=(-.03,1.005),ncol=2,frameon=False,fontsize=9,handletextpad=.25,columnspacing=.8)
age_norm=Normalize(mds.recorded_age_score.min(),90)
age_mesh=axs[1].scatter(mds.MDS1,mds.MDS2,c=mds.recorded_age_score,cmap='viridis',norm=age_norm,s=38,edgecolors='white',linewidths=.5)
age_bar=fig.add_axes([.405,.91,.245,.013]);legend_cb=fig.colorbar(age_mesh,cax=age_bar,orientation='horizontal',ticks=[75,80,85,90])
legend_cb.solids.set_edgecolor('face')
age_bar.set_xticklabels(['75','80','85','90+']);age_bar.tick_params(labelsize=8.5,length=2,pad=1)
sex_colors={'F':AD,'M':CONTROL}
for sex,marker,label in [('F','o','Female'),('M','^','Male')]:
    part=mds[mds.sex==sex]
    axs[2].scatter(part.MDS1,part.MDS2,c=sex_colors[sex],marker=marker,s=40,label=label,edgecolors='white',linewidths=.5)
axs[2].legend(loc='lower left',bbox_to_anchor=(-.03,1.005),ncol=2,frameon=False,fontsize=9,handletextpad=.25,columnspacing=.8)
batch_colors=['#2F6970','#AD7436','#835579']
for axis,field in [(axs[3],'GEO_batch'),(axs[4],'supplementary_library_batch')]:
    levels=sorted(mds[field].dropna().unique())
    for level,color in zip(levels,batch_colors):
        part=mds[mds[field]==level]
        axis.scatter(part.MDS1,part.MDS2,c=color,s=38,label=f'{int(level)}',edgecolors='white',linewidths=.5)
    axis.legend(loc='lower left',bbox_to_anchor=(-.03,1.005),ncol=len(levels),frameon=False,fontsize=9,handletextpad=.15,columnspacing=.8)
axs[5].axis('off')
axs[5].text(0,1,'18 independent donors\n\nSame coordinates in A-E\n500 common high-variance\ngenes; TMM logCPM\n\n90+ is top-coded.\nBatch fields retain their\nseparate source identities.\n\nNo donor is excluded.',va='top',fontsize=10.5,linespacing=1.35,transform=axs[5].transAxes)
fig.text(.08,.045,'Distances summarize expression structure; visual grouping cannot establish or exclude confounding.',fontsize=10.3)
save(fig,'figure_s1_donor_structure',[mp,mjp],[data('figure_s1_donor_structure',mds)],'Common donor MDS coordinates annotated by diagnosis and recorded covariates')

# Fixed primary discoveries: full-data effects, every omission and per-fit FDR.
gsp='results/targeted_validation/gene_set_lodo.tsv'
gss='results/targeted_validation/gene_set_lodo_summary.tsv'
gl=load(gsp);gs=load(gss)
five=en[en.FDR<FDR].sort_values(['PValue','gene_set']).copy()
gsummary=gs.set_index('gene_set').loc[five.gene_set].reset_index()
assert len(five)==5 and gl.shape[0]==49*18
gheat=gl.pivot(index='gene_set',columns='omitted_donor',values='FDR').loc[five.gene_set,order]
assert not gheat.isna().any().any()
gnames=[display.get(g.removeprefix('HALLMARK_'),g.removeprefix('HALLMARK_').replace('_',' ').lower()) for g in five.gene_set]
fig=plt.figure(figsize=(8.7,4.8));a=fig.add_axes([.27,.34,.26,.45]);b=fig.add_axes([.60,.34,.37,.45])
panel(fig,.025,'Effects across donor omissions','A');panel(fig,.59,'BH FDR after each omission','B')
for i,r in enumerate(gsummary.itertuples()):
    part=gl[gl.gene_set==r.gene_set].set_index('omitted_donor').loc[order]
    a.hlines(i,r.min_effect,r.max_effect,lw=1.3,color='#AAA0AE')
    a.scatter(part.rank_biserial,i+np.linspace(-.17,.17,18),s=17,c=part.omitted_diagnosis.map({'AD':AD,'Control':CONTROL}),edgecolors='none',alpha=.9)
    a.scatter([r.primary_rank_biserial],[i],s=34,marker='D',c=DARK,edgecolor='white',lw=.4,zorder=4)
a.set_yticks(range(5),gnames,fontsize=10.5);a.set_ylim(4.5,-.5);a.axvline(0,color=GRAY,lw=.7,ls=':')
bound_effect=float(np.ceil(max(abs(gsummary.min_effect).max(),abs(gsummary.max_effect).max())*10)/10)
a.set_xlim(-bound_effect-.02,bound_effect+.02);a.tick_params(labelsize=10)
a.set_xlabel('Rank-biserial effect\nControl ← 0 → AD',fontsize=10.5)
q=gheat.to_numpy();qlog=-np.log10(np.maximum(q,np.nextafter(0.,1.)))
qbound=max(2.,float(np.ceil(qlog.max())))
qmesh=b.pcolormesh(np.arange(19)-.5,np.arange(6)-.5,qlog,cmap=cmap,norm=Normalize(0,qbound),edgecolors='white',lw=.4)
for i in range(5):
    for j in range(18):
        if q[i,j]<FDR:b.plot(j,i,marker='o',markersize=3.5,mec=DARK,mfc='white',mew=.45,lw=0)
b.set_ylim(4.5,-.5);b.set_yticks([]);b.set_xticks(range(18),[x.removeprefix('Sample-') for x in order],rotation=90,fontsize=9)
b.tick_params(length=0);b.set_xlabel('Omitted donor ID (Sample-)',fontsize=10.5)
b.spines[['top','right','bottom','left']].set_visible(False);b.axvline(10.5,color=DARK,lw=.8)
for start,end,label,color in [(0,11,'AD omitted',AD),(11,18,'Control omitted',CONTROL)]:
    b.plot([start/18,end/18],[1.035,1.035],transform=b.transAxes,color=color,lw=3,clip_on=False)
    b.text((start+end)/36,1.11,label,transform=b.transAxes,ha='center',fontsize=10,color=color)
fig.legend(handles=[Line2D([],[],marker='D',ls='',color=DARK,label='Full data')]+handles(),loc='upper left',bbox_to_anchor=(.27,.91),frameon=False,ncol=3,fontsize=9)
bar=fig.add_axes([.66,.14,.25,.025]);legend_cb=fig.colorbar(qmesh,cax=bar,orientation='horizontal',ticks=np.arange(qbound+1))
legend_cb.solids.set_edgecolor('face')
bar.set_title('-log10(FDR); white dots: FDR < 0.05',fontsize=9.5,pad=5);bar.tick_params(labelsize=9,length=2)
fig.text(.025,.065,'Five sets selected by full-data FDR. BH covers all 49 sets separately within each omission.',fontsize=10)
fig.text(.025,.02,'Fixed correlation 0.01. Effect spans are influence ranges; omissions are not independent replications.',fontsize=10)
save(fig,'figure_s2_program_sensitivity',[gsp,gss,ep,sp,cp],[data('figure_s2_program_omissions',gl[gl.gene_set.isin(five.gene_set)]),data('figure_s2_program_ranges',gsummary)],'Primary-five gene-set effects and complete donor-omission FDR pattern')
print('PASS: four main and two supplementary figures regenerated; PDF/SVG/400-dpi PNG and plotting tables exported.')

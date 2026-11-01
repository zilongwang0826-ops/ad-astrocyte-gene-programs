# Astrocyte expression in Alzheimer's disease: a donor-level reanalysis

**Zilong Wang** | Université Paris-Saclay, Master 2 Bioinformatics & Biostatistics

Research project summary | Project period: August 2026

## Research question and approach

Using public GSE174367 data, I combined counts from [[astrocytes]] astrocyte nuclei into profiles for [[ad_donors]] AD and [[control_donors]] control donors. I compared the groups with edgeR, then checked how age adjustment, gene correlation assumptions and donor omissions affected the results. The donor, rather than each nucleus, was the biological replicate.

## Findings

No individual gene passed correction among [[tested_genes]] tested genes. Inflammatory and protein-stress gene sets ranked toward AD. At the default correlation of [[set_correlation]], [[significant_sets]] of [[tested_sets]] Hallmark sets met BH FDR < [[fdr_threshold]]. Age adjustment retained [[age_retained_original]] of the original five; none remained significant at correlations [[rho_mid]] or [[rho_high]]. All five kept their direction across [[lodo_iterations]] donor omissions, but each met FDR < [[fdr_threshold]] in only [[program_min_fdr_retained]]-[[program_max_fdr_retained]] omissions.

![Figure 3](../../results/figures/figure_3_gene_set_enrichment.png)

**Figure | Program ranking and statistical support.** A, Primary rank-biserial effects; filled points meet default FDR < [[fdr_threshold]]. B, FDR across four assumed correlations, holding ranks fixed. C, The recorded-age model. Outlined cells meet the same FDR threshold; color shows evidence, not pathway activity.

**Interpretation.** The gene-set rankings describe this cohort, with statistical support that depends on the model. Age and sex imbalance, capped ages and conflicting batch labels limit the comparison; donor-level expression also mixes cell-state proportions with expression within states. The overlapping donor omissions do not provide independent validation.

**Data and code.** Morabito et al., Nature Genetics (2021), doi:10.1038/s41588-021-00894-z; GEO GSE174367. Code, complete tables and supplementary methods accompany the report.

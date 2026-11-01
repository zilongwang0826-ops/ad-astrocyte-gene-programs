# Supplementary methods and statistical tables

Project period: August 2026.

## Data provenance and primary specification

The reanalysis uses the two processed GSE174367 snRNA-seq files listed in `data/input_manifest.json`, acquired on 29 October 2026 according to the retained acquisition record. The original reference label is `GRCh38.p12.premrna2`; original methods identify Cell Ranger 3.0.2 and pre-mRNA GRCh38.p12, without an exact GTF release. Counts are integers after the authors' aggregation and read-depth subsampling.

The HDF5 contains [[input_genes]] features and [[matrix_nuclei]] barcodes. All [[metadata_nuclei]] final metadata nuclei match; [[matrix_only]] HDF5-only barcodes are listed explicitly in the audit and excluded because final donor/cell-type annotations are absent. Sorting stored sparse row indices changes representation without changing counts. Aggregating annotated ASC nuclei conserves [[astrocyte_umis]] UMIs. Gene IDs remain distinct in differential expression; enrichment resolves duplicate symbols by largest primary logCPM, then gene ID.

The primary model is ~ diagnosis, with Control as reference. The coverage rules were fixed before differential-expression results were inspected; the local pre-analysis document is not a registered protocol. R [[r_version]], edgeR [[edger_version]] and limma [[limma_version]] were used. `filterByExpr` used min.count=10, min.total.count=15, large.n=10 and min.prop=0.7, with the primary design. Filtered library sizes were recomputed. TMM used logratioTrim=0.3, sumTrim=0.05, doWeighting=TRUE and Acutoff=-1e10. The QL fit used robust=TRUE, legacy=FALSE, abundance.trend=TRUE, winsor.tail.p=(0.05,0.1) and prior.count=0.125; NB dispersion was estimated internally from high-abundance genes. The tested coefficient was diagnosisAD. Complete statistics are retained for all [[tested_genes]] tested genes; BH correction used that full family.

The original Hallmark [[hallmark_version]] collection contains [[source_sets]] sets. Mapping retains [[ranked_symbols]] unique symbols and [[tested_sets]] sets with [[min_set_size]]–[[max_set_size]] genes. Pancreas beta cells is excluded because only [[excluded_set_size]] genes map. Ranking is sign(logFC) × sqrt(QL F). cameraPR uses ranks, directional two-sided competitive testing, and BH across the eligible sets. Its default within-set correlation is [[set_correlation]]; this is a software default rather than an estimate from these donors. Rank-biserial effect uses average ranks for ties and equals Pr(set gene ranks above a background gene) minus Pr(it ranks below).

## Sensitivity specification

Correlation scenarios and one recorded-age-score model follow the plan in `results/sensitivity/sensitivity_plan.md`, written before those results were inspected. At correlations [[rho_zero]], [[set_correlation]], [[rho_mid]] and [[rho_high]], the original ranking, symbol mapping and set membership are frozen; BH is applied separately across the same [[tested_sets]] sets in each scenario. All scenarios are reported as fixed assumptions rather than estimates of within-set correlation.

The supplementary model is ~ diagnosis + recorded_age_centered. A value of 90 denotes the source's 90-or-older category; centering the deposited values does not recover unknown ages. The linear capped-score assumption has limited overlap: the sole exactly shared recorded age is the upper category. The design has rank [[age_design_rank]], [[age_residual_df]] residual degrees of freedom and predictor VIF [[age_vif]]. Detailed conditioning and range checks are in `results/sensitivity/design_preflight.json`. Age has no missing values, and no donor is excluded. Sex and batch are tabulated but no further model is fitted. The original gene-filtering decisions, TMM factors, canonical symbols and mapped sets are reused; QL inference is refitted with the additional predictor. Adjusted ranked testing uses the default correlation [[set_correlation]].

The LODO analysis keeps the full-data filtered universe and refits TMM plus dispersion after each donor omission. Its [[displayed_genes]]-gene display is fixed by the smallest full-data raw P values. The maximum displayed change is [[max_lodo_change]] log2FC; EMX2 ranges from [[emx2_lodo_min]] to [[emx2_lodo_max]]. These are dependent influence summaries of the fixed display panel. Complete LODO tables remain in `results/tables/`.

## Supplementary Table S1. Donor metadata

Values are taken from deposited metadata and the previously verified source supplement. The downloaded metadata and source table do not explicitly state the PMI unit, so its values are reproduced without assigning an inferred unit. NA denotes missing PMI; 90+ denotes top coding, not an exact age. Original supplementary Library.Batch and GEO Batch are retained separately because their assay/stage relationship is unresolved. The source supplement lists a different earlier nucleus total; all nucleus counts here use final deposited annotations.

### Group-level characteristics

[[COHORT_TABLE]]

### Individual donors

[[DONOR_TABLE]]

Group distributions and batch cross-tabulations are also provided in `results/sensitivity/clinical_covariate_summary.tsv`, `sex_by_diagnosis.tsv` and `batch_by_diagnosis.tsv`.

## Supplementary Table S2. Correlation sensitivity

Counts across all eligible sets:

[[CORRELATION_COUNTS]]

BH FDR for the five sets identified in the original analysis:

[[CORRELATION_TABLE]]

Directions are inherited from the fixed ranking. Complete statistics for every set and scenario are in `results/sensitivity/correlation_all_sets.tsv`; no nonsignificant set was removed from that output.

## Supplementary Table S3. Recorded-age-score sensitivity

EMX2 diagnosis-effect comparison:

[[EMX2_TABLE]]

Original five Hallmark sets at the default correlation:

[[AGE_TABLE]]

The adjusted model has [[age_gene_discoveries]] gene discoveries and [[age_all_set_discoveries]] significant sets overall. [[age_retained_original]] of the original five sets remain below FDR [[fdr_threshold]]. PI3K-AKT-mTOR signaling and protein secretion are the additional sets crossing the threshold; complete adjusted results are in `results/sensitivity/age_score_differential_expression.tsv` and `age_score_enrichment.tsv`.

## Targeted donor and gene-set validation

The added validation starts from the included pseudobulk counts and preserves the historical primary and sensitivity outputs. The full-data edgeR fit is checked against the saved primary table before donor omissions. This reconstructs inference from the supplied donor counts; it does not independently reconstruct the original raw-data aggregation or resolve the original reference annotation.

The donor embedding uses TMM logCPM with prior count 2 and the same [[mds_top_genes]] high-variance genes for every distance (`gene.selection="common"`, `method="logFC"`). Within the primary retained universe, variable-gene selection does not use diagnosis labels. The first two dimensions explain [[mds1_variance]]% and [[mds2_variance]]% of common-gene variation. Diagnosis, age, sex and two source batch fields annotate one unchanged embedding. No donor is removed based on the plot. Age 90+ remains top-coded, and the two batch fields are not substituted for one another.

![Figure S1](../../results/figures/figure_s1_donor_structure.png)

**Figure S1 | Donor expression structure and recorded covariates.** Each point is one of [[donors]] donors. A-E show identical MDS coordinates colored by diagnosis, capped recorded age, sex, GEO Batch and original supplementary Library.Batch, respectively. F summarizes the diagnostic scope. These descriptive views do not test a causal diagnosis effect and cannot demonstrate absence of confounding. Complete coordinates and covariates are in `results/targeted_validation/donor_mds.tsv`.

The gene-set LODO analysis refits TMM and the edgeR QL model after each single-donor omission while retaining the full-data filtered gene universe and canonical gene-id-to-symbol mapping. Signed square-root QL F statistics define the new ranking. Frozen memberships identify the same [[tested_sets]] eligible sets in every iteration; cameraPR uses ranks, directional testing and fixed correlation [[set_correlation]]. BH correction is performed separately across all [[tested_sets]] sets within each omission, not across only the displayed discoveries. All 882 set-by-omission records are retained. No across-omission probability or independent replication claim is assigned to threshold-retention counts.

![Figure S2](../../results/figures/figure_s2_program_sensitivity.png)

**Figure S2 | Influence on the five original Hallmark discoveries.** A, Diamonds show full-data rank-biserial effects, points show every donor omission, and spans give the observed minimum and maximum. B, Color represents -log10(FDR), with white dots marking FDR < [[fdr_threshold]] after correction over all [[tested_sets]] sets in that omission. The fixed display panel was selected by full-data FDR; no later nonsignificant result is removed. These are influence summaries under one correlation assumption, not confidence intervals or external validation.

## Supplementary Table S4. Gene-set donor influence

[[PROGRAM_LODO_TABLE]]

Direction and threshold counts use [[lodo_iterations]] dependent omissions. Complete results for every set are in `results/targeted_validation/gene_set_lodo.tsv`, with summaries in `gene_set_lodo_summary.tsv`. The original-five maximum absolute effect change is [[program_max_effect_change]]. The analysis plan, runtime records and baseline comparisons accompany the targeted validation outputs.

## Reproduction and interpretation

`scripts/07_sensitivity.R` produces the sensitivity analyses, and `tests/test_sensitivity.R` checks synthetic cases and empirical output consistency. `results/audit/sensitivity_validation.json` is the machine-readable result summary. Original analysis parameters and outputs remain unchanged. Computation was sequential with numerical thread limits of two. The sensitivity checks quantify dependence on the specified assumptions. Batch-source definitions and the observed donor metadata are retained explicitly. Reproduction commands and validation scope are described in `VALIDATION.md`.

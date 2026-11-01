# Pre-analysis specification

Recorded before expression files or differential-expression results were available.

## Scientific question and unit
Which astrocyte transcriptional programs differ between Alzheimer's disease (AD) and control donors, and how stable are the estimated differences to removal of individual donors? The independent biological unit is the donor. Nuclei and technical libraries are not independent replicates.

## Fixed scope and feasibility gate
- Dataset: GSE174367, processed snRNA-seq counts and author cell metadata only.
- Population: the author-annotated major astrocyte population; no reclustering or new annotation.
- Comparison: AD versus control; positive log2 fold change means higher in AD.
- Screen: at least 20 annotated astrocyte nuclei per donor; retain at least five independent donors per diagnosis.
- Required: reliable donor/diagnosis mapping, raw nonnegative integer counts, documented barcode reconciliation and no unresolved fatal design confounding.
- A failed gate produces BLOCKED_DATA_AUDIT. Dataset substitution and fabricated results are prohibited.

## Primary analysis
Sum counts within donor. Use edgeR filtering and TMM normalization, a negative-binomial quasi-likelihood model with diagnosis as the primary predictor, and BH FDR < 0.05. Inspect available age, sex, tissue quality and technical variables before finalizing the model; justify any covariate deviation before inference. Report the complete tested-gene table. Absolute log2FC >= 0.25 may describe magnitude but does not alter the primary significance criterion.

## Functional analysis and sensitivity
Use one documented ranked gene-set test, one versioned collection, unique mapped identifiers, a signed primary-model statistic and BH correction across all tested sets. Final choice depends on verified resource availability and identifier coverage and must be recorded before viewing enrichment results. LODO refits the same model, normalization and dispersion estimation after removing each donor. The display panel is fixed from the full fit: at most ten genes with smallest raw P values, ties resolved by stable gene identifier; do not reselect per iteration. Report sign consistency, effect ranges and the most influential omission. LODO is internal sensitivity analysis, not external replication.

## Acceptance criteria
Real files pass the feasibility gate; donor-level aggregation and comparison orientation pass synthetic and empirical checks; all analysis stages execute; every reported numerical claim links to machine-readable results; exactly four main figures are rendered and visually inspected; report and one-page summary PDFs are searchable and every page is visually inspected; dependency and input provenance are recorded. No external publishing is authorized.

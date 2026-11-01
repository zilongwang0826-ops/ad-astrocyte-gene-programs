# Targeted validation plan, 31 October 2026

This local plan implements the approved Phase 1 scope. It was recorded before inspecting the new results; it is not a preregistration. The earlier `results/sensitivity/sensitivity_plan.md` describes the preceding age/correlation revision and remains a historical record.

## Question and unit

Do donor expression structure and individual donor omissions alter interpretation of the observed AD/control astrocyte program rankings? The independent unit is the donor (11 AD, 7 control). Both checks are exploratory within-cohort diagnostics, not external replication or causal identification.

## Inputs and baseline gate

Use the distributed pseudobulk counts, sample/annotation tables, primary filtering and gene results, canonical symbol map, frozen mapped Hallmark memberships, and donor covariates. Verify donor and identifier alignment and recorded package versions. Refit the unchanged primary model with the original expression filter, TMM, edgeR QL settings and contrast. Compare filtering, normalization, logFC, logCPM, F, PValue and FDR to the supplied outputs before omissions. Numeric tolerance is abs(error) <= 1e-10 + 1e-8 * abs(reference); identifiers, directions of primary gene-set effects and discovery decisions must agree exactly. Stop on failure rather than replace the baseline.

This reconstructs inference from existing donor counts. It does not rerun raw barcode alignment, aggregation, original annotation or the authors' study. A historical dataset audit is a prerequisite, not evidence that raw data were rechecked in this run.

## Donor structure

Use edgeR logCPM (TMM, prior.count=2) and limma/edgeR MDS with top=500, gene.selection="common", method="logFC". Selection uses expression variance, not diagnostic separation or P values. Save the selected genes and coordinates; use identical coordinates for diagnosis, capped age, sex and the two separately labelled batch sources. Record missingness and individual diagnosis-plus-covariate design rank; do not fit additional adjusted expression models or merge inconsistent batch fields. Mark Age=90 as 90+. No donor is excluded using the display. A two-dimensional view neither captures all variation nor proves absence of confounding.

## Program influence

Omit each donor once. Freeze the full-data tested-gene universe, canonical gene IDs per symbol and the 49 eligible Hallmark memberships. Refit TMM and NB/QL dispersion for each remaining cohort using the unchanged diagnosis model. Rank by sign(logFC)*sqrt(F), and run cameraPR at the original fixed correlation 0.01 with use.ranks=TRUE and directional=TRUE. Apply BH across all 49 sets separately in each omission. Compare the recomputed gene logFC/P/FDR with the already supplied gene-level LODO outputs.

Store every gene ranking and all 882 set results. Summarize effect direction and range, rank and FDR range, maximum effect change and its omitted donor. The original five full-data discoveries are a fixed display panel; all 49 sets remain available. Counts of threshold crossings are descriptive, not binomial success rates, confidence intervals, independent replications or new discovery claims. No joint omission-by-correlation grid or additional adjustment models will be added.

## Execution and acceptance

Run small synthetic tests first, then one baseline fit and 18 serial omissions with two numerical threads. Each omission must pass identifiers, finite-value, BH and saved-gene agreement checks before continuing. Record versions, commands, input/code SHA256, elapsed time and numerical agreement. Preserve the baseline and all scientific source outputs; write additions under results/targeted_validation. Verify generated tables and plots independently and inspect every final PDF page. Acceptance depends on correct execution and transparent interpretation, not stronger significance.

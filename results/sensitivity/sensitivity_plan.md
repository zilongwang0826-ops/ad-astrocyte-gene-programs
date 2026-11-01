# Limited sensitivity analysis plan

The revision adds two checks while preserving the existing primary analysis and its complete scientific outputs. This plan was recorded before the new sensitivity results were inspected. It is a local analysis record, not a preregistration.

## Correlation scenarios

Use the original stored signed quasi-likelihood ranking and the original mapped Hallmark gene sets. Confirm that correlation 0.01 reproduces all original P values, FDR values, directions and mapped set sizes. Then evaluate the fixed scenarios 0, 0.01, 0.03 and 0.05 using limma cameraPR with use.ranks=TRUE and directional=TRUE. Use BH correction across the same 49 eligible sets separately within each scenario. No scenario is an estimate of true biological correlation. Display all scenarios and the five sets selected by the original analysis, regardless of later significance.

The installed cameraPR.default method has inter.gene.cor=0.01 as its function default. The rank version calls rankSumTestWithCorrelation with infinite degrees of freedom and reports a two-sided competitive rank test. Effect directions are a property of the preserved signed ranking; changing the variance adjustment must not be described as changing pathway activity.

## Covariate sensitivity specification

The sole alternative model is ~ diagnosis + recorded_age_centered, selected after the metadata design check and before fitting new gene-level results. The deposited numeric Age value is centered without imputing values above the source's top-coded upper category. A recorded value of 90 represents the source's 90-or-older category. The covariate is a capped observed-age score, not exact chronological age for those donors. A linear association with this recorded score is a working sensitivity assumption; the model cannot resolve residual within-category age variation, sex or technical batch confounding.

Before fitting, assess missingness, group overlap, rank, residual degrees of freedom, diagnosis-age association and design conditioning. Preserve the primary tested-gene universe, normalization factors and canonical gene-symbol mapping. Refit the same edgeR v4 QL model with only this additional predictor; compare its diagnosis coefficient with the original unadjusted coefficient. Repeat the original Hallmark test at correlation 0.01 on the adjusted signed statistics, using the same mapped sets and BH family. This adjusted fit is a sensitivity analysis and does not replace the primary model.

No power analysis, estimated gene-correlation procedure, additional LODO analysis, new dataset or new cell type is included.

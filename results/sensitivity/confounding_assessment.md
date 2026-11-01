# Confounding assessment

Donor identifiers and clinical covariates were checked directly against deposited cell metadata. No donor was excluded. A recorded age of 90 represents the source upper category (90 or older); exact ages above that threshold are unavailable.

The sole supplementary model, ~ diagnosis + recorded_age_centered, has rank 3 and 15 residual degrees of freedom. Age is centered at 85.611111 recorded years. The diagnosis-age correlation is 0.493753, the predictor VIF 1.322388 and the centered design condition number 13.219012.

Recorded age ranges overlap from 80 to 90, but the only exactly shared recorded value is 90. That upper category contains 5 AD and 2 control donors. The fit therefore checks a linear capped-score specification; it cannot identify unobserved within-category age differences or establish adequate control of age confounding.

Sex and both source batch fields were tabulated but not included in additional fits. Female/male counts are 6/5 in AD and 2/5 in controls. GEO batch has both diagnoses in every level; the alternative Library.Batch assignment differs and its assay attribution remains unresolved. Neither table is silently substituted for the other. PMI is missing for Sample-100; age, sex, RIN and GEO batch are complete. Full donor values and distribution summaries accompany this assessment.

The primary diagnosis-only analysis, gene-filtering decision, TMM factors and original symbol/set mapping are preserved. The adjusted coefficient and enrichment ranking are a supplementary sensitivity result, not a replacement primary finding.

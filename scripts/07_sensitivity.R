#!/usr/bin/env Rscript
# Two limited sensitivity checks; original scientific outputs are never replaced.
source("scripts/r_common.R")
correlation_scenarios <- c(0, 0.01, 0.03, 0.05)

make_age_design <- function(samples) {
  assert(all(c("donor_id", "diagnosis", "Age") %in% names(samples)), "Age-design fields missing.")
  assert(!anyDuplicated(samples$donor_id), "Duplicate donors in age design.")
  assert(all(is.finite(samples$Age)) && all(samples$Age <= 90),
         "Recorded capped age must be finite and at most 90 for this deposited cohort.")
  samples$diagnosis <- factor(samples$diagnosis, levels = c("Control", "AD"))
  assert(!anyNA(samples$diagnosis), "Unexpected diagnosis in age design.")
  samples$recorded_age_centered <- samples$Age - mean(samples$Age)
  design <- model.matrix(~ diagnosis + recorded_age_centered, samples)
  rownames(design) <- samples$donor_id
  assert(qr(design)$rank == ncol(design), "Capped-age sensitivity design is not identifiable.")
  assert(nrow(design) - ncol(design) >= 5, "Too few residual degrees of freedom.")
  design
}

fit_capped_age <- function(y, samples) {
  design <- make_age_design(samples)
  assert(identical(colnames(y), samples$donor_id), "Age-fit samples do not match count columns.")
  fit <- glmQLFit(y, design = design, dispersion = NULL, robust = TRUE,
                 abundance.trend = TRUE, legacy = FALSE,
                 winsor.tail.p = c(0.05, 0.1), prior.count = 0.125)
  test <- glmQLFTest(fit, coef = "diagnosisAD", poisson.bound = TRUE)
  tab <- topTags(test, n = Inf, sort.by = "none", adjust.method = "BH")$table
  assert(all(is.finite(as.matrix(tab[, c("logFC", "F", "PValue", "FDR")]))), "Invalid age-fit output.")
  assert(isTRUE(all.equal(tab$FDR, p.adjust(tab$PValue, "BH"), tolerance = 1e-13)), "Age-fit BH mismatch.")
  list(design = design, fit = fit, table = tab)
}

ranked_test <- function(statistic, indices, rho) {
  tab <- cameraPR(statistic, index = indices, use.ranks = TRUE,
                  inter.gene.cor = rho, sort = FALSE, directional = TRUE)
  tab$FDR <- p.adjust(tab$PValue, "BH")
  data.frame(gene_set = rownames(tab), tab, row.names = NULL,
             rank_biserial = vapply(indices, function(ix) rank_effect(statistic, ix), numeric(1)))
}

main <- function() {
  require_inferential_gate()
  out <- "results/sensitivity"
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  de <- read.delim("results/tables/differential_expression.tsv", check.names = FALSE)
  original <- read.delim("results/tables/enrichment_results.tsv", check.names = FALSE)
  ranked <- read.delim("results/tables/enrichment_ranking.tsv", check.names = FALSE)
  coverage <- read.delim("results/tables/enrichment_set_coverage.tsv", check.names = FALSE)
  samples <- read.delim("results/tables/pseudobulk_samples.tsv", check.names = FALSE)
  source_meta <- read.csv(gzfile("data/raw/GSE174367_snRNA-seq_cell_meta.csv.gz"),
                          check.names = FALSE, na.strings = c("NA", ""))
  source_donors <- unique(source_meta[, c("SampleID", "Diagnosis", "Age", "Sex", "PMI", "RIN", "Batch")])
  assert(nrow(source_donors) == nrow(samples) && !anyDuplicated(source_donors$SampleID),
         "Deposited donor clinical metadata are not unique.")
  source_donors <- source_donors[match(samples$donor_id, source_donors$SampleID), ]
  for (field in c("Age", "Sex", "PMI", "RIN", "Batch"))
    assert(isTRUE(all.equal(samples[[field]], source_donors[[field]], check.attributes = FALSE)),
           paste("Deposited metadata mismatch:", field))
  assert(identical(samples$diagnosis, source_donors$Diagnosis), "Diagnosis mismatch against deposit.")
  supplement <- read.delim("data/reference/supplementary_donor_metadata.tsv", check.names = FALSE)
  supplement <- supplement[match(samples$donor_id, supplement$Sample.ID), ]
  assert(!anyNA(supplement$Sample.ID), "Supplementary donor mapping missing.")
  donor_table <- data.frame(donor_id = samples$donor_id, diagnosis = samples$diagnosis,
    recorded_age_score = samples$Age, age_display = ifelse(samples$Age == 90, "90+", as.character(samples$Age)),
    age_top_coded = samples$Age == 90, sex = samples$Sex, PMI = samples$PMI, RIN = samples$RIN,
    astrocyte_nuclei = samples$astrocyte_nuclei, GEO_batch = samples$Batch,
    supplementary_library_batch = supplement$Library.Batch)
  write_tsv(donor_table, file.path(out, "donor_covariates.tsv"))
  covariates <- do.call(rbind, lapply(c("AD", "Control"), function(dx) {
    part <- samples[samples$diagnosis == dx, ]
    do.call(rbind, lapply(c("Age", "PMI", "RIN"), function(field) {
      v <- part[[field]]
      data.frame(diagnosis = dx, covariate = field, n = nrow(part), observed = sum(!is.na(v)),
                 missing = sum(is.na(v)), min = min(v, na.rm = TRUE), q1 = unname(quantile(v, .25, na.rm = TRUE)),
                 median = median(v, na.rm = TRUE), q3 = unname(quantile(v, .75, na.rm = TRUE)),
                 max = max(v, na.rm = TRUE), mean = mean(v, na.rm = TRUE))
    }))
  }))
  write_tsv(covariates, file.path(out, "clinical_covariate_summary.tsv"))
  sex_counts <- as.data.frame(table(diagnosis = samples$diagnosis, sex = samples$Sex), responseName = "donors")
  write_tsv(sex_counts, file.path(out, "sex_by_diagnosis.tsv"))
  batch_tables <- do.call(rbind, lapply(c("GEO_batch", "supplementary_library_batch"), function(nm) {
    z <- as.data.frame(table(diagnosis = donor_table$diagnosis, batch = donor_table[[nm]]), responseName = "donors")
    data.frame(source = nm, z)
  }))
  write_tsv(batch_tables, file.path(out, "batch_by_diagnosis.tsv"))

  # Save the data-design assessment before any new adjusted gene result is fitted.
  design <- make_age_design(samples)
  age_cor <- cor(as.numeric(samples$diagnosis == "AD"), samples$Age)
  shared_values <- intersect(samples$Age[samples$diagnosis == "AD"], samples$Age[samples$diagnosis == "Control"])
  ranges <- lapply(split(samples$Age, samples$diagnosis), range)
  common_range <- c(max(vapply(ranges, `[`, numeric(1), 1)), min(vapply(ranges, `[`, numeric(1), 2)))
  design_check <- list(status = "IDENTIFIABLE_WITH_LIMITED_OVERLAP", model = "~ diagnosis + recorded_age_centered",
    unit = "independent donor", n = nrow(samples), design_columns = colnames(design),
    rank = qr(design)$rank, residual_df = nrow(design) - qr(design)$rank,
    condition_number_centered_year_units = kappa(design, exact = TRUE),
    diagnosis_recorded_age_correlation = age_cor, predictor_VIF = 1/(1-age_cor^2),
    age_missing = sum(is.na(samples$Age)), age_center = mean(samples$Age),
    age_ranges = ranges, overlapping_recorded_range = common_range,
    exact_shared_recorded_values = shared_values,
    donors_in_exact_shared_values = as.list(table(factor(samples$diagnosis[samples$Age %in% shared_values], levels = c("Control", "AD")))),
    donors_in_overlapping_recorded_range = as.list(table(factor(samples$diagnosis[samples$Age >= common_range[1] & samples$Age <= common_range[2]], levels = c("Control", "AD")))),
    top_coded_donors = as.list(table(factor(samples$diagnosis[samples$Age == 90], levels = c("Control", "AD")))),
    missing_PMI_donors = samples$donor_id[is.na(samples$PMI)],
    representation = "Age=90 denotes the source 90-or-older category. The observed capped score is centered; no unknown chronological ages are imputed.",
    assumption = "Linear association with the observed capped score. This is a recorded-age-score sensitivity, not complete adjustment for true age.",
    limitations = "Exact shared age values consist only of the censored upper category; sex and conflicting batch identifiers remain unadjusted.")
  write_json_file(design_check, file.path(out, "design_preflight.json"))
  write_tsv(data.frame(donor_id = rownames(design), design, check.names = FALSE), file.path(out, "age_design_matrix.tsv"))
  implementation_notes <- c("# Installed cameraPR implementation check", "",
    paste("Inspected limma", as.character(packageVersion("limma")), "cameraPR.default and rankSumTestWithCorrelation, together with the installed camera help topic."), "",
    "The method default is inter.gene.cor=0.01. With use.ranks=TRUE, cameraPR calls rankSumTestWithCorrelation using infinite degrees of freedom, then forms a two-sided competitive P value from its two directional tails. Direction indicates the relative signed ranking of set genes against non-set genes. BH correction covers the tested collection. This workflow uses an externally fixed correlation scenario and does not estimate correlation from donor expression.", "",
    "The function default is the reason for the original 0.01 setting; it is not evidence that this cohort's true within-set correlation equals 0.01. Rank-biserial effect and direction remain fixed across the four correlation scenarios because the ranking and membership do not change.", "",
    "Source documentation: installed limma help(camera); Wu and Smyth (2012), DOI 10.1093/nar/gks461. Full implementation inspection output is retained with internal revision records, not redistributed as original project code.")
  writeLines(implementation_notes, file.path(out, "cameraPR_implementation_notes.md"))
  assert(identical(eval(formals(limma:::cameraPR.default)$inter.gene.cor), 0.01), "Installed cameraPR default differs from documented value.")

  # Freeze the original identifiers, ranking and set membership for all scenarios.
  assert(!anyDuplicated(ranked$gene_id) && !anyDuplicated(ranked$gene_symbol), "Nonunique frozen ranking identifiers.")
  statistic <- setNames(ranked$ranking_statistic, ranked$gene_symbol)
  assert(isTRUE(all.equal(statistic, setNames(sign(ranked$logFC)*sqrt(ranked$F), ranked$gene_symbol), tolerance = 1e-12)), "Frozen ranking does not match primary statistics.")
  sets <- read_gmt("data/raw/h.all.v2026.1.Hs.symbols.gmt")
  indices <- lapply(sets, function(s) which(names(statistic) %in% s))
  eligible_names <- coverage$gene_set[coverage$tested]
  indices <- indices[eligible_names]
  assert(length(indices) == nrow(original) && all(lengths(indices) == coverage$mapped_genes[match(names(indices), coverage$gene_set)]), "Frozen mapped gene sets differ from primary analysis.")
  reproduced <- ranked_test(statistic, indices, 0.01)
  matched <- reproduced[match(original$gene_set, reproduced$gene_set), ]
  for (field in c("NGenes", "PValue", "FDR", "rank_biserial"))
    assert(isTRUE(all.equal(matched[[field]], original[[field]], tolerance = 1e-11)), paste("Primary enrichment reproduction mismatch", field))
  assert(identical(matched$Direction, original$Direction), "Primary enrichment direction mismatch.")
  original_five <- original$gene_set[original$FDR < .05]
  assert(length(original_five) == 5, "Frozen primary discovery count differs from baseline.")
  all_rho <- do.call(rbind, lapply(correlation_scenarios, function(rho) {
    z <- ranked_test(statistic, indices, rho); z$correlation <- rho
    z$significant <- z$FDR < .05; z$original_five <- z$gene_set %in% original_five; z
  }))
  assert(all(is.finite(all_rho$PValue)) && all(all_rho$PValue >= 0 & all_rho$PValue <= 1), "Empirical scenario P value outside [0,1].")
  write_tsv(all_rho, file.path(out, "correlation_all_sets.tsv"))
  corr_counts <- do.call(rbind, lapply(correlation_scenarios, function(rho) {
    z <- all_rho[all_rho$correlation == rho, ]
    data.frame(correlation = rho, tested_sets = nrow(z), significant_sets = sum(z$significant),
               original_five_significant = sum(z$significant & z$original_five),
               significant_up = sum(z$significant & z$Direction == "Up"), significant_down = sum(z$significant & z$Direction == "Down"))
  }))
  compact <- data.frame(gene_set = original_five, Direction = original$Direction[match(original_five, original$gene_set)],
                         rank_biserial = original$rank_biserial[match(original_five, original$gene_set)])
  for (rho in correlation_scenarios) {
    z <- all_rho[all_rho$correlation == rho, ]
    compact[[paste0("FDR_rho_", format(rho, trim = TRUE))]] <- z$FDR[match(original_five, z$gene_set)]
  }
  write_tsv(corr_counts, file.path(out, "correlation_summary.tsv"))
  write_tsv(compact, file.path(out, "correlation_original_five.tsv"))
  write_tsv(data.frame(gene_set = rep(names(indices), lengths(indices)), gene_symbol = unlist(lapply(indices, function(ix) names(statistic)[ix]), use.names = FALSE)),
            file.path(out, "frozen_gene_set_membership.tsv"))

  # The only new model. TMM, expression filtering and symbol mapping stay frozen.
  saved <- readRDS("results/tables/primary_fit.rds")
  assert(identical(saved$input$samples$donor_id, samples$donor_id), "Saved primary sample order differs.")
  y <- saved$analysis$y
  assert(setequal(rownames(y), de$gene_id), "Saved primary gene universe differs.")
  norm_before <- y$samples$norm.factors
  adjusted <- fit_capped_age(y, samples)
  assert(identical(norm_before, y$samples$norm.factors), "Sensitivity changed TMM factors.")
  ade <- adjusted$table
  ade$n_control <- sum(samples$diagnosis == "Control"); ade$n_ad <- sum(samples$diagnosis == "AD")
  ade$contrast <- "AD-Control"; ade$model <- "~ diagnosis + recorded_age_centered"
  ade <- ade[order(ade$PValue, ade$gene_id, method = "radix"), ]
  write_tsv(ade, file.path(out, "age_score_differential_expression.tsv"))
  canonical <- ade[match(ranked$gene_id, ade$gene_id), ]
  assert(identical(canonical$gene_id, ranked$gene_id) && identical(canonical$gene_symbol, ranked$gene_symbol), "Canonical symbol mapping changed.")
  adjusted_statistic <- setNames(sign(canonical$logFC)*sqrt(canonical$F), ranked$gene_symbol)
  ae <- ranked_test(adjusted_statistic, indices, .01)
  ae <- ae[order(ae$PValue, ae$gene_set, method = "radix"), ]
  write_tsv(ae, file.path(out, "age_score_enrichment.tsv"))
  write_tsv(data.frame(gene_id = canonical$gene_id, gene_symbol = canonical$gene_symbol, ranking_statistic = adjusted_statistic),
            file.path(out, "age_score_ranking.tsv"))
  comparison <- data.frame(gene_set = original_five,
    primary_direction = original$Direction[match(original_five, original$gene_set)],
    primary_PValue = original$PValue[match(original_five, original$gene_set)],
    primary_FDR = original$FDR[match(original_five, original$gene_set)],
    primary_rank_biserial = original$rank_biserial[match(original_five, original$gene_set)],
    adjusted_direction = ae$Direction[match(original_five, ae$gene_set)],
    adjusted_PValue = ae$PValue[match(original_five, ae$gene_set)],
    adjusted_FDR = ae$FDR[match(original_five, ae$gene_set)],
    adjusted_rank_biserial = ae$rank_biserial[match(original_five, ae$gene_set)])
  write_tsv(comparison, file.path(out, "age_score_original_five.tsv"))
  emx_primary <- de[de$gene_symbol == "EMX2", ]; emx_adjusted <- ade[ade$gene_symbol == "EMX2", ]
  assert(nrow(emx_primary) == 1 && nrow(emx_adjusted) == 1, "EMX2 mapping ambiguous.")
  emx <- data.frame(model = c("Primary diagnosis-only", "Recorded capped-age score sensitivity"),
    gene_id = c(emx_primary$gene_id, emx_adjusted$gene_id),
    logFC = c(emx_primary$logFC, emx_adjusted$logFC), F = c(emx_primary$F, emx_adjusted$F),
    PValue = c(emx_primary$PValue, emx_adjusted$PValue), FDR = c(emx_primary$FDR, emx_adjusted$FDR))
  write_tsv(emx, file.path(out, "emx2_model_comparison.tsv"))
  all_directions_consistent <- all(vapply(split(all_rho$Direction, all_rho$gene_set), function(v) length(unique(v)) == 1, logical(1)))
  summary <- list(status = "EXECUTED", original_method_unchanged = TRUE,
    source = "scripts/07_sensitivity.R", correlation_scenarios = correlation_scenarios,
    primary_enrichment_reproduced = TRUE, frozen_ranked_genes = length(statistic), frozen_gene_sets = length(indices),
    cameraPR_default_correlation = .01, scenario_interpretation = "Sensitivity settings, not estimated correlations",
    correlation_results = corr_counts, original_five_correlation_results = compact,
    correlation_direction_consistency_all_sets = all_directions_consistent,
    adjusted_model = design_check, adjusted_tested_genes = nrow(ade), adjusted_significant_genes = sum(ade$FDR < .05),
    adjusted_min_gene_fdr = min(ade$FDR), adjusted_significant_sets = sum(ae$FDR < .05),
    original_five_significant_after_age_score_adjustment = sum(comparison$adjusted_FDR < .05),
    original_five_directions_preserved_after_age_score_adjustment = sum(comparison$primary_direction == comparison$adjusted_direction),
    original_five_age_score_results = comparison, emx2_model_comparison = emx,
    adjusted_NB_dispersion = as.numeric(adjusted$fit$dispersion),
    adjusted_NB_top_proportion = adjusted$fit$top.proportion,
    versions = list(R = as.character(getRversion()), edgeR = as.character(packageVersion("edgeR")), limma = as.character(packageVersion("limma"))))
  write_json_file(summary, "results/audit/sensitivity_validation.json")
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  assessment <- c("# Confounding assessment", "",
    "Donor identifiers and clinical covariates were checked directly against deposited cell metadata. No donor was excluded. A recorded age of 90 represents the source upper category (90 or older); exact ages above that threshold are unavailable.", "",
    sprintf("The sole supplementary model, ~ diagnosis + recorded_age_centered, has rank %d and %d residual degrees of freedom. Age is centered at %.6f recorded years. The diagnosis-age correlation is %.6f, the predictor VIF %.6f and the centered design condition number %.6f.", qr(design)$rank, nrow(design)-ncol(design), mean(samples$Age), age_cor, 1/(1-age_cor^2), kappa(design, exact=TRUE)), "",
    sprintf("Recorded age ranges overlap from %g to %g, but the only exactly shared recorded value is %s. That upper category contains %d AD and %d control donors. The fit therefore checks a linear capped-score specification; it cannot identify unobserved within-category age differences or establish adequate control of age confounding.", common_range[1],common_range[2],paste(shared_values,collapse=", "),sum(samples$Age==90 & samples$diagnosis=="AD"),sum(samples$Age==90 & samples$diagnosis=="Control")), "",
    sprintf("Sex and both source batch fields were tabulated but not included in additional fits. Female/male counts are %d/%d in AD and %d/%d in controls. GEO batch has both diagnoses in every level; the alternative Library.Batch assignment differs and its assay attribution remains unresolved. Neither table is silently substituted for the other. PMI is missing for %s; age, sex, RIN and GEO batch are complete. Full donor values and distribution summaries accompany this assessment.", sum(samples$diagnosis=="AD" & samples$Sex=="F"),sum(samples$diagnosis=="AD" & samples$Sex=="M"),sum(samples$diagnosis=="Control" & samples$Sex=="F"),sum(samples$diagnosis=="Control" & samples$Sex=="M"),paste(samples$donor_id[is.na(samples$PMI)],collapse=", ")), "",
    "The primary diagnosis-only analysis, gene-filtering decision, TMM factors and original symbol/set mapping are preserved. The adjusted coefficient and enrichment ranking are a supplementary sensitivity result, not a replacement primary finding.")
  writeLines(assessment, file.path(out, "confounding_assessment.md"))
  cat("Correlation scenario counts:\n"); print(corr_counts, row.names=FALSE)
  cat(sprintf("Capped-age-score sensitivity: %d genes tested; %d gene discoveries; %d/%d Hallmark sets below BH .05; %d/5 original sets remain below BH .05.\n", nrow(ade),sum(ade$FDR<.05),sum(ae$FDR<.05),nrow(ae),sum(comparison$adjusted_FDR<.05)))
  print(emx, row.names=FALSE)
}

if (sys.nframe() == 0L) main()

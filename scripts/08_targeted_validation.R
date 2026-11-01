#!/usr/bin/env Rscript
# Bounded, donor-level diagnostics from distributed pseudobulk inputs.
source("scripts/r_common.R")

agree_numeric <- function(actual, reference, label, atol = 1e-10, rtol = 1e-8) {
  assert(length(actual) == length(reference) && all(is.finite(actual)) &&
           all(is.finite(reference)), paste("Invalid comparison:", label))
  error <- abs(actual - reference)
  assert(all(error <= atol + rtol * abs(reference)), paste("Baseline mismatch:", label))
  data.frame(comparison = label, n = length(error), max_absolute_error = max(error),
             atol = atol, rtol = rtol)
}

align_unique <- function(tab, key, ids) {
  assert(!anyNA(tab[[key]]) && !anyDuplicated(tab[[key]]) &&
           setequal(tab[[key]], ids), paste("Identifier mismatch:", key))
  tab[match(ids, tab[[key]]), , drop = FALSE]
}

fixed_ranking <- function(tab, canonical) {
  assert(!anyDuplicated(tab$gene_id) && !anyDuplicated(canonical$gene_id) &&
           !anyDuplicated(canonical$gene_symbol), "Nonunique canonical mapping.")
  j <- match(canonical$gene_id, tab$gene_id)
  assert(!anyNA(j), "Canonical gene missing from omission fit.")
  z <- tab[j, , drop = FALSE]
  assert(identical(z$gene_symbol, canonical$gene_symbol), "Canonical symbol changed.")
  assert(all(is.finite(z$F)) && all(z$F >= 0) && all(is.finite(z$logFC)),
         "Nonfinite or negative ranking input.")
  setNames(sign(z$logFC) * sqrt(z$F), canonical$gene_symbol)
}

frozen_indices <- function(canonical, membership, original) {
  assert(!anyDuplicated(paste(membership$gene_set, membership$gene_symbol)),
         "Duplicate frozen set membership.")
  assert(setequal(unique(membership$gene_set), original$gene_set), "Gene-set family changed.")
  ix <- lapply(original$gene_set, function(g) {
    j <- match(membership$gene_symbol[membership$gene_set == g], canonical$gene_symbol)
    assert(!anyNA(j), "Unmapped frozen member.")
    j
  })
  names(ix) <- original$gene_set
  assert(all(lengths(ix) == original$NGenes) && all(lengths(ix) >= 15 & lengths(ix) <= 500),
         "Frozen set sizes changed.")
  ix
}

program_test <- function(statistic, indices) {
  z <- cameraPR(statistic, index = indices, use.ranks = TRUE,
                inter.gene.cor = 0.01, sort = FALSE, directional = TRUE)
  z$FDR <- p.adjust(z$PValue, "BH")
  ans <- data.frame(gene_set = rownames(z), z, row.names = NULL,
                    rank_biserial = vapply(indices, function(i) rank_effect(statistic, i), numeric(1)))
  assert(all(is.finite(ans$PValue)) && all(ans$PValue >= 0 & ans$PValue <= 1), "Invalid set P values.")
  assert(all(abs(ans$rank_biserial) <= 1), "Invalid set effect.")
  assert(all(abs(ans$rank_biserial) < 1e-14 |
               (ans$Direction == "Up") == (ans$rank_biserial > 0)), "Set direction mismatch.")
  ans$iteration_rank <- match(seq_len(nrow(ans)), order(ans$PValue, ans$gene_set, method = "radix"))
  ans
}

summarize_programs <- function(long) {
  out <- lapply(split(long, long$gene_set), function(z) {
    j <- order(-z$absolute_effect_change, z$omitted_donor, method = "radix")[1]
    data.frame(gene_set = z$gene_set[1], primary_original_five = z$primary_FDR[1] < .05,
      primary_rank_biserial = z$primary_rank_biserial[1], primary_FDR = z$primary_FDR[1],
      primary_rank = z$primary_rank[1], n_iterations = nrow(z),
      direction_consistent_n = sum(z$direction_matches_primary), FDR_below_n = sum(z$FDR < .05),
      min_effect = min(z$rank_biserial), max_effect = max(z$rank_biserial),
      min_FDR = min(z$FDR), max_FDR = max(z$FDR),
      min_rank = min(z$iteration_rank), max_rank = max(z$iteration_rank),
      max_absolute_effect_change = z$absolute_effect_change[j],
      most_influential_omission = z$omitted_donor[j])
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out[order(out$primary_rank), ]
}

design_diagnostics <- function(donors) {
  fields <- c("recorded_age_score", "sex", "GEO_batch", "supplementary_library_batch", "RIN", "PMI")
  do.call(rbind, lapply(fields, function(field) {
    use <- !is.na(donors[[field]])
    d <- donors[use, , drop = FALSE]
    value <- d[[field]]
    categorical <- field %in% c("sex", "GEO_batch", "supplementary_library_batch")
    value <- if (categorical) factor(value) else as.numeric(scale(value))
    dx <- factor(d$diagnosis, levels = c("Control", "AD"))
    estimable <- length(unique(value)) > 1 && all(is.finite(as.numeric(value))) && nlevels(droplevels(dx)) == 2
    design <- if (estimable) model.matrix(~ dx + value) else NULL
    rank <- if (estimable) qr(design)$rank else NA_integer_
    tab <- if (categorical) table(dx, value) else NULL
    data.frame(covariate = field, n = nrow(d), missing = sum(!use),
      n_ad = sum(dx == "AD"), n_control = sum(dx == "Control"),
      design_columns = if (estimable) ncol(design) else NA_integer_, rank = rank,
      residual_df = if (estimable) nrow(design) - rank else NA_integer_,
      full_rank = if (estimable) rank == ncol(design) else FALSE,
      condition_number = if (estimable) kappa(design, exact = TRUE) else NA_real_,
      levels_shared_by_both_groups = if (categorical) sum(colSums(tab > 0) == 2) else NA_integer_,
      numeric_standardized_for_condition_number = !categorical)
  }))
}

donor_structure <- function(y, donors) {
  assert(identical(colnames(y), donors$donor_id), "MDS donor order mismatch.")
  m <- plotMDS(y, top = 500, gene.selection = "common", method = "logFC", prior.count = 2, plot = FALSE)
  coords <- cbind(m$x, m$y)
  # Axis signs are arbitrary; orient deterministically without using diagnosis.
  for (k in 1:2) if (coords[which.max(abs(coords[, k])), k] < 0) coords[, k] <- -coords[, k]
  lcpm <- cpm(y, log = TRUE, prior.count = 2)
  variance <- rowMeans((lcpm - rowMeans(lcpm))^2)
  ids <- order(variance, decreasing = TRUE)[seq_len(min(500, nrow(lcpm)))]
  list(coordinates = data.frame(donors, MDS1 = coords[, 1], MDS2 = coords[, 2]),
       genes = data.frame(gene_id = rownames(lcpm)[ids], variance = variance[ids], selection_rank = seq_along(ids)),
       summary = list(n_donors = ncol(y), top_genes = length(ids), gene_selection = "common",
         method = "logFC", normalization = "TMM", prior_count = 2,
         variance_percent = unname(100 * m$var.explained[1:2]),
         interpretation = "Descriptive common-gene MDS (PCA geometry), not a diagnosis test or a confounder adjustment"))
}

check_runtime <- function() {
  lock <- jsonlite::fromJSON("environment/renv.lock", simplifyVector = FALSE)
  assert(as.character(getRversion()) == lock$R$Version, "R version differs from analysis lock.")
  for (p in names(lock$Packages))
    assert(packageDescription(p)$Version == lock$Packages[[p]]$Version, paste("Version mismatch:", p))
}

main <- function() {
  started <- Sys.time()
  check_runtime(); require_inferential_gate()
  out <- "results/targeted_validation"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_json_file(list(status = "RUNNING", started_utc = format(started, tz = "UTC", usetz = TRUE)),
                  file.path(out, "validation_summary.json"))
  input_paths <- c("results/tables/pseudobulk_counts.tsv.gz", "results/tables/pseudobulk_samples.tsv",
    "results/tables/gene_annotation.tsv", "results/tables/gene_filtering.tsv", "results/tables/normalization.tsv",
    "results/tables/differential_expression.tsv", "results/tables/enrichment_ranking.tsv",
    "results/tables/enrichment_results.tsv", "results/tables/lodo_all_genes.tsv.gz",
    "results/sensitivity/frozen_gene_set_membership.tsv", "results/sensitivity/donor_covariates.tsv",
    "data/reference/supplementary_donor_metadata.tsv", "results/audit/dataset_audit.json",
    "environment/renv.lock", "scripts/r_common.R", "scripts/08_targeted_validation.R",
    "docs/TARGETED_VALIDATION_PLAN.md")
  before_hash <- tools::sha256sum(input_paths)
  write_tsv(data.frame(path = input_paths, sha256 = unname(before_hash)), file.path(out, "input_code_checksums.tsv"))
  x <- read_inputs()
  original <- read.delim("results/tables/differential_expression.tsv", check.names = FALSE)
  canonical <- read.delim("results/tables/enrichment_ranking.tsv", check.names = FALSE)
  primary_sets <- read.delim("results/tables/enrichment_results.tsv", check.names = FALSE)
  primary_sets <- primary_sets[order(primary_sets$PValue, primary_sets$gene_set, method = "radix"), ]
  membership <- read.delim("results/sensitivity/frozen_gene_set_membership.tsv", check.names = FALSE)
  indices <- frozen_indices(canonical, membership, primary_sets)
  old_lodo <- read.delim(gzfile("results/tables/lodo_all_genes.tsv.gz"), check.names = FALSE)
  donors <- align_unique(read.delim("results/sensitivity/donor_covariates.tsv"), "donor_id", x$samples$donor_id)
  assert(identical(donors$diagnosis, x$samples$diagnosis) && identical(donors$sex, x$samples$Sex),
         "Clinical labels disagree with pseudobulk metadata.")
  for (field in c("PMI", "RIN", "astrocyte_nuclei"))
    assert(isTRUE(all.equal(donors[[field]], x$samples[[field]], check.attributes = FALSE)), paste("Clinical mismatch:", field))
  assert(identical(donors$recorded_age_score, x$samples$Age) &&
           identical(donors$GEO_batch, x$samples$Batch), "Age/batch metadata mismatch.")
  supp <- align_unique(read.delim("data/reference/supplementary_donor_metadata.tsv"), "Sample.ID", x$samples$donor_id)
  assert(identical(donors$supplementary_library_batch, supp$Library.Batch), "Supplementary batch mismatch.")
  cat("Refitting unchanged full-data primary model from distributed pseudobulk counts.\n")
  a <- fit_primary(x$counts, x$samples, x$annotation)
  filtering <- align_unique(read.delim("results/tables/gene_filtering.tsv"), "gene_id", x$annotation$gene_id)
  assert(identical(unname(a$keep), filtering$retained), "Primary gene filter changed.")
  ref <- align_unique(original, "gene_id", a$table$gene_id)
  comparisons <- lapply(c("logFC", "logCPM", "F", "PValue", "FDR"),
    function(f) agree_numeric(a$table[[f]], ref[[f]], paste("primary", f)))
  assert(identical(a$table$FDR < .05, ref$FDR < .05), "Primary discovery decisions changed.")
  norm <- align_unique(read.delim("results/tables/normalization.tsv"), "donor_id", x$samples$donor_id)
  comparisons[[length(comparisons)+1]] <- agree_numeric(a$y$samples$norm.factors, norm$TMM_factor, "primary TMM")
  mapped <- unique_ranked_genes(a$table)
  assert(identical(mapped$gene_id, canonical$gene_id), "Primary canonical gene selection changed.")
  base_stat <- fixed_ranking(a$table, canonical)
  comparisons[[length(comparisons)+1]] <- agree_numeric(unname(base_stat), canonical$ranking_statistic, "primary ranking")
  base_sets <- program_test(base_stat, indices)
  for (f in c("NGenes", "PValue", "FDR", "rank_biserial"))
    comparisons[[length(comparisons)+1]] <- agree_numeric(base_sets[[f]], primary_sets[[f]], paste("primary sets", f))
  assert(identical(base_sets$Direction, primary_sets$Direction) &&
           identical(base_sets$FDR < .05, primary_sets$FDR < .05), "Primary set decisions changed.")
  structure <- donor_structure(a$y, donors)
  write_tsv(structure$coordinates, file.path(out, "donor_mds.tsv"))
  write_tsv(structure$genes, file.path(out, "mds_genes.tsv"))
  write_json_file(structure$summary, file.path(out, "mds_summary.json"))
  write_tsv(design_diagnostics(donors), file.path(out, "design_diagnostics.tsv"))
  write_tsv(base_sets, file.path(out, "primary_gene_sets_recomputed.tsv"))
  long <- vector("list", nrow(donors)); rankings <- vector("list", nrow(donors))
  for (i in seq_len(nrow(donors))) {
    donor <- donors$donor_id[i]
    fit <- fit_primary(x$counts[, -i, drop = FALSE], x$samples[-i, , drop = FALSE], x$annotation, keep = a$keep)
    z <- fit$table
    assert(identical(z$gene_id, a$table$gene_id), "Omission gene universe changed.")
    old <- align_unique(old_lodo[old_lodo$omitted_donor == donor, ], "gene_id", z$gene_id)
    for (f in c("logFC", "PValue", "FDR"))
      comparisons[[length(comparisons)+1]] <- agree_numeric(z[[f]], old[[f]], paste(donor, f))
    statistic <- fixed_ranking(z, canonical)
    sets <- program_test(statistic, indices)
    sets$omitted_donor <- donor; sets$omitted_diagnosis <- donors$diagnosis[i]
    sets$primary_rank_biserial <- primary_sets$rank_biserial
    sets$primary_FDR <- primary_sets$FDR; sets$primary_rank <- seq_len(nrow(primary_sets))
    sets$direction_matches_primary <- sign(sets$rank_biserial) == sign(sets$primary_rank_biserial)
    sets$absolute_effect_change <- abs(sets$rank_biserial - sets$primary_rank_biserial)
    long[[i]] <- sets
    zz <- z[match(canonical$gene_id, z$gene_id), c("gene_id", "gene_symbol", "logFC", "F", "PValue", "FDR")]
    rankings[[i]] <- data.frame(omitted_donor = donor, zz, ranking_statistic = unname(statistic))
    cat(sprintf("Validated omission %02d/%02d: %s; %d genes; %d sets; %d sets below FDR .05.\n",
                i, nrow(donors), donor, nrow(z), nrow(sets), sum(sets$FDR < .05)))
  }
  long <- do.call(rbind, long); rownames(long) <- NULL
  assert(nrow(long) == nrow(donors)*nrow(primary_sets) && !anyDuplicated(paste(long$gene_set, long$omitted_donor)),
         "Incomplete/duplicated set-by-omission results.")
  summary <- summarize_programs(long)
  write_tsv(long, file.path(out, "gene_set_lodo.tsv"))
  write_tsv(summary, file.path(out, "gene_set_lodo_summary.tsv"))
  con <- gzfile(file.path(out, "lodo_gene_rankings.tsv.gz"), "wt")
  write_tsv(do.call(rbind, rankings), con); close(con)
  comparisons <- do.call(rbind, comparisons)
  write_tsv(comparisons, file.path(out, "baseline_agreement.tsv"))
  assert(identical(before_hash, tools::sha256sum(input_paths)), "A source input changed during validation.")
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  write_json_file(list(status = "PASS", started_utc = format(started, tz = "UTC", usetz = TRUE),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    command = "Rscript --vanilla scripts/08_targeted_validation.R", unit = "independent donor",
    model = "~ diagnosis", contrast = "AD-Control", full_data_fits = 1, omission_fits = nrow(donors),
    tested_genes = nrow(a$table), canonical_symbols = nrow(canonical), tested_sets = nrow(primary_sets),
    gene_set_omission_rows = nrow(long), correlation = .01, multiple_testing = "BH over all 49 sets separately per omission",
    baseline_agreement_max_absolute_error = max(comparisons$max_absolute_error),
    original_five = summary[summary$primary_original_five, ],
    mds = structure$summary, numerical_threads = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))),
    raw_data_reaggregated = FALSE, external_validation = FALSE,
    limitations = "Conditioned on fixed full-data filtering and mapping and correlation 0.01. Descriptive internal influence, not independent validation, confidence intervals, or causal adjustment."),
    file.path(out, "validation_summary.json"))
  cat("PASS: unchanged baseline and 18 donor omissions validated from supplied donor counts.\n")
}

if (sys.nframe() == 0L) main()

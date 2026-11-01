#!/usr/bin/env Rscript
source("scripts/r_common.R")
x <- read_inputs()
a <- fit_primary(x$counts, x$samples, x$annotation)
de <- a$table
# Explicit sample counts and direction make the standalone result interpretable.
de$n_control <- sum(x$samples$diagnosis == "Control")
de$n_ad <- sum(x$samples$diagnosis == "AD")
de$contrast <- "AD-Control"
de <- de[order(de$PValue, de$gene_id, method = "radix"), ]
write_tsv(de, "results/tables/differential_expression.tsv")
write_tsv(data.frame(gene_id = x$annotation$gene_id,
                     gene_symbol = x$annotation$gene_symbol,
                     retained = a$keep,
                     total_count = rowSums(x$counts)),
          "results/tables/gene_filtering.tsv")
normalization <- data.frame(donor_id = colnames(a$y),
                            diagnosis = x$samples$diagnosis,
                            astrocyte_nuclei = x$samples$astrocyte_nuclei,
                            unfiltered_library_size = colSums(x$counts),
                            filtered_library_size = a$y$samples$lib.size,
                            TMM_factor = a$y$samples$norm.factors,
                            effective_library_size = a$y$samples$lib.size * a$y$samples$norm.factors)
write_tsv(normalization, "results/tables/normalization.tsv")
write_tsv(data.frame(donor_id = rownames(a$design), a$design, check.names = FALSE),
          "results/tables/design_matrix.tsv")
write_tsv(data.frame(gene_id = rownames(a$fit),
                     NB_dispersion = rep_len(a$fit$dispersion, nrow(a$fit)),
                     QL_posterior_dispersion = a$fit$s2.post,
                     adjusted_residual_df = a$fit$df.residual.adj,
                     prior_df = a$fit$df.prior),
          "results/tables/dispersion_diagnostics.tsv")
selected <- de[seq_len(min(10L, nrow(de))), c("gene_id", "gene_symbol", "logFC", "PValue", "FDR")]
selected$selection_rank <- seq_len(nrow(selected))
write_tsv(selected, "results/tables/lodo_selected_genes.tsv")
summary <- list(status = "EXECUTED", analysis = "donor-level astrocyte pseudobulk",
  model = "~ diagnosis", reference_level = "Control", contrast = "AD-Control",
  positive_logFC = "higher expression in AD", unit = "independent donor",
  donors = ncol(x$counts), n_ad = sum(x$samples$diagnosis == "AD"),
  n_control = sum(x$samples$diagnosis == "Control"),
  astrocyte_nuclei = sum(x$samples$astrocyte_nuclei),
  input_genes = nrow(x$counts), tested_genes = nrow(de),
  significant_genes = sum(de$FDR < 0.05),
  significant_up = sum(de$FDR < 0.05 & de$logFC > 0),
  significant_down = sum(de$FDR < 0.05 & de$logFC < 0),
  significant_abs_logFC_ge_0_25 = sum(de$FDR < 0.05 & abs(de$logFC) >= 0.25),
  min_raw_p = min(de$PValue), min_fdr = min(de$FDR),
  filtering = list(function_name = "edgeR::filterByExpr", min_count = 10,
                   min_total_count = 15, large_n = 10, min_prop = 0.7,
                   design = "~ diagnosis", reset_library_sizes = TRUE),
  normalization = list(method = "TMM", logratioTrim = 0.3, sumTrim = 0.05,
                       doWeighting = TRUE, Acutoff = -1e10),
  inference = list(method = "edgeR v4 negative-binomial quasi-likelihood F-test",
                   coefficient = "diagnosisAD", numerator_df = 1,
                   nominal_residual_df = ncol(x$counts) - ncol(a$design),
                   robust = TRUE, legacy = FALSE, abundance_trend = TRUE,
                   winsor_tail_p = c(0.05, 0.1), prior_count_for_logFC = 0.125,
                   NB_dispersion_estimation = "glmQLFit automatic from high-abundance genes",
                   NB_dispersion = as.numeric(a$fit$dispersion),
                   NB_top_proportion = a$fit$top.proportion,
                   QL_average_dispersion = a$fit$average.ql.dispersion,
                   poisson_bound = TRUE, p_value_alternative = "two-sided (one-coefficient F test)",
                   multiple_testing = "Benjamini-Hochberg across all tested genes",
                   fdr_threshold = 0.05),
  selected_gene_rule = "10 smallest full-data raw P values; ties ordered by gene_id",
  interpretation = "Unadjusted observational group association; age, sex and technical confounding may remain",
  versions = list(R = as.character(getRversion()), edgeR = as.character(packageVersion("edgeR")),
                  limma = as.character(packageVersion("limma")), statmod = as.character(packageVersion("statmod"))))
write_json_file(summary, "results/audit/analysis_summary.json")
saveRDS(list(input = x, analysis = a), "results/tables/primary_fit.rds", compress = "gzip")
writeLines(capture.output(sessionInfo()), "results/audit/R_sessionInfo.txt")
cat(sprintf("DE complete: %d genes, %d AD and %d control donors; %d genes at BH FDR < 0.05.\n",
            nrow(de), summary$n_ad, summary$n_control, summary$significant_genes))

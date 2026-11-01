#!/usr/bin/env Rscript
source("scripts/r_common.R")
# Synthetic counts test implementation only; never report these as real results.
set.seed(174367)
ng <- 1500L; ns <- 12L
samples <- data.frame(donor_id = sprintf("synthetic_%02d", seq_len(ns)),
                       diagnosis = rep(c("Control", "AD"), each = 6), astrocyte_nuclei = 100)
annotation <- data.frame(gene_id = sprintf("test_gene_%04d", seq_len(ng)),
                          gene_symbol = sprintf("TEST%04d", seq_len(ng)))
mu <- matrix(rep(exp(rnorm(ng, log(100), 0.5)), ns), nrow = ng)
mu[1:20, 7:12] <- mu[1:20, 7:12] * 8
mu[21:40, 7:12] <- mu[21:40, 7:12] / 8
counts <- matrix(rnbinom(ng * ns, mu = mu, size = 30), nrow = ng,
                 dimnames = list(annotation$gene_id, samples$donor_id))
validate_inputs(counts, samples, annotation)
a <- fit_primary(counts, samples, annotation)
assert(nrow(a$table) == ng, "Synthetic filtering unexpectedly dropped adequately expressed genes.")
assert(all(a$table$logFC[1:20] > 0) && all(a$table$logFC[21:40] < 0), "Contrast direction failed.")
assert(sum(a$table$FDR[1:40] < 0.05) >= 38, "Synthetic strong signals not recovered.")
assert(isTRUE(all.equal(a$table$FDR, p.adjust(a$table$PValue, "BH"))), "BH FDR mismatch.")
assert(inherits(try(validate_inputs(counts[, ns:1], samples, annotation), silent = TRUE), "try-error"),
       "Sample misalignment was not rejected.")
dup <- samples; dup$donor_id[2] <- dup$donor_id[1]
assert(inherits(try(validate_inputs(counts, dup, annotation), silent = TRUE), "try-error"),
       "Duplicate biological donor was not rejected.")
negative <- counts; negative[1, 1] <- -1
assert(inherits(try(validate_inputs(negative, samples, annotation), silent = TRUE), "try-error"),
       "Negative counts were not rejected.")
b <- fit_primary(counts[, -1], samples[-1, ], annotation, keep = a$keep)
assert(identical(a$table$gene_id, b$table$gene_id), "Fixed LODO universe changed.")
assert(all(b$table$logFC[1:20] > 0), "LODO contrast orientation failed.")
# Verify deterministic duplicate handling does not select by significance.
z <- a$table[c(1, 2), ]; z$gene_symbol <- "DUPLICATE"; z$logCPM <- c(1, 10); z$PValue <- c(1e-20, 0.5)
u <- unique_ranked_genes(z)
assert(nrow(u) == 1 && u$gene_id == z$gene_id[2], "Duplicate symbols not selected by mean abundance.")
statistic <- seq(-5, 5, length.out = 1000)
ix <- list(positive = 901:1000, negative = 1:100)
cam <- cameraPR(statistic, ix, use.ranks = TRUE, inter.gene.cor = 0.01, sort = FALSE)
assert(identical(cam$Direction, c("Up", "Down")) && all(cam$FDR < 0.05), "Rank-test directions failed.")
assert(abs(rank_effect(statistic, ix$positive) - 1) < 1e-12 &&
         abs(rank_effect(statistic, ix$negative) + 1) < 1e-12, "Rank-biserial effect failed.")
cat("PASS: 12 statistical checks (synthetic fit, sign orientation, signal recovery, BH correction, alignment rejection, duplicate rejection, count validity, fixed LODO universe, LODO orientation, symbol selection, enrichment direction, rank effect).\n")

if ("--empirical" %in% commandArgs(trailingOnly = TRUE)) {
  # Independent arithmetic checks of the exported empirical artifacts.
  de <- read.delim("results/tables/differential_expression.tsv", check.names = FALSE)
  en <- read.delim("results/tables/enrichment_results.tsv", check.names = FALSE)
  lo <- read.delim("results/tables/lodo_selected_results.tsv", check.names = FALSE)
  selected <- read.delim("results/tables/lodo_selected_genes.tsv", check.names = FALSE)
  ls <- read.delim("results/tables/lodo_selected_summary.tsv", check.names = FALSE)
  iterations <- read.delim("results/tables/lodo_iterations.tsv", check.names = FALSE)
  norm <- read.delim("results/tables/normalization.tsv", check.names = FALSE)
  input <- read_inputs()
  primary <- readRDS("results/tables/primary_fit.rds")
  bh_manual <- function(p) {
    o <- order(p); n <- length(p)
    adjusted <- pmin(1, rev(cummin(rev(p[o] * n/seq_len(n)))))
    out <- numeric(n); out[o] <- adjusted; out
  }
  assert(isTRUE(all.equal(de$FDR, bh_manual(de$PValue), tolerance = 1e-12)), "Empirical gene BH mismatch.")
  assert(isTRUE(all.equal(en$FDR, bh_manual(en$PValue), tolerance = 1e-12)), "Empirical set BH mismatch.")
  assert(nrow(de) == sum(primary$analysis$keep), "Complete tested gene table missing rows.")
  assert(identical(selected$gene_id, head(de$gene_id[order(de$PValue, de$gene_id)], 10)),
         "Displayed genes differ from fixed selection rule.")
  assert(nrow(iterations) == nrow(input$samples) && nrow(lo) == 10 * nrow(input$samples),
         "Missing LODO iterations or selected-gene rows.")
  assert(all(table(lo$gene_id) == nrow(input$samples)) && all(table(lo$omitted_donor) == 10),
         "Incomplete crossed LODO output.")
  assert(all(iterations$tested_genes == nrow(de)), "LODO testing universe differs.")
  assert(all(abs(log(norm$TMM_factor)) < Inf) && abs(prod(norm$TMM_factor) - 1) < 1e-10,
         "Invalid centered TMM factors.")
  assert(all(norm$unfiltered_library_size == colSums(input$counts)), "Library sizes mismatch pseudobulks.")
  assert(all(ls$direction_consistent_n == vapply(ls$gene_id, function(g) sum(lo$direction_matches_full[lo$gene_id == g]), numeric(1))),
         "LODO sign summaries mismatch output.")
  assert(all(abs(ls$min_logFC - vapply(ls$gene_id, function(g) min(lo$logFC[lo$gene_id == g]), numeric(1))) < 1e-12),
         "LODO lower ranges mismatch.")
  assert(all(abs(ls$max_logFC - vapply(ls$gene_id, function(g) max(lo$logFC[lo$gene_id == g]), numeric(1))) < 1e-12),
         "LODO upper ranges mismatch.")
  assert(all(de$n_ad == sum(input$samples$diagnosis == "AD")) &&
           all(de$n_control == sum(input$samples$diagnosis == "Control")), "Exported sample counts mismatch.")
  cat("PASS: 13 empirical statistical consistency checks.\n")
}

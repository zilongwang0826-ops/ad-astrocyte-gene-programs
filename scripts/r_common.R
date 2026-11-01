# Shared statistical utilities. Execute scripts from the project root.
local_library <- normalizePath(Sys.getenv("AD_ASTRO_R_LIBRARY", "environment/R-library"), mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(jsonlite))
options(stringsAsFactors = FALSE, digits = 16)

assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "NA")
}
write_json_file <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
                       digits = NA, null = "null", na = "null")
}
require_inferential_gate <- function() {
  p <- "results/audit/dataset_audit.json"
  assert(file.exists(p), "Dataset audit JSON is missing; run the audit first.")
  audit <- jsonlite::fromJSON(p)
  assert(identical(audit$status, "PASS"),
         "BLOCKED_DATA_AUDIT: dataset audit status must equal PASS.")
}
make_design <- function(samples) {
  samples$diagnosis <- factor(samples$diagnosis, levels = c("Control", "AD"))
  assert(!anyNA(samples$diagnosis), "Unexpected or missing diagnosis.")
  design <- model.matrix(~diagnosis, data = samples)
  rownames(design) <- samples$donor_id
  assert(qr(design)$rank == ncol(design), "Diagnosis design is not identifiable.")
  assert(identical(colnames(design), c("(Intercept)", "diagnosisAD")),
         "Unexpected coefficient orientation.")
  design
}
validate_inputs <- function(counts, samples, annotation) {
  assert(is.matrix(counts) && is.numeric(counts), "Counts must be a numeric matrix.")
  assert(all(is.finite(counts)) && all(counts >= 0) && all(counts == floor(counts)),
         "Pseudobulk values must be finite nonnegative integer counts.")
  assert(!anyDuplicated(rownames(counts)) && !anyDuplicated(colnames(counts)),
         "Duplicate gene or donor identifiers in counts.")
  assert(!anyDuplicated(samples$donor_id), "Duplicated biological donor.")
  assert(identical(colnames(counts), samples$donor_id),
         "Sample metadata must exactly match count columns in order.")
  assert(identical(rownames(counts), annotation$gene_id),
         "Gene annotation must exactly match count rows in order.")
  assert(all(colSums(counts) > 0), "A donor has zero pseudobulk library size.")
  assert(all(samples$diagnosis %in% c("Control", "AD")), "Unknown diagnosis.")
  assert(all(table(factor(samples$diagnosis, levels = c("Control", "AD"))) >= 5),
         "At least five independent donors per group are required.")
  assert(all(samples$astrocyte_nuclei >= 20), "Donor below prespecified coverage screen.")
  invisible(TRUE)
}
read_inputs <- function() {
  require_inferential_gate()
  z <- read.delim(gzfile("results/tables/pseudobulk_counts.tsv.gz"), check.names = FALSE)
  assert(names(z)[1] == "gene_id", "Expected gene_id as first count column.")
  counts <- as.matrix(z[, -1, drop = FALSE]); rownames(counts) <- z$gene_id
  storage.mode(counts) <- "double"
  samples <- read.delim("results/tables/pseudobulk_samples.tsv", check.names = FALSE)
  annotation <- read.delim("results/tables/gene_annotation.tsv", check.names = FALSE)
  assert(!anyDuplicated(annotation$gene_id), "Duplicate identifiers in gene annotation.")
  annotation <- annotation[match(rownames(counts), annotation$gene_id), , drop = FALSE]
  validate_inputs(counts, samples, annotation)
  list(counts = counts, samples = samples, annotation = annotation)
}
fit_primary <- function(counts, samples, annotation, keep = NULL) {
  design <- make_design(samples)
  y <- DGEList(counts = counts, genes = annotation, group = samples$diagnosis)
  if (is.null(keep)) {
    keep <- filterByExpr(y, design = design, min.count = 10, min.total.count = 15,
                         large.n = 10, min.prop = 0.7)
  }
  assert(length(keep) == nrow(y) && sum(keep) >= 100, "Insufficient retained genes.")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM", logratioTrim = 0.3,
                      sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10)
  # In edgeR v4 the NB dispersion is estimated inside glmQLFit. Do not supply
  # gene-wise NB dispersions, or describe a separate estimateDisp step.
  fit <- glmQLFit(y, design = design, dispersion = NULL, robust = TRUE,
                  abundance.trend = TRUE, legacy = FALSE,
                  winsor.tail.p = c(0.05, 0.1), prior.count = 0.125)
  test <- glmQLFTest(fit, coef = "diagnosisAD", poisson.bound = TRUE)
  tab <- topTags(test, n = Inf, sort.by = "none", adjust.method = "BH")$table
  assert(nrow(tab) == sum(keep), "Incomplete DE result table.")
  assert(all(is.finite(as.matrix(tab[, c("logFC", "logCPM", "F", "PValue", "FDR")]))),
         "Non-finite differential-expression result.")
  assert(all(tab$F >= 0) && all(tab$PValue >= 0 & tab$PValue <= 1),
         "Invalid F statistic or raw P value.")
  assert(isTRUE(all.equal(tab$FDR, p.adjust(tab$PValue, "BH"), tolerance = 1e-14)),
         "BH correction mismatch.")
  list(y = y, fit = fit, test = test, table = tab, keep = keep, design = design)
}
unique_ranked_genes <- function(de) {
  valid <- !is.na(de$gene_symbol) & nzchar(de$gene_symbol)
  out <- de[valid, , drop = FALSE]
  # Prefer the highest expression annotation per symbol, never the best P value.
  out <- out[order(out$gene_symbol, -out$logCPM, out$gene_id, method = "radix"), ]
  out <- out[!duplicated(out$gene_symbol), ]
  out$ranking_statistic <- sign(out$logFC) * sqrt(out$F)
  assert(!anyDuplicated(out$gene_symbol), "Gene symbols must be unique for enrichment.")
  assert(all(is.finite(out$ranking_statistic)), "Invalid ranking statistic.")
  out
}
read_gmt <- function(path) {
  lines <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  assert(length(lines) > 0 && all(lengths(lines) >= 3), "Malformed gene-set GMT.")
  names_ <- vapply(lines, `[`, character(1), 1)
  assert(!anyDuplicated(names_), "Duplicate gene-set names.")
  sets <- lapply(lines, function(x) unique(x[-c(1, 2)]))
  names(sets) <- names_; sets
}
rank_effect <- function(statistic, index) {
  # Rank-biserial effect: 2*Pr(a set gene outranks a background gene)-1,
  # with average ranks for ties. This descriptive effect is not an NES.
  r <- rank(statistic, ties.method = "average")
  n <- length(r); m <- length(index)
  2 * (sum(r[index]) - m * (m + 1)/2)/(m * (n - m)) - 1
}

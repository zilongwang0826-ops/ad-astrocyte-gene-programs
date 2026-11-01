#!/usr/bin/env Rscript
source("scripts/r_common.R")
require_inferential_gate()
assert(file.exists("results/audit/analysis_summary.json"), "Run the primary DE analysis first.")
de <- read.delim("results/tables/differential_expression.tsv", check.names = FALSE)
ranked <- unique_ranked_genes(de)
statistic <- setNames(ranked$ranking_statistic, ranked$gene_symbol)
sets <- read_gmt("data/raw/h.all.v2026.1.Hs.symbols.gmt")
assert(all(startsWith(names(sets), "HALLMARK_")), "Expected only the Hallmark collection.")
indices <- lapply(sets, function(s) which(names(statistic) %in% s))
eligible <- lengths(indices) >= 15 & lengths(indices) <= 500
assert(any(eligible), "No gene sets pass mapped-size limits 15 to 500.")
set_audit <- data.frame(gene_set = names(sets), source_genes = lengths(sets),
                        mapped_genes = lengths(indices), tested = eligible)
write_tsv(set_audit, "results/tables/enrichment_set_coverage.tsv")
tab <- cameraPR(statistic, index = indices[eligible], use.ranks = TRUE,
                inter.gene.cor = 0.01, sort = FALSE, directional = TRUE)
tab$FDR <- p.adjust(tab$PValue, method = "BH")
results <- data.frame(gene_set = rownames(tab), tab, row.names = NULL,
                      rank_biserial = vapply(indices[eligible], function(ix) rank_effect(statistic, ix), numeric(1)))
assert(all(results$PValue >= 0 & results$PValue <= 1), "Invalid gene-set P value.")
assert(all(abs(results$rank_biserial) <= 1), "Invalid rank-biserial effect.")
assert(all(abs(results$rank_biserial) < 1e-14 |
           (results$Direction == "Up") == (results$rank_biserial > 0)),
       "Direction and rank-biserial effect disagree.")
results <- results[order(results$PValue, results$gene_set, method = "radix"), ]
write_tsv(results, "results/tables/enrichment_results.tsv")
write_tsv(ranked[, c("gene_id", "gene_symbol", "logFC", "logCPM", "F", "PValue", "FDR", "ranking_statistic")],
          "results/tables/enrichment_ranking.tsv")
write_json_file(list(status = "EXECUTED", method = "limma::cameraPR", collection = "MSigDB Hallmark 2026.1.Hs",
  source_file = "data/raw/h.all.v2026.1.Hs.symbols.gmt", provenance = "data/input_manifest.json",
  ranking = "sign(logFC) * sqrt(QL F statistic)", unique_symbols = length(statistic),
  tested_genes_before_symbol_mapping = nrow(de),
  genes_removed_by_symbol_mapping = nrow(de) - nrow(ranked),
  duplicate_symbol_rule = "highest primary logCPM per symbol, ties by gene_id; no P-value-based selection",
  source_gene_sets = length(sets), tested_gene_sets = nrow(results),
  significant_gene_sets = sum(results$FDR < 0.05),
  significant_up = sum(results$FDR < 0.05 & results$Direction == "Up"),
  significant_down = sum(results$FDR < 0.05 & results$Direction == "Down"),
  min_mapped_size = 15, max_mapped_size = 500, use_ranks = TRUE, directional = TRUE,
  inter_gene_correlation = 0.01, correlation_source = "prespecified fixed assumption, not estimated from donors",
  test = "two-sided competitive rank-sum test with correlation variance adjustment; analytical, no permutations",
  effect = "rank-biserial effect = 2*Pr(set gene outranks non-set gene)-1, average ranks for ties; not NES",
  positive_direction = "set genes rank toward higher expression in AD relative to other tested genes",
  multiple_testing = "Benjamini-Hochberg across all eligible gene sets", fdr_threshold = 0.05,
  interpretation = "Exploratory program-level association; significance is conditional on assumed gene correlation",
  software = list(limma = as.character(packageVersion("limma")))),
  "results/audit/enrichment_summary.json")
cat(sprintf("Enrichment complete: %d mapped unique genes, %d Hallmark sets, %d at BH FDR < 0.05.\n",
            length(statistic), nrow(results), sum(results$FDR < 0.05)))

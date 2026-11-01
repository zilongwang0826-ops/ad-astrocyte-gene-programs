#!/usr/bin/env Rscript
source("scripts/r_common.R")
require_inferential_gate()
saved <- readRDS("results/tables/primary_fit.rds")
x <- saved$input; primary <- saved$analysis
selected <- read.delim("results/tables/lodo_selected_genes.tsv", check.names = FALSE)
full <- primary$table
long <- vector("list", nrow(x$samples)); iterations <- vector("list", nrow(x$samples))
for (i in seq_len(nrow(x$samples))) {
  donor <- x$samples$donor_id[i]
  s <- x$samples[-i, , drop = FALSE]
  a <- fit_primary(x$counts[, -i, drop = FALSE], s, x$annotation, keep = primary$keep)
  z <- a$table
  assert(identical(z$gene_id, full$gene_id), "LODO gene universe changed.")
  long[[i]] <- data.frame(omitted_donor = donor, omitted_diagnosis = x$samples$diagnosis[i],
                          gene_id = z$gene_id, gene_symbol = z$gene_symbol,
                          logFC = z$logFC, PValue = z$PValue, FDR = z$FDR,
                          full_logFC = full$logFC,
                          direction_matches_full = sign(z$logFC) == sign(full$logFC),
                          absolute_effect_change = abs(z$logFC - full$logFC))
  iterations[[i]] <- data.frame(omitted_donor = donor, omitted_diagnosis = x$samples$diagnosis[i],
                                n_control = sum(s$diagnosis == "Control"), n_ad = sum(s$diagnosis == "AD"),
                                tested_genes = nrow(z), significant_genes = sum(z$FDR < 0.05),
                                NB_dispersion = as.numeric(a$fit$dispersion),
                                NB_top_proportion = a$fit$top.proportion)
  cat(sprintf("LODO %02d/%02d: removed %s (%s); %d genes retained.\n", i, nrow(x$samples), donor,
              x$samples$diagnosis[i], nrow(z)))
}
long <- do.call(rbind, long)
con <- gzfile("results/tables/lodo_all_genes.tsv.gz", "wt")
write_tsv(long, con); close(con)
iterations <- do.call(rbind, iterations)
write_tsv(iterations, "results/tables/lodo_iterations.tsv")
by_gene <- split(long, long$gene_id)
summaries <- lapply(by_gene, function(z) {
  j <- which.max(z$absolute_effect_change)
  data.frame(gene_id = z$gene_id[1], gene_symbol = z$gene_symbol[1], full_logFC = z$full_logFC[1],
              n_iterations = nrow(z), direction_consistent_n = sum(z$direction_matches_full),
              direction_consistency = mean(z$direction_matches_full),
              min_logFC = min(z$logFC), max_logFC = max(z$logFC), median_logFC = median(z$logFC),
              effect_range = diff(range(z$logFC)), max_absolute_effect_change = z$absolute_effect_change[j],
              most_influential_omission = z$omitted_donor[j])
})
all_summary <- do.call(rbind, summaries); rownames(all_summary) <- NULL
all_summary <- all_summary[match(full$gene_id, all_summary$gene_id), ]
write_tsv(all_summary, "results/tables/lodo_gene_summary.tsv")
panel <- long[long$gene_id %in% selected$gene_id, ]
panel$selection_rank <- selected$selection_rank[match(panel$gene_id, selected$gene_id)]
panel <- panel[order(panel$selection_rank, panel$omitted_donor, method = "radix"), ]
write_tsv(panel, "results/tables/lodo_selected_results.tsv")
panel_summary <- all_summary[match(selected$gene_id, all_summary$gene_id), ]
panel_summary$selection_rank <- selected$selection_rank
write_tsv(panel_summary, "results/tables/lodo_selected_summary.tsv")
write_json_file(list(status = "EXECUTED", method = "leave-one-donor-out internal sensitivity analysis",
  external_replication = FALSE, iterations = nrow(x$samples), model = "~ diagnosis",
  gene_universe = "fixed full-data filterByExpr universe", tested_genes = nrow(full),
  normalization_and_dispersion = "TMM and edgeR v4 NB/QL dispersion refitted after each omission",
  selected_rule = "10 smallest full-data raw P values, ties by gene_id; same genes in every iteration",
  selected_genes = nrow(selected), selected_all_directions_consistent = sum(panel_summary$direction_consistency == 1),
  selected_min_direction_consistency = min(panel_summary$direction_consistency),
  selected_max_absolute_change = max(panel_summary$max_absolute_effect_change),
  all_genes_all_directions_consistent = sum(all_summary$direction_consistency == 1),
  significant_full_data_genes_all_directions_consistent = sum(all_summary$direction_consistency[full$FDR < 0.05] == 1),
  selected_summary = panel_summary,
  uncertainty_note = "LODO ranges are sensitivity ranges, not confidence intervals; full-data selection favors strong signals"),
  "results/audit/lodo_summary.json")
cat(sprintf("LODO complete: %d/%d displayed genes retained their effect direction in every omission.\n",
            sum(panel_summary$direction_consistency == 1), nrow(panel_summary)))

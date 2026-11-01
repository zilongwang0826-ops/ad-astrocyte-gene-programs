#!/usr/bin/env Rscript
source("scripts/08_targeted_validation.R")
check_runtime()
checks <- 0L
check <- function(ok, label) { assert(ok, label); checks <<- checks + 1L }
rejects <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
manual_bh <- function(p) {
  o <- order(p); n <- length(p); z <- numeric(n)
  z[o] <- pmin(1, rev(cummin(rev(p[o]*n/seq_len(n))))); z
}

set.seed(174367)
s <- data.frame(donor_id = paste0("fixture_", 1:12), diagnosis = rep(c("Control", "AD"), each=6),
                astrocyte_nuclei = 100)
ann <- data.frame(gene_id = paste0("fixture_gene_",1:600), gene_symbol = paste0("FIXTURE",1:600))
mu <- matrix(150,600,12); mu[1:30,7:12] <- 600; mu[31:60,7:12] <- 40
counts <- matrix(rnbinom(7200,mu=mu,size=25),600,12,dimnames=list(ann$gene_id,s$donor_id))
a <- fit_primary(counts,s,ann)
canonical <- unique_ranked_genes(a$table)
f <- fit_primary(counts[,-1],s[-1,],ann,keep=a$keep)
stat <- fixed_ranking(f$table,canonical)
check(identical(names(stat),canonical$gene_symbol), "Canonical order lost.")
check(identical(stat,fixed_ranking(f$table[nrow(f$table):1,],canonical)), "Ranking depends on table order.")
check(rejects(fixed_ranking(f$table[-1,],canonical)), "Missing canonical gene accepted.")
bad <- f$table; bad$F[1] <- -1
check(rejects(fixed_ranking(bad,canonical)), "Negative F accepted.")
ix <- list(high=match(paste0("FIXTURE",1:30),names(stat)), low=match(paste0("FIXTURE",31:60),names(stat)), null=match(paste0("FIXTURE",101:150),names(stat)))
z <- program_test(stat,ix)
check(z$rank_biserial[1]>0 && z$rank_biserial[2]<0, "Synthetic contrast orientation wrong.")
check(isTRUE(all.equal(z$FDR,manual_bh(z$PValue),tolerance=1e-13)), "Synthetic BH family wrong.")
check(setequal(z$iteration_rank,seq_len(nrow(z))), "Set ranks are not a permutation.")
member <- data.frame(gene_set=rep(names(ix),lengths(ix)),gene_symbol=unlist(lapply(ix,function(j)names(stat)[j])))
orig <- data.frame(gene_set=names(ix),NGenes=lengths(ix))
check(identical(frozen_indices(canonical,member,orig),ix), "Frozen memberships changed.")
check(rejects(frozen_indices(canonical,rbind(member,member[1,]),orig)), "Duplicate membership accepted.")
bad <- member;bad$gene_symbol[1] <- "UNKNOWN"
check(rejects(frozen_indices(canonical,bad,orig)), "Unknown membership accepted.")
check(rejects(align_unique(rbind(s,s[1,]),"donor_id",s$donor_id)), "Duplicate donor accepted.")
check(rejects(agree_numeric(c(1,2),c(1,2.1),"fixture")), "Material baseline deviation accepted.")

donors <- data.frame(s,recorded_age_score=rep(71:76,2),sex=rep(c("F","M"),6),
                     GEO_batch=rep(1:3,4),supplementary_library_batch=rep(1:2,6),RIN=seq(6,9,length.out=12),PMI=seq(2,8,length.out=12))
diag <- design_diagnostics(donors)
check(diag$full_rank[diag$covariate=="recorded_age_score"], "Independent age design rejected.")
bad <- donors;bad$sex <- bad$diagnosis
check(!design_diagnostics(bad)$full_rank[2], "Perfect confounding not detected.")
bad$PMI[1] <- NA_real_
check(design_diagnostics(bad)$missing[6]==1, "Missing PMI not reported.")
mds <- donor_structure(a$y,donors)
lcpm <- cpm(a$y,log=TRUE,prior.count=2)[mds$genes$gene_id,]
pc <- prcomp(t(lcpm),center=TRUE,scale.=FALSE)
check(max(abs(as.vector(dist(mds$coordinates[,c("MDS1","MDS2")]))-
              as.vector(dist(pc$x[,1:2]/sqrt(nrow(lcpm))))))<1e-10,
      "Common-gene MDS disagrees with independent PCA geometry.")
bad <- donors;bad$diagnosis <- rev(bad$diagnosis)
other <- donor_structure(a$y,bad)
check(identical(mds$coordinates[,c("MDS1","MDS2")],other$coordinates[,c("MDS1","MDS2")]),
      "Diagnosis labels changed MDS coordinates.")
cat(sprintf("PASS: %d synthetic targeted-validation checks.\n",checks))

if ("--empirical" %in% commandArgs(trailingOnly=TRUE)) {
  out <- "results/targeted_validation"
  summary <- jsonlite::fromJSON(file.path(out,"validation_summary.json"))
  check(summary$status=="PASS", "Targeted run incomplete.")
  long <- read.delim(file.path(out,"gene_set_lodo.tsv"))
  sums <- read.delim(file.path(out,"gene_set_lodo_summary.tsv"))
  genes <- read.delim(gzfile(file.path(out,"lodo_gene_rankings.tsv.gz")))
  canonical <- read.delim("results/tables/enrichment_ranking.tsv")
  original <- read.delim("results/tables/enrichment_results.tsv")
  membership <- read.delim("results/sensitivity/frozen_gene_set_membership.tsv")
  donors <- read.delim("results/tables/pseudobulk_samples.tsv")
  check(nrow(long)==nrow(original)*nrow(donors) && nrow(long)==882,
        "Unexpected set/omission coverage.")
  check(!anyDuplicated(paste(long$gene_set,long$omitted_donor)), "Duplicate set/omission row.")
  for (d in donors$donor_id) {
    z <- long[long$omitted_donor==d,];g <- genes[genes$omitted_donor==d,]
    check(setequal(z$gene_set,original$gene_set), "Set family changed.")
    check(isTRUE(all.equal(z$FDR,manual_bh(z$PValue),tolerance=1e-12)), "Empirical omission BH wrong.")
    check(identical(g$gene_id,canonical$gene_id) && identical(g$gene_symbol,canonical$gene_symbol), "Canonical mapping changed.")
    check(max(abs(g$ranking_statistic-sign(g$logFC)*sqrt(g$F)))<1e-12, "Omission ranking wrong.")
    # Independent Mann-Whitney U expression for the saved rank-biserial effects.
    r <- rank(g$ranking_statistic,ties.method="average"); n <- length(r)
    effects <- vapply(z$gene_set,function(set) {
      inside <- g$gene_symbol %in% membership$gene_symbol[membership$gene_set==set]
      m <- sum(inside); U <- sum(r[inside])-m*(m+1)/2
      (2*U-m*(n-m))/(m*(n-m))
    },numeric(1))
    check(max(abs(effects-z$rank_biserial))<1e-12, "Saved set effect mismatch.")
  }
  for (g in sums$gene_set) {
    z <- long[long$gene_set==g,];s <- sums[sums$gene_set==g,]
    check(nrow(z)==18 && s$n_iterations==18 && s$direction_consistent_n==sum(sign(z$rank_biserial)==sign(s$primary_rank_biserial)) &&
      s$FDR_below_n==sum(z$FDR<.05) && abs(s$min_effect-min(z$rank_biserial))<1e-12 &&
      abs(s$max_effect-max(z$rank_biserial))<1e-12 && s$min_rank==min(z$iteration_rank) && s$max_rank==max(z$iteration_rank),
      "Omission summary differs from complete table.")
  }
  md <- read.delim(file.path(out,"donor_mds.tsv"))
  check(identical(md$donor_id,donors$donor_id) && all(is.finite(as.matrix(md[,c("MDS1","MDS2")]))), "MDS donor coverage invalid.")
  check(sum(md$age_top_coded)==7 && sum(is.na(md$PMI))==1, "Clinical metadata altered.")
  hashes <- read.delim(file.path(out,"input_code_checksums.tsv"))
  check(identical(unname(tools::sha256sum(hashes$path)),hashes$sha256), "Scientific inputs/code changed since execution.")
  cat(sprintf("PASS: %d total targeted-validation assertions (synthetic + empirical).\n",checks))
}

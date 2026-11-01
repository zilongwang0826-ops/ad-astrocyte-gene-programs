#!/usr/bin/env Rscript
# Small implementation fixtures; synthetic values never enter scientific outputs.
source("scripts/07_sensitivity.R")
assert(identical(correlation_scenarios, c(0, .01, .03, .05)), "Scenario grid changed.")
s <- data.frame(donor_id=paste0("fixture_",1:12),diagnosis=rep(c("Control","AD"),each=6),
                Age=rep(c(79,83,85,87,90,90),2))
design <- make_age_design(s)
assert(qr(design)$rank==3 && nrow(design)-ncol(design)==9, "Synthetic age design wrong.")
assert(all(design[s$Age==90,"recorded_age_centered"]==90-mean(s$Age)), "Top-coded age values were imputed.")
bad <- s; bad$Age[1] <- NA_real_
assert(inherits(try(make_age_design(bad),silent=TRUE),"try-error"),"Missing age was not rejected.")
bad <- s;bad$Age <- ifelse(bad$diagnosis=="AD",90,79)
assert(inherits(try(make_age_design(bad),silent=TRUE),"try-error"),"Perfectly collinear age was not rejected.")
bad <- s;bad$donor_id[2] <- bad$donor_id[1]
assert(inherits(try(make_age_design(bad),silent=TRUE),"try-error"),"Duplicate donors were not rejected.")
set.seed(174367)
ng <- 800L
base <- exp(rnorm(ng,log(150),.3))
mu <- matrix(rep(base,12),ng,12)
mu[1:15,7:12] <- mu[1:15,7:12]*5
mu[16:30,7:12] <- mu[16:30,7:12]/5
mu[31:100,] <- mu[31:100,]*matrix(rep(exp((s$Age-85)*.08),each=70),70,12)
counts <- matrix(rnbinom(ng*12,mu=mu,size=35),ng,12,
                  dimnames=list(paste0("fixture_gene_",1:ng),s$donor_id))
annotation <- data.frame(gene_id=rownames(counts),gene_symbol=paste0("FIXTURE",1:ng))
y <- normLibSizes(DGEList(counts=counts,genes=annotation),method="TMM")
before <- y$samples$norm.factors
fit <- fit_capped_age(y,s)
assert(all(fit$table$logFC[1:15]>0) && all(fit$table$logFC[16:30]<0),"Adjusted contrast sign wrong.")
assert(sum(fit$table$FDR[1:30]<.05)>=28,"Synthetic large signals not recovered.")
assert(identical(before,y$samples$norm.factors) && nrow(fit$table)==ng,"Normalized universe changed.")
stats <- setNames(seq(-3,3,length.out=1000),paste0("G",1:1000))
ix <- list(high=901:1000,low=1:100,mild_shift=501:600)
tabs <- lapply(correlation_scenarios,function(rho) ranked_test(stats,ix,rho))
assert(all(vapply(tabs,function(t) identical(t$Direction,tabs[[1]]$Direction),logical(1))),"Correlation altered rank direction.")
assert(all(vapply(tabs,function(t) identical(t$rank_biserial,tabs[[1]]$rank_biserial),logical(1))),"Correlation altered rank effect.")
assert(all(apply(sapply(tabs,function(t)t$PValue),1,function(v)all(diff(v)>=-1e-14))),"P values shrink with greater positive correlation.")
bh_manual <- function(p) {o<-order(p);n<-length(p);a<-pmin(1,rev(cummin(rev(p[o]*n/seq_len(n)))));z<-numeric(n);z[o]<-a;z}
assert(all(vapply(tabs,function(t)isTRUE(all.equal(t$FDR,bh_manual(t$PValue),tolerance=1e-12)),logical(1))),"Synthetic BH mismatch.")
cat("PASS: 13 synthetic sensitivity checks.\n")

if ("--empirical" %in% commandArgs(trailingOnly=TRUE)) {
  out <- "results/sensitivity"
  corr <- read.delim(file.path(out,"correlation_all_sets.tsv"),check.names=FALSE)
  old <- read.delim("results/tables/enrichment_results.tsv",check.names=FALSE)
  rank <- read.delim("results/tables/enrichment_ranking.tsv",check.names=FALSE)
  ade <- read.delim(file.path(out,"age_score_differential_expression.tsv"),check.names=FALSE)
  ae <- read.delim(file.path(out,"age_score_enrichment.tsv"),check.names=FALSE)
  ar <- read.delim(file.path(out,"age_score_ranking.tsv"),check.names=FALSE)
  original_de <- read.delim("results/tables/differential_expression.tsv",check.names=FALSE)
  pre <- jsonlite::fromJSON(file.path(out,"design_preflight.json"))
  assert(nrow(corr)==49*4 && identical(sort(unique(corr$correlation)),correlation_scenarios),"Missing correlation scenarios/sets.")
  by_set <- split(corr,corr$gene_set)
  assert(all(vapply(by_set,function(z)length(unique(z$Direction))==1,logical(1))),"Empirical direction changed across correlations.")
  assert(all(vapply(by_set,function(z)length(unique(z$rank_biserial))==1,logical(1))),"Empirical rank effect changed across correlations.")
  assert(all(vapply(by_set,function(z)all(diff(z$PValue[order(z$correlation)])>=-1e-12),logical(1))),"Empirical correlation/P monotonicity failed.")
  for (rho in correlation_scenarios) {
    z<-corr[corr$correlation==rho,]
    assert(isTRUE(all.equal(z$FDR,bh_manual(z$PValue),tolerance=1e-12)),"Empirical scenario BH mismatch.")
  }
  base<-corr[corr$correlation==.01,];base<-base[match(old$gene_set,base$gene_set),]
  assert(isTRUE(all.equal(base$PValue,old$PValue,tolerance=1e-11)) && isTRUE(all.equal(base$FDR,old$FDR,tolerance=1e-11)),"Baseline scenario mismatch.")
  assert(setequal(ade$gene_id,original_de$gene_id) && !anyDuplicated(ade$gene_id),"Adjusted gene universe changed.")
  assert(identical(ar$gene_id,rank$gene_id) && identical(ar$gene_symbol,rank$gene_symbol),"Adjusted canonical mapping changed.")
  z<-ade[match(ar$gene_id,ade$gene_id),]
  assert(isTRUE(all.equal(ar$ranking_statistic,sign(z$logFC)*sqrt(z$F),tolerance=1e-12)),"Adjusted statistic mismatch.")
  assert(isTRUE(all.equal(ade$FDR,bh_manual(ade$PValue),tolerance=1e-12)),"Adjusted gene BH mismatch.")
  assert(isTRUE(all.equal(ae$FDR,bh_manual(ae$PValue),tolerance=1e-12)) && setequal(ae$gene_set,old$gene_set),"Adjusted gene-set family/BH mismatch.")
  assert(pre$rank==3 && pre$residual_df==15 && pre$age_missing==0 && pre$predictor_VIF<2,"Design preflight mismatch.")
  donors<-read.delim(file.path(out,"donor_covariates.tsv"))
  assert(nrow(donors)==18 && sum(donors$age_top_coded)==7 && sum(is.na(donors$PMI))==1,"Clinical audit counts mismatch.")
  comp<-read.delim(file.path(out,"age_score_original_five.tsv"))
  assert(identical(comp$gene_set,old$gene_set[old$FDR<.05]) && all(comp$adjusted_FDR==ae$FDR[match(comp$gene_set,ae$gene_set)]),"Original-five comparison mismatch.")
  cat("PASS: 14 empirical sensitivity check categories, including BH checks in all four scenarios.\n")
}

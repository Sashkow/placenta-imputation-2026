#!/usr/bin/env Rscript
## Routes without imputation, 2 of 4: adds the arm "batch_in_limma": no imputation and
## no ComBat; limma on the measured cells with the dataset in the design
## (trimester + dataset + fetal sex). Same data, same stored masks, same scoring
## as routes_dissection_and_harmonizr.R. Run from the repository root after it.
##
## A gene is tested only if at least one visible dataset measures it in both
## trimesters; otherwise the trimester effect cannot be separated from the
## dataset effect and the gene is reported as not testable.
HOLDOUT_DEFAULT_OUT <- "data/pipeline/routes"
source("scripts/validation/_holdout_common.R")
HOLDOUT_DIR <- "data/pipeline/holdout/main"
ARM <- "batch_in_limma"
D <- prepare_data(config); tools <- make_de_tools(D, config)
ph <- config$phenotype; covariates <- ph$covariates %||% character(0)

de_bil <- function(X) {
  E <- X[, tools$keep_samples, drop = FALSE]
  pd_de <- D$pd[colnames(E), , drop = FALSE]
  grp <- pd_de[[ph$group_column]]; ds <- D$ds_of_sample[colnames(E)]
  both <- sapply(D$datasets, function(d)
    rowSums(!is.na(E[, ds == d & grp == ph$baseline, drop = FALSE])) > 0 &
    rowSums(!is.na(E[, ds == d & grp == ph$contrast, drop = FALSE])) > 0)
  testable <- rowSums(both) > 0
  design <- de_design(pd_de, config, informative_covariates(pd_de, covariates), batch_in_design = TRUE)
  de <- suppressWarnings(fit_limma(E[testable, , drop = FALSE], design))
  list(de = data.frame(gene = de$gene, logFC = de$logFC, adj.P.Val = de$adj.P.Val, stringsAsFactors = FALSE),
       n_not_testable = sum(!testable))
}

## ---- Part A
pg <- read.csv(file.path(OUT_DIR, "partA_per_gene.csv"), stringsAsFactors = FALSE)
r <- de_bil(D$X); i <- match(pg$gene, r$de$gene)
pg[[paste0("logFC_", ARM)]] <- r$de$logFC[i]; pg[[paste0("padj_", ARM)]] <- r$de$adj.P.Val[i]
pg[[paste0("sig_", ARM)]] <- sig_call(r$de$logFC[i], r$de$adj.P.Val[i])
write.csv(pg, file.path(OUT_DIR, "partA_per_gene.csv"), row.names = FALSE)
strata <- list(all = rep(TRUE, nrow(pg)), shared_all6 = pg$k == 6,
               partial_ref_visible = pg$k < 6 & pg$ref_visible, partial_ref_hidden = pg$k < 6 & !pg$ref_visible)
A <- read.csv(file.path(OUT_DIR, "partA_summary.csv"), stringsAsFactors = FALSE); A <- A[A$arm != ARM, ]
for (s in names(strata)) {
  p <- pg[strata[[s]], ]; l <- p[[paste0("logFC_", ARM)]]
  sa <- p[[paste0("sig_", ARM)]] %in% TRUE; sr <- p$sig_softimpute %in% TRUE
  ok <- is.finite(l) & is.finite(p$logFC_softimpute)
  A <- rbind(A, data.frame(arm = ARM, stratum = s, genes = nrow(p), tested = sum(!is.na(l)),
    deg_arm = sum(sa), deg_softimpute = sum(sr), both = sum(sa & sr), only_softimpute = sum(sr & !sa),
    only_arm = sum(sa & !sr), jaccard = sum(sa & sr) / max(1, sum(sa | sr)),
    logFC_r = cor(l[ok], p$logFC_softimpute[ok]),
    logFC_slope = unname(coef(lm(l[ok] ~ p$logFC_softimpute[ok]))[2]),
    opposite_sign_among_both = sum(sa & sr & sign(l) != sign(p$logFC_softimpute))))
}
write.csv(A, file.path(OUT_DIR, "partA_summary.csv"), row.names = FALSE)
print(A[A$arm == ARM, ], digits = 3)

## ---- Part B
B <- read.csv(file.path(OUT_DIR, "partB_by_repeat.csv"), stringsAsFactors = FALSE); B <- B[B$arm != ARM, ]
B$not_tested <- B$not_tested %||% NA
for (scheme in c("gene_dataset_block", "progressive_tax_block")) for (rep in 1:3) {
  tag <- sprintf("%s_rep%d", scheme, rep)
  mk <- readRDS(file.path(HOLDOUT_DIR, "masks", paste0(tag, ".rds")))
  vis <- mk$gene_visibility
  Xm <- D$X; Xm[mk$mask_idx] <- NA_real_
  f <- file.path(OUT_DIR, paste0("partB_calls_", tag, ".csv"))
  calls <- read.csv(f, stringsAsFactors = FALSE)
  de <- de_bil(Xm[!vis$lost_group, , drop = FALSE])$de
  i <- match(calls$gene, de$gene)
  calls[[paste0("logFC_", ARM)]] <- de$logFC[i]
  calls[[paste0("sig_", ARM)]] <- sig_call(de$logFC[i], de$adj.P.Val[i])
  write.csv(calls, f, row.names = FALSE)
  for (s in c("overall", "ref_visible", "ref_hidden")) {
    p <- calls[switch(s, overall = TRUE, ref_visible = calls$ref_visible, ref_hidden = !calls$ref_visible), ]
    l <- p[[paste0("logFC_", ARM)]]; sa <- p[[paste0("sig_", ARM)]] %in% TRUE; tr <- p$sig_full %in% TRUE
    ok <- is.finite(l) & is.finite(p$logFC_full)
    B <- rbind(B, data.frame(scheme = scheme, rep = rep, arm = ARM, stratum = s, genes = nrow(p),
      true_deg = sum(tr), recovered = sum(sa & tr), lost = sum(tr & !sa), false_new = sum(sa & !tr),
      logFC_ccc = if (sum(ok) > 2) lin_ccc(l[ok], p$logFC_full[ok]) else NA,
      logFC_slope = if (sum(ok) > 2) unname(coef(lm(l[ok] ~ p$logFC_full[ok]))[2]) else NA,
      not_tested = sum(is.na(l))))
  }
  log_msg(tag, " done")
}
write.csv(B, file.path(OUT_DIR, "partB_by_repeat.csv"), row.names = FALSE)
num <- c("genes", "true_deg", "recovered", "lost", "false_new", "logFC_ccc", "logFC_slope", "not_tested")
S <- aggregate(B[, num], B[, c("arm", "stratum", "scheme")], mean, na.action = na.pass)
S$lost_pct <- 100 * S$lost / pmax(S$true_deg, 1e-9)
S <- S[order(S$scheme, S$stratum, S$arm), ]
write.csv(S, file.path(OUT_DIR, "partB_summary.csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)

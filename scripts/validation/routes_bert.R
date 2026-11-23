#!/usr/bin/env Rscript
## Routes without imputation, 3 of 4: adds BERT arms (Batch-Effect Reduction Trees, Bioconductor package BERT,
## the successor of HarmonizR) to the route comparison. No imputation: BERT
## adjusts each feature on the batches that measure it. Same data, same stored
## masks and same scoring as routes_dissection_and_harmonizr.R.
##
## BERT takes categorical covariates only, and its "references" are bridge
## samples, not a reference batch, and cannot be combined with covariates. So
## the article's ComBat model (gestational age in weeks, linear, + fetal sex,
## reference batch) cannot be given to it. The arms are:
##   bert_trim_sex  covariates trimester + fetal sex
##   bert_trim      covariate trimester
##   bert_nocov     no covariates (comparable to the HarmonizR package)
##
## It also records why the covariate arms fail (sample counts, dataset orders),
## what the limma method of BERT returns, and the two-dataset case in which
## trimester can be given as a covariate.
##
## Run from the repository root after the other arms (BERT 1.4.0, in renv.lock):
##   Rscript scripts/validation/routes_bert.R [--no_holdout]
HOLDOUT_DEFAULT_OUT <- "data/pipeline/routes"
source("scripts/validation/_holdout_common.R")
suppressPackageStartupMessages(library(BERT))
HOLDOUT_DIR <- "data/pipeline/holdout/main"
D <- prepare_data(config); tools <- make_de_tools(D, config)
ph <- config$phenotype; covariates <- ph$covariates %||% character(0)
cat("BERT version:", as.character(packageVersion("BERT")), "\n")

bert_run <- function(X, covs) {
  pd <- D$pd[colnames(X), , drop = FALSE]
  df <- as.data.frame(t(X)); colnames(df) <- paste0("g_", rownames(X))
  df$Batch <- as.integer(factor(pd[[DATASET_ID]]))
  if ("trimester" %in% covs) df$Cov_1 <- as.integer(factor(pd[[ph$group_column]]))
  if ("sex" %in% covs)       df$Cov_2 <- as.integer(factor(pd[[covariates[1]]]))
  out <- suppressMessages(BERT(df, cores = 1, method = "ComBat", combatmode = 1,
                               qualitycontrol = FALSE, verify = TRUE))
  M <- t(as.matrix(out[, paste0("g_", rownames(X)), drop = FALSE]))
  dimnames(M) <- dimnames(X)
  list(matrix = M, lost_cells = sum(!is.na(X) & is.na(M)))
}
ARMS <- list(bert_trim_sex = c("trimester", "sex"), bert_trim = "trimester", bert_nocov = character(0))

de_na <- function(mat) {
  E <- mat[, tools$keep_samples, drop = FALSE]
  pd_de <- D$pd[colnames(E), , drop = FALSE]
  E <- E[rowSums(!is.na(E)) >= 3, , drop = FALSE]
  de <- suppressWarnings(fit_limma(E, de_design(pd_de, config, informative_covariates(pd_de, covariates))))
  data.frame(gene = de$gene, logFC = de$logFC, adj.P.Val = de$adj.P.Val, stringsAsFactors = FALSE)
}

## ---- Part A
pg <- read.csv(file.path(OUT_DIR, "partA_per_gene.csv"), stringsAsFactors = FALSE)
A <- read.csv(file.path(OUT_DIR, "partA_summary.csv"), stringsAsFactors = FALSE)
strata <- list(all = rep(TRUE, nrow(pg)), shared_all6 = pg$k == 6,
               partial_ref_visible = pg$k < 6 & pg$ref_visible, partial_ref_hidden = pg$k < 6 & !pg$ref_visible)
lost <- list(); failed <- character(0)
for (ARM in names(ARMS)) {
  t0 <- Sys.time()
  r <- tryCatch(bert_run(D$X, ARMS[[ARM]]), error = function(e) e)
  if (inherits(r, "error")) {   # e.g. "Not enough samples at batch/covariate level."
    log_msg(ARM, " FAILED on the full data: ", conditionMessage(r))
    failed <- c(failed, ARM); next
  }
  saveRDS(r$matrix, file.path(OUT_DIR, paste0("corrected_", ARM, ".rds")))
  de <- de_na(r$matrix); i <- match(pg$gene, de$gene)
  pg[[paste0("logFC_", ARM)]] <- de$logFC[i]; pg[[paste0("padj_", ARM)]] <- de$adj.P.Val[i]
  pg[[paste0("sig_", ARM)]] <- sig_call(de$logFC[i], de$adj.P.Val[i])
  lost[[ARM]] <- r$lost_cells
  log_msg(ARM, ": ", sum(pg[[paste0("sig_", ARM)]], na.rm = TRUE), " DEGs, ", sum(!is.na(de$logFC[i])),
          " genes tested, ", r$lost_cells, " measured cells set to NA (",
          round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s)")
  A <- A[A$arm != ARM, ]
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
}
ARMS <- ARMS[setdiff(names(ARMS), failed)]
write.csv(pg, file.path(OUT_DIR, "partA_per_gene.csv"), row.names = FALSE)
write.csv(A, file.path(OUT_DIR, "partA_summary.csv"), row.names = FALSE)
print(A[A$arm %in% names(ARMS), ], digits = 3, row.names = FALSE)

## BERT against our dissection route (hz_mod_ref), full data
cmp <- do.call(rbind, lapply(names(ARMS), function(ARM) do.call(rbind, lapply(names(strata), function(s) {
  p <- pg[strata[[s]], ]; x <- p$logFC_hz_mod_ref; y <- p[[paste0("logFC_", ARM)]]
  ok <- is.finite(x) & is.finite(y); sx <- p$sig_hz_mod_ref %in% TRUE; sy <- p[[paste0("sig_", ARM)]] %in% TRUE
  data.frame(arm = ARM, stratum = s, genes = sum(ok), logFC_r = cor(x[ok], y[ok]),
             logFC_slope = unname(coef(lm(y[ok] ~ x[ok]))[2]), deg_dissection = sum(sx), deg_bert = sum(sy),
             both = sum(sx & sy), jaccard = sum(sx & sy) / max(1, sum(sx | sy)))
}))))
cmp$measured_cells_set_NA <- unlist(lost)[cmp$arm]
write.csv(cmp, file.path(OUT_DIR, "bert_vs_dissection.csv"), row.names = FALSE)
cat("\nBERT against the dissection route:\n"); print(cmp, digits = 3, row.names = FALSE)

## ---- Why the covariate arms fail: sample counts, and the order of datasets in the tree
pd_all <- D$pd; ds_all <- pd_all[[DATASET_ID]]
counts <- as.data.frame(table(dataset = ds_all, trimester = pd_all[[ph$group_column]],
                              sex = pd_all[[covariates[1]]]))
write.csv(counts, file.path(OUT_DIR, "bert_sample_counts.csv"), row.names = FALSE)
bert_try <- function(order, covs) {
  df <- as.data.frame(t(D$X)); colnames(df) <- paste0("g_", rownames(D$X))
  df$Batch <- match(ds_all, order)
  if ("trimester" %in% covs) df$Cov_1 <- as.integer(factor(pd_all[[ph$group_column]]))
  if ("sex" %in% covs)       df$Cov_2 <- as.integer(factor(pd_all[[covariates[1]]]))
  tryCatch({ suppressMessages(BERT(df, cores = 1, qualitycontrol = FALSE)); "ok" },
           error = function(e) conditionMessage(e))
}
orders <- list(alphabetical = sort(unique(ds_all)),
  both_paired_with_single = c("GSE100051", "GSE37901", "GSE9984", "GSE122214", "GSE28551", "GSE93520"),
  both_first = c("GSE100051", "GSE9984", "GSE37901", "GSE122214", "GSE28551", "GSE93520"),
  both_last  = c("GSE28551", "GSE93520", "GSE122214", "GSE37901", "GSE9984", "GSE100051"))
ord <- do.call(rbind, lapply(names(orders), function(n) do.call(rbind, lapply(
  list("trimester", c("trimester", "sex")), function(cv)
    data.frame(order = n, datasets = paste(orders[[n]], collapse = " "), covariates = paste(cv, collapse = "+"),
               result = bert_try(orders[[n]], cv), stringsAsFactors = FALSE)))))
write.csv(ord, file.path(OUT_DIR, "bert_dataset_order_test.csv"), row.names = FALSE)
cat("
BERT with covariates, by order of datasets:
"); print(ord[, c("order", "covariates", "result")], row.names = FALSE)

## ---- BERT with its limma method (removeBatchEffect) instead of ComBat
## limma does not stop on a confounded covariate: the dataset effect of a
## confounded pair is not estimable and stays in the data, so it is then
## counted as a trimester effect. Recorded to show that the run "succeeds"
## with a wrong result.
full_de <- read.delim("data/pipeline/main/difexp_softimpute_combat_ref.tsv", stringsAsFactors = FALSE)
full_de$sig <- sig_call(full_de$logFC, full_de$adj.P.Val)
lim <- do.call(rbind, lapply(list("trimester", c("trimester", "sex"), character(0)), function(cv) {
  df <- as.data.frame(t(D$X)); colnames(df) <- paste0("g_", rownames(D$X))
  df$Batch <- as.integer(factor(ds_all))
  if ("trimester" %in% cv) df$Cov_1 <- as.integer(factor(pd_all[[ph$group_column]]))
  if ("sex" %in% cv)       df$Cov_2 <- as.integer(factor(pd_all[[covariates[1]]]))
  r <- tryCatch(suppressWarnings(suppressMessages(BERT(df, cores = 1, method = "limma", qualitycontrol = FALSE))),
                error = function(e) e)
  row <- data.frame(method = "limma", covariates = paste(cv, collapse = "+"), result = "ok", deg = NA, shared_with_541 = NA,
                    logFC_r = NA, logFC_slope = NA, stringsAsFactors = FALSE)
  if (inherits(r, "error")) { row$result <- conditionMessage(r); return(row) }
  Mx <- t(as.matrix(r[, paste0("g_", rownames(D$X))])); dimnames(Mx) <- dimnames(D$X)
  d <- de_na(Mx); d$sig <- sig_call(d$logFC, d$adj.P.Val)
  g <- intersect(d$gene, full_de$gene); x <- d[match(g, d$gene), ]; y <- full_de[match(g, full_de$gene), ]
  row$deg <- sum(d$sig); row$shared_with_541 <- sum(x$sig & y$sig)
  row$logFC_r <- cor(x$logFC, y$logFC, use = "complete.obs"); row$logFC_slope <- unname(coef(lm(x$logFC ~ y$logFC))[2])
  row
}))
write.csv(lim, file.path(OUT_DIR, "bert_limma_method.csv"), row.names = FALSE)
cat("
BERT with method limma:
"); print(lim, digits = 3, row.names = FALSE)

## ---- BERT on the two datasets that hold both trimesters, covariate trimester
two <- c("GSE100051", "GSE9984")
X2 <- D$X[, D$ds_of_sample %in% two]; X2 <- X2[rowSums(!is.na(X2)) > 0, ]
pd2 <- D$pd[colnames(X2), , drop = FALSE]
df2 <- as.data.frame(t(X2)); colnames(df2) <- paste0("g_", rownames(X2))
df2$Batch <- as.integer(factor(pd2[[DATASET_ID]])); df2$Cov_1 <- as.integer(factor(pd2[[ph$group_column]]))
o2 <- suppressMessages(BERT(df2, cores = 1, qualitycontrol = FALSE))
M2 <- t(as.matrix(o2[, paste0("g_", rownames(X2))])); dimnames(M2) <- dimnames(X2)
E2 <- M2[rowSums(!is.na(M2)) >= 3, ]
d2 <- suppressWarnings(fit_limma(E2, de_design(pd2, config, informative_covariates(pd2, covariates))))
d2 <- data.frame(gene = d2$gene, logFC = d2$logFC, adj.P.Val = d2$adj.P.Val, stringsAsFactors = FALSE)
d2$sig <- sig_call(d2$logFC, d2$adj.P.Val)
in_both <- rowSums(!is.na(X2[, pd2[[DATASET_ID]] == two[1]])) > 0 & rowSums(!is.na(X2[, pd2[[DATASET_ID]] == two[2]])) > 0
d2$measured_by_both <- in_both[match(d2$gene, rownames(X2))]
write.csv(d2, file.path(OUT_DIR, "bert_two_datasets_difexp.csv"), row.names = FALSE)
refs <- list(balanced_two_datasets = "data/pipeline/balanced/difexp_softimpute_combat_ref.tsv",
             full_six_datasets = "data/pipeline/main/difexp_softimpute_combat_ref.tsv")
two_cmp <- do.call(rbind, lapply(names(refs), function(n) {
  r <- read.delim(refs[[n]], stringsAsFactors = FALSE); r$sig <- sig_call(r$logFC, r$adj.P.Val)
  g <- intersect(d2$gene, r$gene); x <- d2[match(g, d2$gene), ]; y <- r[match(g, r$gene), ]
  data.frame(reference = n, genes = length(g), samples_bert = ncol(X2), deg_bert = sum(x$sig), deg_reference = sum(y$sig),
             shared = sum(x$sig & y$sig), logFC_r = cor(x$logFC, y$logFC),
             logFC_slope = unname(coef(lm(x$logFC ~ y$logFC))[2]),
             deg_bert_measured_by_both = sum(x$sig & x$measured_by_both), deg_bert_measured_by_one = sum(x$sig & !x$measured_by_both))
}))
write.csv(two_cmp, file.path(OUT_DIR, "bert_two_datasets_summary.csv"), row.names = FALSE)
cat("
BERT on the two datasets with both trimesters:
"); print(two_cmp, digits = 3, row.names = FALSE)

if (has_flag("no_holdout")) quit(save = "no")

## ---- Part B
B <- read.csv(file.path(OUT_DIR, "partB_by_repeat.csv"), stringsAsFactors = FALSE)
B <- B[!B$arm %in% names(ARMS), ]; if (is.null(B$not_tested)) B$not_tested <- NA
for (scheme in c("gene_dataset_block", "progressive_tax_block")) for (rep in 1:3) {
  tag <- sprintf("%s_rep%d", scheme, rep)
  mk <- readRDS(file.path(HOLDOUT_DIR, "masks", paste0(tag, ".rds")))
  vis <- mk$gene_visibility
  Xm <- D$X; Xm[mk$mask_idx] <- NA_real_
  f <- file.path(OUT_DIR, paste0("partB_calls_", tag, ".csv"))
  calls <- read.csv(f, stringsAsFactors = FALSE)
  for (ARM in names(ARMS)) {
    r <- tryCatch(bert_run(Xm[!vis$lost_group, , drop = FALSE], ARMS[[ARM]]), error = function(e) e)
    if (inherits(r, "error")) { log_msg(tag, " ", ARM, " FAILED: ", conditionMessage(r)); next }
    de <- de_na(r$matrix); i <- match(calls$gene, de$gene)
    calls[[paste0("logFC_", ARM)]] <- de$logFC[i]
    calls[[paste0("sig_", ARM)]] <- sig_call(de$logFC[i], de$adj.P.Val[i])
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
  }
  write.csv(calls, f, row.names = FALSE)
  log_msg(tag, " done")
}
write.csv(B, file.path(OUT_DIR, "partB_by_repeat.csv"), row.names = FALSE)
num <- c("genes", "true_deg", "recovered", "lost", "false_new", "logFC_ccc", "logFC_slope", "not_tested")
S <- aggregate(B[, num], B[, c("arm", "stratum", "scheme")], mean, na.action = na.pass)
S$lost_pct <- 100 * S$lost / pmax(S$true_deg, 1e-9)
S <- S[order(S$scheme, S$stratum, S$arm), ]
write.csv(S, file.path(OUT_DIR, "partB_summary.csv"), row.names = FALSE)
print(S[S$arm %in% c(names(ARMS), "hz_mod_ref", "softimpute"), ], digits = 3, row.names = FALSE)

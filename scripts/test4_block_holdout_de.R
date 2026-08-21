#!/usr/bin/env Rscript
#
# Test 4: block-holdout differential-expression validation.
#
# THE QUESTION. The pipeline reports 262 "gained" DEGs -- genes it can only
# test because imputation filled the datasets that never measured them. Nothing
# in the existing validation checks whether the DE answer for such a gene is
# right, because for a real gained gene there is no ground truth to check
# against. This test manufactures ground truth: it takes genes measured on all
# six platforms (whose correct answer IS known, from the full run), hides their
# values in the same coverage patterns real gained genes have, and asks whether
# the pipeline still reaches the known answer.
#
# THREE ARMS per masked gene. Each is a limma fit; they differ in exactly two
# ingredients, which is what lets the total error be split into its two causes:
#
#                        ComBat parameters estimated on...
#                        the full matrix        the masked matrix
#                        (never saw this        (this gene's hidden
#                         gene imputed)          blocks were imputed)
#   limma uses...
#   all 117 samples      FULL  (ground truth)   --
#   only the samples     ORACLE                 MASKED (= production)
#   left visible
#
#   ORACLE vs FULL  -- same batch correction, fewer real samples.
#                      Isolates loss of power / coverage.
#   MASKED vs ORACLE -- same samples and weights, different batch correction.
#                      Isolates distortion leaking in through ComBat, which
#                      with weight-0 imputed cells is the ONLY route by which
#                      an imputed value can reach a p-value.
#   MASKED vs FULL  -- both effects together: what a reader of the published
#                      DEG list actually experiences for a gained gene.
#
# ORACLE is built by running limma on the FULL run's ComBat-corrected matrix
# while applying the MASKED run's weights. The corrected matrix embodies the
# full-data ComBat parameters, so no parameter surgery is needed -- the D2
# fallback in design.md is not required.
#
# Usage:
#   Rscript scripts/test4_block_holdout_de.R
#   Rscript scripts/test4_block_holdout_de.R --n_genes=200 --seeds=4001

suppressPackageStartupMessages({
  library(limma)
})

script_dir <- {
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "scripts"
}
source(file.path(script_dir, "pipeline", "subsampling_helpers.R"))

args <- commandArgs(TRUE)
get_arg <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", n, "="), "", h[1]) else d }

config <- parse_config_arg("config/config_test4_block_holdout.yaml")
tcfg <- config$test4

n_genes <- as.integer(get_arg("n_genes", tcfg$n_genes_per_replicate %||% 800))
seeds <- if (!is.null(get_arg("seeds"))) {
  as.integer(strsplit(get_arg("seeds"), ",")[[1]])
} else as.integer(unlist(tcfg$seeds))
REF_BATCH <- tcfg$ref_batch %||% "GSE100051"

out_dir <- config$paths$output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR <- config$thresholds$fdr %||% 0.05
LFC <- config$thresholds$logfc %||% 1.0

cat("=== Test 4: block-holdout DE validation ===\n")
cat("Genes per replicate:", n_genes, " Replicates:", length(seeds), "\n\n")

## ---------------------------------------------------------------- load data
data <- load_pipeline_data(config)
exprs_list <- data$exprs_list
pdata <- data$phenodata
dsnames <- names(exprs_list)
n_ds <- length(dsnames)

incomplete <- create_incomplete_matrix(exprs_list, min_datasets = 1L)
X_full <- incomplete$matrix
orig_mask <- is.na(X_full)
genes <- rownames(X_full); samples <- colnames(X_full)
cat("Merged matrix:", nrow(X_full), "genes x", ncol(X_full), "samples;",
    sprintf("%.1f%% missing\n", 100 * mean(orig_mask)))

## which datasets measured each gene (computed after the drop above)
present <- sapply(exprs_list, function(m) genes %in% rownames(m))
rownames(present) <- genes
k <- rowSums(present)
intersection_genes <- genes[k == n_ds]
cat("Intersection genes available to mask:", length(intersection_genes), "\n")

## sample-to-dataset map
ds_of_sample <- setNames(as.character(pdata$secondaryaccession),
                         pdata$arraydatafile_exprscolumnnames)[samples]

## Match production: the pipeline drops genes with no real measurement in one
## of the two groups, because limma cannot estimate a fold change for them.
## Keeping them would let 298 unfittable genes into the eBayes variance
## moderation that the primary analysis never sees.
grp_of_sample <- setNames(as.character(pdata[[config$phenotype$group_column]]),
                          pdata$arraydatafile_exprscolumnnames)[samples]
real_base <- rowSums(!orig_mask[, grp_of_sample == config$phenotype$baseline, drop = FALSE])
real_cont <- rowSums(!orig_mask[, grp_of_sample == config$phenotype$contrast, drop = FALSE])
drop_genes <- genes[real_base == 0 | real_cont == 0]
if (length(drop_genes)) {
  cat("Dropping", length(drop_genes), "genes with no real values in one group",
      "(as the primary pipeline does)\n")
  keep_g <- !(genes %in% drop_genes)
  X_full <- X_full[keep_g, , drop = FALSE]
  orig_mask <- orig_mask[keep_g, , drop = FALSE]
  incomplete$matrix <- X_full
  if (!is.null(incomplete$gene_info)) incomplete$gene_info <- incomplete$gene_info[keep_g, , drop = FALSE]
  genes <- rownames(X_full)
}

## ------------------------------------------- empirical patterns to mask into
non_int <- present[k < n_ds, , drop = FALSE]
pat_str <- apply(non_int, 1, function(r) paste(dsnames[r], collapse = "+"))
pat_tab <- sort(table(pat_str), decreasing = TRUE)
pat_names <- names(pat_tab)
pat_prob <- as.numeric(pat_tab) / sum(pat_tab)
cat("Empirical presence patterns among non-intersection genes:",
    length(pat_names), "\n")
for (i in seq_along(pat_names))
  cat(sprintf("  %-58s %5d (%.1f%%)\n", pat_names[i], pat_tab[i], 100 * pat_prob[i]))

## ------------------------------------------------------- shared limma runner
group_col <- config$phenotype$group_column
baseline  <- config$phenotype$baseline
contrast  <- config$phenotype$contrast
covariates <- config$phenotype$covariates %||% character(0)

pd <- pdata[match(samples, pdata$arraydatafile_exprscolumnnames), ]
rownames(pd) <- samples
keep_samples <- pd[[group_col]] %in% c(baseline, contrast)

build_design <- function(pdsub) {
  grp <- factor(pdsub[[group_col]], levels = c(baseline, contrast))
  dd <- data.frame(group = grp, stringsAsFactors = FALSE)
  terms <- "group"
  for (cov in covariates) {
    if (length(unique(na.omit(pdsub[[cov]]))) >= 2) {
      dd[[cov]] <- pdsub[[cov]]; terms <- c(terms, cov)
    }
  }
  ok <- stats::complete.cases(dd)
  list(design = model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))),
                             data = dd[ok, , drop = FALSE]),
       keep = ok)
}

de_fit <- function(corrected, weight_mask) {
  ## corrected: gene x sample corrected expression
  ## weight_mask: TRUE where the cell must get weight 0
  E <- corrected[, keep_samples, drop = FALSE]
  W_mask <- weight_mask[rownames(E), colnames(E), drop = FALSE]
  bd <- build_design(pd[keep_samples, , drop = FALSE])
  E <- E[, bd$keep, drop = FALSE]
  W_mask <- W_mask[, bd$keep, drop = FALSE]
  W <- matrix(1.0, nrow(E), ncol(E), dimnames = dimnames(E))
  W[W_mask] <- 0.0
  ## limma cannot fit a gene with no positive-weight observation in a group;
  ## such genes are exactly the ones the pipeline drops, and are recorded as
  ## dropped rather than silently returning NA.
  fit <- eBayes(lmFit(E, bd$design, weights = W))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC,
             adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
}

combat_correct <- function(mat) {
  merged_pdata <- pd[colnames(mat), , drop = FALSE]
  batch <- as.factor(merged_pdata$secondaryaccession)
  cbc <- config$phenotype$combat_bio_covariate
  if (!is.null(cbc)) {
    if (!cbc %in% colnames(merged_pdata))
      merged_pdata[[cbc]] <- compute_ga_weeks(merged_pdata)
    md <- data.frame(ga_linear = merged_pdata[[cbc]])
    terms <- "ga_linear"
  } else {
    md <- data.frame(bio_group = merged_pdata[[group_col]], stringsAsFactors = FALSE)
    terms <- "bio_group"
  }
  for (cov in covariates)
    if (length(unique(na.omit(merged_pdata[[cov]]))) >= 2) {
      md[[cov]] <- merged_pdata[[cov]]; terms <- c(terms, cov)
    }
  mod <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = md)
  rb <- config$normalization$combat_ref$ref_batch %||% REF_BATCH
  tryCatch(normalize_combat_ref(mat, batch, mod = mod, ref_batch = rb),
           error = function(e) normalize_combat(mat, batch, mod = mod))
}

## --------------------------------------------------------------- FULL arm
cat("\n--- FULL arm (ground truth) ---\n")
ref_corr_path <- config$paths$reference_corrected
if (!is.null(ref_corr_path) && file.exists(ref_corr_path)) {
  cat("Using shipped corrected matrix:", ref_corr_path, "\n")
  corrected_full <- as.matrix(read.delim(ref_corr_path, row.names = 1, check.names = FALSE))
  corrected_full <- corrected_full[intersect(genes, rownames(corrected_full)),
                                   intersect(samples, colnames(corrected_full)), drop = FALSE]
} else {
  imp_full <- run_imputer("softimpute", incomplete, config$imputation$softimpute)
  corrected_full <- combat_correct(imp_full$matrix)
}
de_full <- de_fit(corrected_full, orig_mask)
cat("FULL DEGs:", sum(de_full$adj.P.Val < FDR & abs(de_full$logFC) > LFC, na.rm = TRUE), "\n")

## ------------------------------------------------------------- replicates
all_rows <- list()
used <- character(0)

for (ri in seq_along(seeds)) {
  seed <- seeds[ri]
  cat("\n--- Replicate", ri, "(seed", seed, ") ---\n")
  set.seed(seed)

  pool <- setdiff(intersection_genes, used)
  if (length(pool) < n_genes) stop("not enough unused intersection genes")
  picked <- sort(sample(pool, n_genes))
  used <- c(used, picked)

  assigned <- sample(pat_names, n_genes, replace = TRUE, prob = pat_prob)
  names(assigned) <- picked

  ## build the mask: hide every cell of a picked gene in datasets outside its
  ## assigned pattern
  new_mask <- matrix(FALSE, nrow(X_full), ncol(X_full), dimnames = dimnames(X_full))
  for (g in picked) {
    keep_ds <- strsplit(assigned[[g]], "\\+")[[1]]
    new_mask[g, !(ds_of_sample %in% keep_ds)] <- TRUE
  }
  mask_all <- orig_mask | new_mask
  cat("Added missing cells:", sum(new_mask),
      sprintf(" (total missing now %.1f%%)\n", 100 * mean(mask_all)))

  ## MASKED arm: impute the masked matrix, ComBat on it, limma with w=0
  X_masked <- X_full
  X_masked[new_mask] <- NA_real_
  inc_masked <- incomplete
  inc_masked$matrix <- X_masked
  imp_masked <- run_imputer("softimpute", inc_masked, config$imputation$softimpute)
  corrected_masked <- combat_correct(imp_masked$matrix)
  de_masked <- de_fit(corrected_masked, mask_all)

  ## ORACLE arm: FULL's corrected matrix, MASKED's weights
  de_oracle <- de_fit(corrected_full, mask_all)

  ## ------------------------------------------------------------ assemble
  idx <- function(d) match(picked, d$gene)
  fi <- idx(de_full); mi <- idx(de_masked); oi <- idx(de_oracle)

  kept_ds <- sapply(picked, function(g) length(strsplit(assigned[[g]], "\\+")[[1]]))
  has_ref <- sapply(picked, function(g) REF_BATCH %in% strsplit(assigned[[g]], "\\+")[[1]])

  sig <- function(l, p) !is.na(p) & !is.na(l) & p < FDR & abs(l) > LFC

  rows <- data.frame(
    replicate = ri, seed = seed, gene = picked,
    pattern = unname(assigned[picked]), k_visible = unname(kept_ds),
    ref_batch_visible = unname(has_ref),
    logFC_full   = de_full$logFC[fi],   padj_full   = de_full$adj.P.Val[fi],
    logFC_masked = de_masked$logFC[mi], padj_masked = de_masked$adj.P.Val[mi],
    logFC_oracle = de_oracle$logFC[oi], padj_oracle = de_oracle$adj.P.Val[oi],
    stringsAsFactors = FALSE
  )
  rows$sig_full   <- sig(rows$logFC_full,   rows$padj_full)
  rows$sig_masked <- sig(rows$logFC_masked, rows$padj_masked)
  rows$sig_oracle <- sig(rows$logFC_oracle, rows$padj_oracle)
  rows$dropped_masked <- is.na(rows$logFC_masked)

  cat(sprintf("  FULL sig: %d   ORACLE sig: %d   MASKED sig: %d   dropped: %d\n",
              sum(rows$sig_full, na.rm = TRUE), sum(rows$sig_oracle, na.rm = TRUE),
              sum(rows$sig_masked, na.rm = TRUE), sum(rows$dropped_masked)))
  all_rows[[ri]] <- rows
}

res <- do.call(rbind, all_rows)
write.csv(res, file.path(out_dir, "test4_per_gene.csv"), row.names = FALSE)

## ----------------------------------------------------------------- summarise
summarise <- function(d, label) {
  ok <- !is.na(d$logFC_masked) & !is.na(d$logFC_full)
  okO <- !is.na(d$logFC_oracle) & !is.na(d$logFC_full)
  data.frame(
    stratum = label, n_genes = nrow(d),
    dropped = sum(d$dropped_masked),
    logfc_r_MF   = if (sum(ok) > 2) cor(d$logFC_masked[ok], d$logFC_full[ok]) else NA,
    logfc_mae_MF = if (any(ok)) mean(abs(d$logFC_masked[ok] - d$logFC_full[ok])) else NA,
    logfc_r_OF   = if (sum(okO) > 2) cor(d$logFC_oracle[okO], d$logFC_full[okO]) else NA,
    logfc_mae_OF = if (any(okO)) mean(abs(d$logFC_oracle[okO] - d$logFC_full[okO])) else NA,
    logfc_mae_MO = if (any(ok)) mean(abs(d$logFC_masked[ok] - d$logFC_oracle[ok]), na.rm = TRUE) else NA,
    n_sig_full = sum(d$sig_full, na.rm = TRUE),
    sensitivity = { p <- d$sig_full %in% TRUE
      if (sum(p)) mean(d$sig_masked[p] %in% TRUE) else NA },
    specificity = { n <- d$sig_full %in% FALSE
      if (sum(n)) mean(!(d$sig_masked[n] %in% TRUE)) else NA },
    direction_flip = { p <- d$sig_full %in% TRUE & !is.na(d$logFC_masked)
      if (sum(p)) mean(sign(d$logFC_masked[p]) != sign(d$logFC_full[p])) else NA },
    stringsAsFactors = FALSE
  )
}

strata <- list(summarise(res, "ALL"))
for (kk in sort(unique(res$k_visible)))
  strata[[length(strata) + 1]] <- summarise(res[res$k_visible == kk, ], paste0("k=", kk))
for (rb in c(TRUE, FALSE)) {
  d <- res[res$ref_batch_visible == rb, ]
  if (nrow(d)) strata[[length(strata) + 1]] <-
    summarise(d, paste0("ref_batch_", ifelse(rb, "visible", "hidden")))
}
summ <- do.call(rbind, strata)
write.csv(summ, file.path(out_dir, "test4_summary.csv"), row.names = FALSE)

cat("\n=== Test 4 summary ===\n")
print(summ, row.names = FALSE, digits = 3)

## attribution headline
ok <- !is.na(res$logFC_masked) & !is.na(res$logFC_oracle) & !is.na(res$logFC_full)
tot <- mean(abs(res$logFC_masked[ok] - res$logFC_full[ok]))
cov_loss <- mean(abs(res$logFC_oracle[ok] - res$logFC_full[ok]))
combat_path <- mean(abs(res$logFC_masked[ok] - res$logFC_oracle[ok]))
cat(sprintf("\nMean |logFC| discrepancy by comparison (n=%d genes):\n", sum(ok)))
cat(sprintf("  MASKED vs FULL   (both effects)       : %.4f\n", tot))
cat(sprintf("  ORACLE vs FULL   (coverage/power loss): %.4f\n", cov_loss))
cat(sprintf("  MASKED vs ORACLE (ComBat pathway)     : %.4f\n", combat_path))
cat(sprintf("  ratio ComBat-pathway : coverage-loss  : %.2f\n", combat_path / cov_loss))
cat("  (these are mean absolute errors of three different differences and do\n",
    "   not sum; the ratio is the interpretable quantity)\n", sep = "")

cat("\nWritten to:", out_dir, "\n")

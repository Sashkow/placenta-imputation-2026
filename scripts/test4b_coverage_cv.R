#!/usr/bin/env Rscript
#
# SUPERSEDED (2026-08): this script belongs to the earlier empirical-pattern /
# uniform-k holdout design and is no longer run by run_everything.R. Its role is
# taken over by holdout_unified.R (per-gene values by GSE100051 presence; holdout_tables.R Table S2), which use one set of masks for both the imputation-
# accuracy and the DE-call readouts. Kept for the record.
#
#
# Test 4b: coverage-stratified imputation cross-validation.
#
# The existing block-mask CV (Supp Table S1) hides ONE dataset's block from
# genes measured in all six, then reports a single Pearson r pooled over every
# masked cell. Two things are missing:
#
#   1. Real imputation-only genes are usually left with far less than five
#      datasets. Accuracy at k=5 says little about a gene that survives in one.
#      Here genes are masked down to a TARGET TIER k = 1, 2, 3, 4, 5.
#   2. A cell-pooled correlation is dominated by genes with many masked cells
#      and by between-gene spread in expression level. A gene-level median
#      answers the question a reader actually has -- "how well is a TYPICAL
#      gene reconstructed?" -- so per-gene r and MAE are reported alongside.
#
# Usage:
#   Rscript scripts/test4b_coverage_cv.R
#   Rscript scripts/test4b_coverage_cv.R --n_genes=300 --tiers=1,3,5 --n_repeats=1

script_dir <- {
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "scripts"
}
source(file.path(script_dir, "pipeline", "subsampling_helpers.R"))

args <- commandArgs(TRUE)
get_arg <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", n, "="), "", h[1]) else d }

config <- parse_config_arg("config/config_test4_block_holdout.yaml")
out_dir <- config$paths$output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_genes   <- as.integer(get_arg("n_genes", 800))
tiers     <- as.integer(strsplit(get_arg("tiers", "1,2,3,4,5"), ",")[[1]])
n_repeats <- as.integer(get_arg("n_repeats", 3))

cat("=== Test 4b: coverage-stratified imputation CV ===\n")
cat("Genes per tier per repeat:", n_genes, " Tiers:", paste(tiers, collapse = ","),
    " Repeats:", n_repeats, "\n\n")

data <- load_pipeline_data(config)
exprs_list <- data$exprs_list
pdata <- data$phenodata
dsnames <- names(exprs_list)
n_ds <- length(dsnames)

incomplete <- create_incomplete_matrix(exprs_list, min_datasets = 1L)
X <- incomplete$matrix
orig_mask <- is.na(X)
genes <- rownames(X); samples <- colnames(X)

present <- sapply(exprs_list, function(m) genes %in% rownames(m))
rownames(present) <- genes
k_gene <- rowSums(present)
inter <- genes[k_gene == n_ds]
ds_of_sample <- setNames(as.character(pdata$secondaryaccession),
                         pdata$arraydatafile_exprscolumnnames)[samples]
cat("Intersection genes:", length(inter), "\n\n")

rows <- list()
per_gene_rows <- list()

for (tier in tiers) {
  for (rep in seq_len(n_repeats)) {
    set.seed(1000 * tier + rep)
    picked <- sample(inter, min(n_genes, length(inter)))

    ## each picked gene keeps `tier` randomly chosen datasets, loses the rest
    mask <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
    for (g in picked) {
      keep_ds <- sample(dsnames, tier)
      mask[g, !(ds_of_sample %in% keep_ds)] <- TRUE
    }
    n_masked <- sum(mask)

    Xm <- X; Xm[mask] <- NA_real_
    inc <- incomplete; inc$matrix <- Xm
    imp <- run_imputer("softimpute", inc, config$imputation$softimpute)

    truth <- X[mask]; filled <- imp$matrix[mask]
    pooled_r   <- cor(truth, filled, use = "complete.obs")
    pooled_rmse <- sqrt(mean((truth - filled)^2, na.rm = TRUE))
    pooled_mae  <- mean(abs(truth - filled), na.rm = TRUE)

    ## per-gene metrics
    pg <- do.call(rbind, lapply(picked, function(g) {
      idx <- mask[g, ]
      if (sum(idx) < 3) return(NULL)
      t <- X[g, idx]; f <- imp$matrix[g, idx]
      data.frame(tier = tier, repeat_n = rep, gene = g, n_cells = sum(idx),
                 r = suppressWarnings(cor(t, f)),
                 mae = mean(abs(t - f)),
                 rmse = sqrt(mean((t - f)^2)), stringsAsFactors = FALSE)
    }))
    per_gene_rows[[length(per_gene_rows) + 1]] <- pg

    rows[[length(rows) + 1]] <- data.frame(
      tier = tier, repeat_n = rep, n_genes = length(picked), n_cells_masked = n_masked,
      pooled_r = pooled_r, pooled_rmse = pooled_rmse, pooled_mae = pooled_mae,
      median_gene_r = median(pg$r, na.rm = TRUE),
      median_gene_mae = median(pg$mae, na.rm = TRUE),
      median_gene_rmse = median(pg$rmse, na.rm = TRUE),
      stringsAsFactors = FALSE)

    cat(sprintf("tier k=%d rep %d: pooled r=%.3f MAE=%.3f | median per-gene r=%.3f MAE=%.3f\n",
                tier, rep, pooled_r, pooled_mae,
                median(pg$r, na.rm = TRUE), median(pg$mae, na.rm = TRUE)))
  }
}

res <- do.call(rbind, rows)
per_gene <- do.call(rbind, per_gene_rows)
write.csv(res, file.path(out_dir, "test4b_coverage_cv.csv"), row.names = FALSE)
write.csv(per_gene, file.path(out_dir, "test4b_coverage_cv_per_gene.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, res$tier), function(d) data.frame(
  tier = d$tier[1],
  pooled_r        = sprintf("%.3f +/- %.3f", mean(d$pooled_r), sd(d$pooled_r)),
  pooled_mae      = sprintf("%.3f +/- %.3f", mean(d$pooled_mae), sd(d$pooled_mae)),
  median_gene_r   = sprintf("%.3f +/- %.3f", mean(d$median_gene_r), sd(d$median_gene_r)),
  median_gene_mae = sprintf("%.3f +/- %.3f", mean(d$median_gene_mae), sd(d$median_gene_mae)),
  stringsAsFactors = FALSE)))

cat("\n=== Accuracy by retained-coverage tier ===\n")
print(summ, row.names = FALSE)
write.csv(summ, file.path(out_dir, "test4b_summary.csv"), row.names = FALSE)
cat("\nWritten to:", out_dir, "\n")

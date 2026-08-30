#!/usr/bin/env Rscript
#
# SUPERSEDED (2026-08): this script belongs to the earlier empirical-pattern /
# uniform-k holdout design and is no longer run by run_everything.R. Its role is
# taken over by holdout_tables.R (Table S2), which use one set of masks for both the imputation-
# accuracy and the DE-call readouts. Kept for the record.
#
#
# Test 4b, post-hoc: split the coverage-stratified CV results by whether the
# ComBat reference batch was among the datasets left visible.
#
# The question this answers is the one that decides how test 4's result is
# read. Test 4 finds that a masked gene is recovered as differentially
# expressed almost always when the reference batch survives and never when it
# does not. Two explanations fit: either genes missing the reference batch are
# imputed worse, or the imputed values are fine and the damage is done later,
# by reference-batch ComBat. Comparing reconstruction error across the same
# split separates them.
#
# test4b_coverage_cv.R does not record which datasets it kept, but the choice
# is deterministic given the seed, so it is reconstructed here rather than
# paying for another 15 softImpute fits. The reconstruction is checked against
# the masked-cell counts the original run recorded; if the RNG stream had
# drifted, those counts would not agree.
#
# Usage:
#   Rscript scripts/test4b_refbatch_split.R

script_dir <- {
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "scripts"
}
source(file.path(script_dir, "pipeline", "subsampling_helpers.R"))

config <- parse_config_arg("config/config_test4_block_holdout.yaml")
out_dir <- config$paths$output
ref_batch <- config$test4$ref_batch
if (is.null(ref_batch)) ref_batch <- config$normalization$combat_ref$ref_batch
if (is.null(ref_batch)) stop("config does not name a ComBat reference batch")

data <- load_pipeline_data(config)
exprs_list <- data$exprs_list
pdata <- data$phenodata
dsnames <- names(exprs_list)
n_ds <- length(dsnames)

incomplete <- create_incomplete_matrix(exprs_list, min_datasets = 1L)
X <- incomplete$matrix
genes <- rownames(X); samples <- colnames(X)
present <- sapply(exprs_list, function(m) genes %in% rownames(m))
inter <- genes[rowSums(present) == n_ds]
ds_of_sample <- setNames(as.character(pdata$secondaryaccession),
                         pdata$arraydatafile_exprscolumnnames)[samples]

cv <- read.csv(file.path(out_dir, "test4b_coverage_cv.csv"), stringsAsFactors = FALSE)
pg <- read.csv(file.path(out_dir, "test4b_coverage_cv_per_gene.csv"), stringsAsFactors = FALSE)

## Replay the mask assignment. Draw order must match test4b_coverage_cv.R
## exactly: one sample() for the gene set, then one per gene for its datasets.
rows <- list()
for (tier in sort(unique(cv$tier))) {
  for (rep in sort(unique(cv$repeat_n[cv$tier == tier]))) {
    n_genes <- cv$n_genes[cv$tier == tier & cv$repeat_n == rep]
    set.seed(1000 * tier + rep)
    picked <- sample(inter, min(n_genes, length(inter)))
    ref_vis <- logical(length(picked))
    n_cells <- 0L
    for (i in seq_along(picked)) {
      keep_ds <- sample(dsnames, tier)
      ref_vis[i] <- ref_batch %in% keep_ds
      n_cells <- n_cells + sum(!(ds_of_sample %in% keep_ds))
    }
    recorded <- cv$n_cells_masked[cv$tier == tier & cv$repeat_n == rep]
    if (n_cells != recorded) {
      stop(sprintf("mask replay mismatch at tier %d repeat %d: %d cells vs %d recorded",
                   tier, rep, n_cells, recorded))
    }
    rows[[length(rows) + 1]] <- data.frame(
      tier = tier, repeat_n = rep, gene = picked,
      ref_visible = ref_vis, stringsAsFactors = FALSE)
  }
}
key <- do.call(rbind, rows)
cat("Mask replay verified against recorded cell counts for all",
    nrow(cv), "tier/repeat combinations.\n\n")

m <- merge(pg, key, by = c("tier", "repeat_n", "gene"))
if (nrow(m) != nrow(pg)) {
  stop(sprintf("join lost rows: %d of %d matched", nrow(m), nrow(pg)))
}

by_tier <- do.call(rbind, lapply(split(m, list(m$tier, m$ref_visible), drop = TRUE),
  function(d) data.frame(tier = d$tier[1], ref_visible = d$ref_visible[1],
                         n_genes = nrow(d),
                         median_mae = median(d$mae, na.rm = TRUE),
                         median_r = median(d$r, na.rm = TRUE),
                         stringsAsFactors = FALSE)))
by_tier <- by_tier[order(by_tier$tier, by_tier$ref_visible), ]

wt <- wilcox.test(mae ~ ref_visible, data = m)
overall <- data.frame(
  tier = "ALL",
  ref_visible = c(FALSE, TRUE),
  n_genes = c(sum(!m$ref_visible), sum(m$ref_visible)),
  median_mae = c(median(m$mae[!m$ref_visible]), median(m$mae[m$ref_visible])),
  median_r = c(median(m$r[!m$ref_visible], na.rm = TRUE),
               median(m$r[m$ref_visible], na.rm = TRUE)),
  stringsAsFactors = FALSE)

res <- rbind(transform(by_tier, tier = as.character(tier)), overall)
res$wilcox_p <- c(rep(NA_real_, nrow(by_tier)), rep(wt$p.value, 2))
write.csv(res, file.path(out_dir, "test4b_refbatch_split.csv"), row.names = FALSE)

cat("=== Imputation accuracy by reference-batch visibility ===\n")
print(res, row.names = FALSE)
cat(sprintf("\nWilcoxon rank-sum on per-gene MAE (visible vs hidden): p = %.3g\n",
            wt$p.value))
cat("\nWritten to:", file.path(out_dir, "test4b_refbatch_split.csv"), "\n")

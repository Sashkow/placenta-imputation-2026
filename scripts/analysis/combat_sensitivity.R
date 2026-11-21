#!/usr/bin/env Rscript
#
# ComBat sensitivity analysis: compare batch parameters (gamma, delta)
# between intersection-only (8,260 genes) and softImpute (17,531 genes) matrices.
#
# gamma.hat = additive batch effect per gene (batch mean shift)
# delta.hat = multiplicative batch effect per gene (batch variance scaling)
#
# If imputation doesn't distort ComBat, the parameters for the ~8,260 shared
# genes should be nearly identical between the two runs.

library(sva)

base_dir <- "data/pipeline/sensitivity"
pheno_path <- "data/phenodata.tsv"
output_dir <- "data/pipeline/sensitivity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pheno <- read.delim(pheno_path, stringsAsFactors = FALSE)
datasets_6ds <- c("GSE100051", "GSE122214", "GSE28551", "GSE37901", "GSE93520", "GSE9984")
pheno <- pheno[pheno$dataset_id %in% datasets_6ds &
               pheno$condition == "healthy" &
               pheno$trimester %in% c("First trimester", "Second trimester"), ]

# -------------------------------------------------------------------
# Extract ComBat batch parameters (gamma.hat, delta.hat) from a matrix
# by replicating the standardization step from sva::ComBat
# -------------------------------------------------------------------
extract_combat_params <- function(exprs, batch, mod = NULL, ref_batch = NULL) {
  batch <- as.factor(batch)
  if (!is.null(ref_batch)) {
    batch <- relevel(batch, ref = ref_batch)
  }

  batchmod <- model.matrix(~ -1 + batch)
  n_batch <- nlevels(batch)
  batches <- lapply(levels(batch), function(b) which(batch == b))
  n_batches <- sapply(batches, length)

  if (is.null(mod)) {
    design <- batchmod
  } else {
    design <- cbind(batchmod, mod[, -1, drop = FALSE])
  }

  # Standardize: fit gene-wise OLS, get residuals
  B.hat <- solve(crossprod(design), t(design) %*% t(exprs))
  grand.mean <- crossprod(n_batches / ncol(exprs), B.hat[1:n_batch, ])
  var.pooled <- ((exprs - t(design %*% B.hat))^2) %*% rep(1 / ncol(exprs), ncol(exprs))

  stand.mean <- crossprod(grand.mean, t(rep(1, ncol(exprs))))
  s.data <- (exprs - stand.mean) / sqrt(var.pooled %*% t(rep(1, ncol(exprs))))

  # gamma.hat: batch mean of standardized data
  gamma.hat <- matrix(NA, n_batch, nrow(exprs))
  rownames(gamma.hat) <- levels(batch)
  colnames(gamma.hat) <- rownames(exprs)
  for (i in seq_len(n_batch)) {
    gamma.hat[i, ] <- rowMeans(s.data[, batches[[i]], drop = FALSE])
  }

  # delta.hat: batch variance of standardized data
  delta.hat <- matrix(NA, n_batch, nrow(exprs))
  rownames(delta.hat) <- levels(batch)
  colnames(delta.hat) <- rownames(exprs)
  for (i in seq_len(n_batch)) {
    delta.hat[i, ] <- apply(s.data[, batches[[i]], drop = FALSE], 1, var)
  }

  list(gamma = gamma.hat, delta = delta.hat,
       var_pooled = as.numeric(var.pooled),
       genes = rownames(exprs))
}

# -------------------------------------------------------------------
# Load matrices: pre-ComBat (imputed) and post-ComBat (corrected)
# -------------------------------------------------------------------
cat("Loading expression matrices...\n")

# Pre-ComBat: for parameter extraction
exprs_inter_pre <- as.matrix(read.delim(
  file.path(base_dir, "exprs_imputed_none.tsv"),
  row.names = 1, check.names = FALSE))
exprs_soft_pre <- as.matrix(read.delim(
  file.path(base_dir, "exprs_imputed_softimpute.tsv"),
  row.names = 1, check.names = FALSE))

# Post-ComBat: to see if correction differs on shared genes
exprs_inter_post <- as.matrix(read.delim(
  file.path(base_dir, "exprs_none_combat_ref.tsv"),
  row.names = 1, check.names = FALSE))
exprs_soft_post <- as.matrix(read.delim(
  file.path(base_dir, "exprs_softimpute_combat_ref.tsv"),
  row.names = 1, check.names = FALSE))

# Alias for backward compat with rest of script
exprs_inter <- exprs_inter_pre
exprs_soft  <- exprs_soft_pre

# Align to phenodata
common_samples <- intersect(colnames(exprs_inter), pheno$sample_id)
exprs_inter <- exprs_inter[, common_samples]
exprs_soft  <- exprs_soft[, common_samples]
pd <- pheno[match(common_samples, pheno$sample_id), ]

batch <- as.factor(pd$dataset_id)
mod <- model.matrix(~ pd$trimester)

cat("Intersection matrix:", nrow(exprs_inter), "genes x", ncol(exprs_inter), "samples\n")
cat("softImpute matrix:", nrow(exprs_soft), "genes x", ncol(exprs_soft), "samples\n")

# Remove zero-variance genes
keep_inter <- apply(exprs_inter, 1, var) > 0
keep_soft  <- apply(exprs_soft, 1, var) > 0
exprs_inter <- exprs_inter[keep_inter, ]
exprs_soft  <- exprs_soft[keep_soft, ]

cat("After variance filter: intersection=", nrow(exprs_inter),
    ", softImpute=", nrow(exprs_soft), "\n")

# -------------------------------------------------------------------
# Extract parameters
# -------------------------------------------------------------------
cat("\nExtracting ComBat parameters (intersection-only)...\n")
params_inter <- extract_combat_params(exprs_inter, batch, mod, ref_batch = "GSE100051")

cat("Extracting ComBat parameters (softImpute)...\n")
params_soft <- extract_combat_params(exprs_soft, batch, mod, ref_batch = "GSE100051")

# -------------------------------------------------------------------
# Compare on shared genes
# -------------------------------------------------------------------
shared_genes <- intersect(params_inter$genes, params_soft$genes)
cat("\nShared genes:", length(shared_genes), "\n")

batches <- rownames(params_inter$gamma)
cat("Batches:", paste(batches, collapse = ", "), "\n\n")

# Per-batch comparison
results <- list()
for (b in batches) {
  g_inter <- params_inter$gamma[b, shared_genes]
  g_soft  <- params_soft$gamma[b, shared_genes]
  d_inter <- params_inter$delta[b, shared_genes]
  d_soft  <- params_soft$delta[b, shared_genes]

  r_gamma <- cor(g_inter, g_soft, use = "complete.obs")
  r_delta <- cor(d_inter, d_soft, use = "complete.obs")
  max_dev_gamma <- max(abs(g_inter - g_soft), na.rm = TRUE)
  max_dev_delta <- max(abs(d_inter - d_soft), na.rm = TRUE)
  mae_gamma <- mean(abs(g_inter - g_soft), na.rm = TRUE)
  mae_delta <- mean(abs(d_inter - d_soft), na.rm = TRUE)

  cat(sprintf("Batch %s:\n", b))
  cat(sprintf("  gamma: r=%.6f  MAE=%.6f  max_dev=%.4f\n", r_gamma, mae_gamma, max_dev_gamma))
  cat(sprintf("  delta: r=%.6f  MAE=%.6f  max_dev=%.4f\n", r_delta, mae_delta, max_dev_delta))

  results[[b]] <- data.frame(
    batch = b,
    gamma_r = r_gamma, gamma_mae = mae_gamma, gamma_max_dev = max_dev_gamma,
    delta_r = r_delta, delta_mae = mae_delta, delta_max_dev = max_dev_delta,
    stringsAsFactors = FALSE
  )
}

results_df <- do.call(rbind, results)
write.csv(results_df, file.path(output_dir, "combat_param_comparison.csv"), row.names = FALSE)
cat("\nSaved:", file.path(output_dir, "combat_param_comparison.csv"), "\n")

# -------------------------------------------------------------------
# Scatter plots: gamma and delta for each batch
# -------------------------------------------------------------------
cat("\nGenerating scatter plots...\n")

png(file.path(output_dir, "combat_gamma_scatter.png"),
    width = 3000, height = 2000, res = 200)
par(mfrow = c(2, 3), mar = c(5, 5, 3, 1))
for (b in batches) {
  g_inter <- params_inter$gamma[b, shared_genes]
  g_soft  <- params_soft$gamma[b, shared_genes]
  r_val <- cor(g_inter, g_soft, use = "complete.obs")
  plot(g_inter, g_soft, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.2),
       xlab = "gamma (intersection)", ylab = "gamma (softImpute)",
       main = sprintf("%s (r=%.4f)", b, r_val),
       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1)
  abline(0, 1, col = "red", lwd = 1.5)
}
dev.off()
cat("Saved: combat_gamma_scatter.png\n")

png(file.path(output_dir, "combat_delta_scatter.png"),
    width = 3000, height = 2000, res = 200)
par(mfrow = c(2, 3), mar = c(5, 5, 3, 1))
for (b in batches) {
  d_inter <- params_inter$delta[b, shared_genes]
  d_soft  <- params_soft$delta[b, shared_genes]
  r_val <- cor(d_inter, d_soft, use = "complete.obs")
  plot(d_inter, d_soft, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.2),
       xlab = "delta (intersection)", ylab = "delta (softImpute)",
       main = sprintf("%s (r=%.4f)", b, r_val),
       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1)
  abline(0, 1, col = "red", lwd = 1.5)
}
dev.off()
cat("Saved: combat_delta_scatter.png\n")

# Summary of pre-ComBat parameter comparison
cat("\n=== Pre-ComBat Parameter Summary ===\n")
cat(sprintf("Overall gamma correlation range: %.6f - %.6f\n",
            min(results_df$gamma_r), max(results_df$gamma_r)))
cat(sprintf("Overall delta correlation range: %.6f - %.6f\n",
            min(results_df$delta_r), max(results_df$delta_r)))
cat(sprintf("Max gamma deviation across all batches: %.4f\n",
            max(results_df$gamma_max_dev)))
cat(sprintf("Max delta deviation across all batches: %.4f\n",
            max(results_df$delta_max_dev)))

# -------------------------------------------------------------------
# Post-ComBat comparison: do the corrected values for shared genes
# differ between intersection-only and softImpute pipelines?
# This is the real question: does adding 9,271 imputed genes to the
# ComBat input change the correction applied to the 8,260 shared genes?
# -------------------------------------------------------------------
cat("\n=== Post-ComBat Corrected Values (shared genes) ===\n")

shared_post <- intersect(rownames(exprs_inter_post), rownames(exprs_soft_post))
common_samp <- intersect(colnames(exprs_inter_post), colnames(exprs_soft_post))
cat("Shared genes post-ComBat:", length(shared_post), "\n")
cat("Shared samples:", length(common_samp), "\n")

inter_post <- exprs_inter_post[shared_post, common_samp]
soft_post  <- exprs_soft_post[shared_post, common_samp]

# Per-gene correlation of corrected values
diff_mat <- inter_post - soft_post
gene_mae <- rowMeans(abs(diff_mat))
gene_max_dev <- apply(abs(diff_mat), 1, max)

overall_r <- cor(as.vector(inter_post), as.vector(soft_post))
overall_mae <- mean(abs(diff_mat))
overall_max <- max(abs(diff_mat))
overall_rmse <- sqrt(mean(diff_mat^2))

cat(sprintf("Overall Pearson r: %.6f\n", overall_r))
cat(sprintf("Overall MAE: %.6f\n", overall_mae))
cat(sprintf("Overall RMSE: %.6f\n", overall_rmse))
cat(sprintf("Max absolute deviation: %.4f\n", overall_max))
cat(sprintf("Mean per-gene MAE: %.6f\n", mean(gene_mae)))
cat(sprintf("Median per-gene MAE: %.6f\n", median(gene_mae)))
cat(sprintf("95th percentile per-gene MAE: %.6f\n", quantile(gene_mae, 0.95)))

post_summary <- data.frame(
  metric = c("pearson_r", "mae", "rmse", "max_dev",
             "gene_mae_mean", "gene_mae_median", "gene_mae_p95"),
  value = c(overall_r, overall_mae, overall_rmse, overall_max,
            mean(gene_mae), median(gene_mae), quantile(gene_mae, 0.95)),
  stringsAsFactors = FALSE
)
write.csv(post_summary, file.path(output_dir, "post_combat_comparison.csv"), row.names = FALSE)
cat("\nSaved:", file.path(output_dir, "post_combat_comparison.csv"), "\n")

# Scatter plot of post-ComBat values (flattened)
set.seed(42)
idx <- sample.int(length(inter_post), min(50000, length(inter_post)))
png(file.path(output_dir, "post_combat_scatter.png"),
    width = 2000, height = 2000, res = 200)
par(mar = c(5, 5, 3, 1))
plot(as.vector(inter_post)[idx], as.vector(soft_post)[idx],
     pch = 16, cex = 0.2, col = rgb(0, 0, 0, 0.1),
     xlab = "Corrected expression (intersection-only ComBat-ref)",
     ylab = "Corrected expression (softImpute ComBat-ref)",
     main = sprintf("Post-ComBat values, %d shared genes (r=%.3f)", length(shared_post), overall_r),
     cex.main = 1.3, cex.lab = 1.2)
abline(0, 1, col = "red", lwd = 1.5)
dev.off()
cat("Saved: post_combat_scatter.png\n")

# Histogram of per-gene MAE
png(file.path(output_dir, "per_gene_mae_hist.png"),
    width = 2000, height = 1200, res = 200)
par(mar = c(5, 5, 3, 1))
hist(gene_mae, breaks = 50, col = "steelblue", border = "white",
     xlab = "Per-gene MAE (intersection vs softImpute post-ComBat)",
     ylab = "Number of genes",
     main = sprintf("Distribution of per-gene correction difference (n=%d)", length(shared_post)),
     cex.main = 1.3, cex.lab = 1.2)
abline(v = median(gene_mae), col = "red", lwd = 2, lty = 2)
legend("topright", sprintf("median = %.4f", median(gene_mae)),
       col = "red", lty = 2, lwd = 2, cex = 1.1, bty = "n")
dev.off()
cat("Saved: per_gene_mae_hist.png\n")

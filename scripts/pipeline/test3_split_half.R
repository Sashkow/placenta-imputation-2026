#!/usr/bin/env Rscript
#
# Test 3: Stratified split-half validation
# Split 148 samples into 2 halves maintaining trimester ratio, run DE on each.
#
# Usage: Rscript test3_split_half.R [--config=config_validation.yaml]

library(parallel)

script_dir <- if (length(grep("--file=", commandArgs(FALSE), value = TRUE)) > 0) {
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
} else {
  "scripts/pipeline"
}
source(file.path(script_dir, "subsampling_helpers.R"))

default_config <- file.path(script_dir, "config_validation.yaml")
config <- parse_config_arg(default_config)
n_iter <- parse_int_arg("n_iter", config$test3$n_iter %||% 50)
n_cores <- parse_int_arg("n_cores", config$parallel$n_cores %||% 10)
if (n_cores == 0) n_cores <- max(1, detectCores() - 2)

output_dir <- config$paths$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Test 3: Stratified split-half ===\n")
cat("Iterations:", n_iter, " Cores:", n_cores, "\n\n")

data <- load_pipeline_data(config)
first_pool <- data$first_trim_samples
second_pool <- data$second_trim_samples
n1 <- length(first_pool)
n2 <- length(second_pool)
half1 <- floor(n1 / 2)
half2 <- floor(n2 / 2)

cat("\n1st-trim pool:", n1, " 2nd-trim pool:", n2, "\n")
cat("Split: half A =", half1, "+", half2, " half B =",
    n1 - half1, "+", n2 - half2, "\n\n")

t0 <- Sys.time()

iter_results <- mclapply(seq_len(n_iter), function(iter) {
  set.seed(iter * 3000)

  perm_1st <- sample(first_pool)
  perm_2nd <- sample(second_pool)

  half_a_ids <- c(perm_1st[seq_len(half1)], perm_2nd[seq_len(half2)])
  half_b_ids <- c(perm_1st[(half1 + 1):n1], perm_2nd[(half2 + 1):n2])

  sub_a <- subsample_exprs_list(data$exprs_list, data$phenodata, half_a_ids)
  sub_b <- subsample_exprs_list(data$exprs_list, data$phenodata, half_b_ids)

  sink(tempfile())
  on.exit(sink(), add = TRUE)
  res_a <- run_lean_pipeline(sub_a$exprs_list, sub_a$phenodata,
                              config, data$ref_imputed)
  res_b <- run_lean_pipeline(sub_b$exprs_list, sub_b$phenodata,
                              config, data$ref_imputed)
  sink()

  metrics_a <- compute_subsample_metrics(res_a, data$ref_de, data$ref_sig_genes, data$ref_sig_genes_fdr_only)
  metrics_b <- compute_subsample_metrics(res_b, data$ref_de, data$ref_sig_genes, data$ref_sig_genes_fdr_only)
  pairwise <- compute_pairwise_metrics(res_a, res_b)

  data.frame(
    iter = iter,
    n_deg_a = res_a$n_sig, n_deg_b = res_b$n_sig,
    jaccard_a_vs_full = metrics_a$jaccard_vs_full,
    jaccard_b_vs_full = metrics_b$jaccard_vs_full,
    logfc_r_a_vs_full = metrics_a$logfc_pearson_vs_full,
    logfc_ccc_a_vs_full = metrics_a$logfc_ccc_vs_full,
    logfc_r_b_vs_full = metrics_b$logfc_pearson_vs_full,
    logfc_ccc_b_vs_full = metrics_b$logfc_ccc_vs_full,
    overlap_a_vs_full = metrics_a$overlap_vs_full,
    overlap_b_vs_full = metrics_b$overlap_vs_full,
    overlap_fdr_a_vs_full = metrics_a$overlap_fdr_only,
    overlap_fdr_b_vs_full = metrics_b$overlap_fdr_only,
    jaccard_a_vs_b = pairwise$jaccard_between,
    logfc_pearson_a_vs_b = pairwise$logfc_pearson_between,
    logfc_ccc_a_vs_b = pairwise$logfc_ccc_between,
    logfc_spearman_a_vs_b = pairwise$logfc_spearman_between,
    same_direction_a_vs_b = pairwise$same_direction_pct,
    n_deg_fdr_a = res_a$n_sig_fdr_only,
    n_deg_fdr_b = res_b$n_sig_fdr_only,
    jaccard_fdr_a_vs_b = pairwise$jaccard_fdr_only_between,
    jaccard_fdr_a_vs_full = metrics_a$jaccard_fdr_only,
    jaccard_fdr_b_vs_full = metrics_b$jaccard_fdr_only,
    n_deg_lfc15_a = res_a$n_sig_lfc15,
    n_deg_lfc15_b = res_b$n_sig_lfc15,
    jaccard_lfc15_a_vs_b = pairwise$jaccard_lfc15_between,
    n_deg_lfc20_a = res_a$n_sig_lfc20,
    n_deg_lfc20_b = res_b$n_sig_lfc20,
    jaccard_lfc20_a_vs_b = pairwise$jaccard_lfc20_between,
    status_a = res_a$status, status_b = res_b$status,
    stringsAsFactors = FALSE
  )
}, mc.cores = n_cores)

results_df <- do.call(rbind, iter_results)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
ok <- sum(results_df$status_a == "ok" & results_df$status_b == "ok")
cat(sprintf("ok: %d/%d  median Jaccard(A,B): %.3f  median logFC_r(A,B): %.3f  (%.1fs)\n",
            ok, n_iter,
            median(results_df$jaccard_a_vs_b, na.rm = TRUE),
            median(results_df$logfc_pearson_a_vs_b, na.rm = TRUE),
            elapsed))

out_file <- file.path(output_dir, "test3_split_half.tsv")
write.table(results_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved:", out_file, "\n")
cat("Test 3 complete.\n")

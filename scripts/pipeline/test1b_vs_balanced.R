#!/usr/bin/env Rscript
#
# Test 1b: Compare subsamples against balanced 2-dataset reference
# Same subsampling as test1 but metrics computed vs GSE100051+GSE9984
# instead of vs the full integration.
#
# Usage: Rscript test1b_vs_balanced.R [--config=config_validation.yaml]

library(parallel)

script_dir <- if (length(grep("--file=", commandArgs(FALSE), value = TRUE)) > 0) {
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
} else {
  "scripts/pipeline"
}
source(file.path(script_dir, "subsampling_helpers.R"))

default_config <- file.path(script_dir, "config_validation.yaml")
config <- parse_config_arg(default_config)
n_iter <- parse_int_arg("n_iter", config$test1$n_iter %||% 30)
n_cores <- parse_int_arg("n_cores", config$parallel$n_cores %||% 10)
if (n_cores == 0) n_cores <- max(1, detectCores() - 2)

output_dir <- config$paths$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bal_de_path <- config$paths$balanced_reference_de
if (is.null(bal_de_path) || !file.exists(bal_de_path))
  stop("balanced_reference_de not set or file not found: ", bal_de_path)

cat("=== Test 1b: vs balanced reference ===\n")
cat("Iterations per size:", n_iter, " Cores:", n_cores, "\n")

data <- load_pipeline_data(config)

bal_de <- read.delim(bal_de_path, stringsAsFactors = FALSE)
fdr_thresh <- config$thresholds$fdr %||% 0.05
logfc_thresh <- config$thresholds$logfc %||% 1.0
bal_sig <- bal_de$gene[
  !is.na(bal_de$adj.P.Val) & !is.na(bal_de$logFC) &
  bal_de$adj.P.Val < fdr_thresh & abs(bal_de$logFC) > logfc_thresh
]
bal_sig_fdr <- bal_de$gene[
  !is.na(bal_de$adj.P.Val) & bal_de$adj.P.Val < fdr_thresh
]

cat("Balanced ref DEGs:", length(bal_sig),
    " FDR-only:", length(bal_sig_fdr), "\n")
cat("Full ref DEGs:", length(data$ref_sig_genes), "\n\n")

sizes <- config$test1$sizes_1st_trim
first_pool <- data$first_trim_samples
second_pool <- data$second_trim_samples

cat("1st-trim pool:", length(first_pool),
    " 2nd-trim pool:", length(second_pool), "\n")
cat("Sizes:", paste(sizes, collapse = ", "), "\n\n")

all_rows <- list()

for (N in sizes) {
  if (N > length(first_pool)) {
    cat("Skipping N=", N, "\n")
    next
  }
  cat("--- N_1st =", N, "---\n")
  t0 <- Sys.time()

  iter_results <- mclapply(seq_len(n_iter), function(iter) {
    set.seed(iter * 1000 + N)
    sampled_1st <- sample(first_pool, N, replace = FALSE)
    keep <- c(sampled_1st, second_pool)
    sub <- subsample_exprs_list(data$exprs_list,
                                data$phenodata, keep)

    sink(tempfile())
    on.exit(sink(), add = TRUE)
    res <- run_lean_pipeline(sub$exprs_list, sub$phenodata,
                              config, data$ref_imputed)
    sink()

    vs_full <- compute_subsample_metrics(
      res, data$ref_de, data$ref_sig_genes,
      data$ref_sig_genes_fdr_only
    )
    vs_bal <- compute_subsample_metrics(
      res, bal_de, bal_sig, bal_sig_fdr
    )

    data.frame(
      N_1st = N, N_total = N + length(second_pool),
      iter = iter,
      n_deg = res$n_sig,
      jaccard_vs_full = vs_full$jaccard_vs_full,
      logfc_r_vs_full = vs_full$logfc_pearson_vs_full,
      logfc_ccc_vs_full = vs_full$logfc_ccc_vs_full,
      jaccard_vs_balanced = vs_bal$jaccard_vs_full,
      logfc_r_vs_balanced = vs_bal$logfc_pearson_vs_full,
      logfc_ccc_vs_balanced = vs_bal$logfc_ccc_vs_full,
      overlap_vs_balanced = vs_bal$overlap_vs_full,
      same_dir_vs_balanced = vs_bal$same_direction_pct,
      jaccard_fdr_vs_full = vs_full$jaccard_fdr_only,
      jaccard_fdr_vs_balanced = vs_bal$jaccard_fdr_only,
      overlap_fdr_vs_balanced = vs_bal$overlap_fdr_only,
      status = res$status,
      stringsAsFactors = FALSE
    )
  }, mc.cores = n_cores)

  rows <- do.call(rbind, iter_results)
  all_rows[[length(all_rows) + 1]] <- rows

  ok <- sum(rows$status == "ok", na.rm = TRUE)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0,
                                        units = "secs")), 1)
  cat(sprintf(
    "  ok: %d/%d  J_full: %.3f  J_bal: %.3f  (%.1fs)\n",
    ok, n_iter,
    median(rows$jaccard_vs_full, na.rm = TRUE),
    median(rows$jaccard_vs_balanced, na.rm = TRUE),
    elapsed
  ))
}

results_df <- do.call(rbind, all_rows)
out_file <- file.path(output_dir, "test1b_vs_balanced.tsv")
write.table(results_df, out_file, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nSaved:", out_file, "\n")
cat("Test 1b complete.\n")

#!/usr/bin/env Rscript
#
# Test 2: Balanced subsampling
# Equal N from each trimester.
#
# Usage: Rscript test2_balanced_subsample.R [--config=config_validation.yaml]

library(parallel)

script_dir <- if (length(grep("--file=", commandArgs(FALSE), value = TRUE)) > 0) {
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
} else {
  "scripts/pipeline"
}
source(file.path(script_dir, "subsampling_helpers.R"))

default_config <- file.path(script_dir, "config_validation.yaml")
config <- parse_config_arg(default_config)
n_iter <- parse_int_arg("n_iter", config$test2$n_iter %||% 30)
n_cores <- parse_int_arg("n_cores", config$parallel$n_cores %||% 10)
if (n_cores == 0) n_cores <- max(1, detectCores() - 2)

output_dir <- config$paths$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Test 2: Balanced subsampling ===\n")
cat("Iterations per size:", n_iter, " Cores:", n_cores, "\n\n")

data <- load_pipeline_data(config)
sizes <- config$test2$sizes_per_trim
first_pool <- data$first_trim_samples
second_pool <- data$second_trim_samples

cat("\n1st-trim pool:", length(first_pool), " 2nd-trim pool:", length(second_pool), "\n")
cat("Sizes per trimester:", paste(sizes, collapse = ", "), "\n\n")

all_rows <- list()

for (n_per in sizes) {
  if (n_per > length(first_pool) || n_per > length(second_pool)) {
    cat("Skipping n_per=", n_per, " (insufficient samples)\n")
    next
  }
  cat("--- n_per_trim =", n_per, "---\n")
  t0 <- Sys.time()

  iter_results <- mclapply(seq_len(n_iter), function(iter) {
    set.seed(iter * 2000 + n_per)
    sampled_1st <- sample(first_pool, n_per, replace = FALSE)
    sampled_2nd <- sample(second_pool, n_per, replace = FALSE)
    keep <- c(sampled_1st, sampled_2nd)
    sub <- subsample_exprs_list(data$exprs_list, data$phenodata, keep)

    sink(tempfile())
    on.exit(sink(), add = TRUE)
    res <- run_lean_pipeline(sub$exprs_list, sub$phenodata, config,
                              data$ref_imputed)
    sink()

    metrics <- compute_subsample_metrics(res, data$ref_de, data$ref_sig_genes, data$ref_sig_genes_fdr_only)
    cbind(
      data.frame(n_per_trim = n_per, N_total = 2 * n_per, iter = iter,
                 stringsAsFactors = FALSE),
      metrics,
      data.frame(status = res$status, stringsAsFactors = FALSE)
    )
  }, mc.cores = n_cores)

  rows <- do.call(rbind, iter_results)
  all_rows[[length(all_rows) + 1]] <- rows

  ok <- sum(rows$status == "ok", na.rm = TRUE)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("  ok: %d/%d  median DEGs: %.0f  median Jaccard: %.3f  (%.1fs)\n",
              ok, n_iter,
              median(rows$n_deg, na.rm = TRUE),
              median(rows$jaccard_vs_full, na.rm = TRUE),
              elapsed))
}

results_df <- do.call(rbind, all_rows)
out_file <- file.path(output_dir, "test2_balanced_subsample.tsv")
write.table(results_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved:", out_file, "\n")
cat("\nTest 2 complete.\n")

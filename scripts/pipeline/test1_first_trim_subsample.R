#!/usr/bin/env Rscript
#
# Test 1: First-trimester subsampling
# Fix all 16 2nd-trimester samples, vary 1st-trimester N.
#
# Usage: Rscript test1_first_trim_subsample.R [--config=config_validation.yaml]

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

cat("=== Test 1: First-trimester subsampling ===\n")
cat("Iterations per size:", n_iter, " Cores:", n_cores, "\n\n")

data <- load_pipeline_data(config)
sizes <- config$test1$sizes_1st_trim
first_pool <- data$first_trim_samples
second_pool <- data$second_trim_samples

cat("\n1st-trim pool:", length(first_pool), " 2nd-trim pool:", length(second_pool), "\n")
cat("Sizes to test:", paste(sizes, collapse = ", "), "\n\n")

all_rows <- list()
within_jaccard_rows <- list()

for (N in sizes) {
  if (N > length(first_pool)) {
    cat("Skipping N=", N, " (only", length(first_pool), "available)\n")
    next
  }
  cat("--- N_1st =", N, "---\n")
  t0 <- Sys.time()

  iter_results <- mclapply(seq_len(n_iter), function(iter) {
    set.seed(iter * 1000 + N)
    sampled_1st <- sample(first_pool, N, replace = FALSE)
    keep <- c(sampled_1st, second_pool)
    sub <- subsample_exprs_list(data$exprs_list, data$phenodata, keep)

    sink(tempfile())
    on.exit(sink(), add = TRUE)
    res <- run_lean_pipeline(sub$exprs_list, sub$phenodata, config,
                              data$ref_imputed)
    sink()

    metrics <- compute_subsample_metrics(res, data$ref_de, data$ref_sig_genes, data$ref_sig_genes_fdr_only)
    row <- cbind(
      data.frame(N_1st = N, N_total = N + length(second_pool), iter = iter,
                 stringsAsFactors = FALSE),
      metrics,
      data.frame(status = res$status, stringsAsFactors = FALSE)
    )
    list(row = row, sig_genes = res$sig_genes)
  }, mc.cores = n_cores)

  rows <- do.call(rbind, lapply(iter_results, `[[`, "row"))
  all_rows[[length(all_rows) + 1]] <- rows

  ok <- sum(rows$status == "ok", na.rm = TRUE)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  cat(sprintf("  ok: %d/%d  median DEGs: %.0f  median Jaccard: %.3f  (%.1fs)\n",
              ok, n_iter,
              median(rows$n_deg, na.rm = TRUE),
              median(rows$jaccard_vs_full, na.rm = TRUE),
              elapsed))

  if (isTRUE(config$test1$compute_within_size_jaccard)) {
    sig_lists <- lapply(iter_results, `[[`, "sig_genes")
    valid <- which(sapply(sig_lists, function(x) length(x) > 0))
    if (length(valid) >= 2) {
      pairs <- combn(valid, 2)
      pw_j <- apply(pairs, 2, function(idx) {
        jaccard(sig_lists[[idx[1]]], sig_lists[[idx[2]]])
      })
      within_jaccard_rows[[length(within_jaccard_rows) + 1]] <- data.frame(
        N_1st = N,
        mean_pairwise_jaccard = mean(pw_j, na.rm = TRUE),
        sd_pairwise_jaccard = sd(pw_j, na.rm = TRUE),
        median_pairwise_jaccard = median(pw_j, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}

results_df <- do.call(rbind, all_rows)
out_file <- file.path(output_dir, "test1_first_trim_subsample.tsv")
write.table(results_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved:", out_file, "\n")

if (length(within_jaccard_rows) > 0) {
  wj_df <- do.call(rbind, within_jaccard_rows)
  wj_file <- file.path(output_dir, "test1_within_size_jaccard.tsv")
  write.table(wj_df, wj_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved:", wj_file, "\n")
}

cat("\nTest 1 complete.\n")

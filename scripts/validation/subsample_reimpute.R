#!/usr/bin/env Rscript
##
## Subsampling with imputation repeated in every iteration.
##
## The article's subsampling test reuses the imputed matrix of the full cohort.
## Here every subsample is run twice: with the reused matrix ("reuse", as in
## the article) and with softImpute fitted on the subsample alone ("reimpute").
## Settings: N first-trimester samples x all 15 or a subsample of second-
## trimester samples. Seeds match scripts/validation/subsample_first_trimester.R,
## so the "reuse" arm with all second-trimester samples repeats that test.
##
## Run from the repository root:
##   Rscript scripts/validation/subsample_reimpute.R [--sizes=10,30,60]
##           [--n_second=15,10] [--n_iter=30] [--n_cores=5]
## Output: data/pipeline/validation/subsample_reimpute.tsv and _summary.csv.

library(parallel)
source("scripts/lib/subsampling.R")

config  <- config_from_cli("config/validation.yaml")
out_dir <- get_arg("out", config$paths$output %||% "data/pipeline/validation")
sizes   <- as.integer(split_arg("sizes", c("10", "30", "60")))
n_seconds <- as.integer(split_arg("n_second", c("15", "10")))
n_iter  <- int_arg("n_iter", 30)
n_cores <- int_arg("n_cores", 5)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data <- load_pipeline_data(config)
first_pool <- data$first_trim_samples; second_pool <- data$second_trim_samples
cat("1st-trim pool:", length(first_pool), " 2nd-trim pool:", length(second_pool), "\n")

one <- function(N, n2, iter) {
  set.seed(iter * 1000 + N)
  s1 <- sample(first_pool, N, replace = FALSE)
  s2 <- if (n2 >= length(second_pool)) second_pool else sample(second_pool, n2, replace = FALSE)
  sub <- subsample_exprs_list(data$exprs_list, data$phenodata, c(s1, s2))
  arms <- list(reuse = data$ref_imputed, reimpute = NULL)
  sink(tempfile())
  res <- lapply(arms, function(ri) run_lean_pipeline(sub$exprs_list, sub$phenodata, config, ri))
  sink()
  rows <- lapply(names(res), function(a) cbind(
    data.frame(N_1st = N, N_2nd = length(s2), iter = iter, arm = a, stringsAsFactors = FALSE),
    compute_subsample_metrics(res[[a]], data$ref_de, data$ref_sig_genes, data$ref_sig_genes_fdr_only),
    data.frame(status = res[[a]]$status, stringsAsFactors = FALSE)))
  rows <- do.call(rbind, rows)
  rows$jaccard_reuse_vs_reimpute <- if (all(rows$status == "ok"))
    jaccard(res$reuse$sig_genes, res$reimpute$sig_genes) else NA
  rows
}

all_rows <- list()
for (n2 in n_seconds) for (N in sizes) {
  t0 <- Sys.time()
  r <- mclapply(seq_len(n_iter), function(i) one(N, n2, i), mc.cores = n_cores)
  bad <- sapply(r, inherits, "try-error")
  if (any(bad)) cat("  failed iterations:", sum(bad), "\n")
  r <- do.call(rbind, r[!bad])
  all_rows[[length(all_rows) + 1]] <- r
  for (a in c("reuse", "reimpute")) {
    x <- r[r$arm == a, ]
    cat(sprintf("N_1st=%3d N_2nd=%2d %-8s ok %d/%d  median DEGs %.0f  retention %.3f  Jaccard %.3f\n",
                N, n2, a, sum(x$status == "ok"), n_iter, median(x$n_deg, na.rm = TRUE),
                median(x$overlap_vs_full, na.rm = TRUE), median(x$jaccard_vs_full, na.rm = TRUE)))
  }
  cat(sprintf("  (%.0f s)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  write.table(do.call(rbind, all_rows), file.path(out_dir, "subsample_reimpute.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

d <- do.call(rbind, all_rows)
summ <- aggregate(cbind(n_deg, overlap_vs_full, jaccard_vs_full, logfc_ccc_vs_full,
                        overlap_fdr_only, jaccard_reuse_vs_reimpute) ~ N_1st + N_2nd + arm,
                  data = d, FUN = function(x) median(x, na.rm = TRUE), na.action = na.pass)
summ$n_ok <- aggregate(status ~ N_1st + N_2nd + arm, data = d, FUN = function(s) sum(s == "ok"))$status
write.csv(summ, file.path(out_dir, "subsample_reimpute_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 3)

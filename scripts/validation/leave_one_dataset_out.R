#!/usr/bin/env Rscript
# Leave-one-dataset-out (LODO): drop each of the six datasets in turn, re-run
# batch correction + DE on the rest, compare with the full-cohort DEG list.
#
# Run from repo root:
#   Rscript scripts/validation/leave_one_dataset_out.R [--config=config/validation.yaml]
#
# Notes
# - Dropping GSE100051 removes the ComBat reference batch; the lean pipeline then
#   falls back to plain ComBat (recorded in the `combat` column).
# - Reference DEGs measured only in the dropped dataset can no longer be tested;
#   `retention_testable` uses only reference DEGs still tested in that run.

source("scripts/validation/_lodo_permutation_common.R")

check_full_run()

ref_sig <- DATA$ref_sig_genes
ref_de  <- DATA$ref_de
rows <- list()
# Extra run: GSE100051 dropped but ComBat-ref kept, with GSE9984 (the other
# dataset holding both trimesters) as reference -- separates losing the data
# from losing the reference batch.
runs <- c(setNames(names(DATA$exprs_list), names(DATA$exprs_list)),
          "GSE100051 (ref GSE9984)" = "GSE100051")
cfg0 <- config
for (run_id in names(runs)) {
  ds <- runs[[run_id]]
  config <- cfg0
  if (grepl("ref GSE9984", run_id)) config$normalization$ref_batch <- "GSE9984"
  keep <- DATA$phenodata[[sample_col]][DATA$phenodata[[ds_col]] != ds]
  sub  <- subsample_exprs_list(DATA$exprs_list, DATA$phenodata, keep)
  t0 <- Sys.time()
  r <- run_on(sub$exprs_list, sub$phenodata)
  if (r$status != "ok") stop(ds, ": ", r$status)
  m <- compute_subsample_metrics(r, ref_de, ref_sig, DATA$ref_sig_genes_fdr_only)
  testable <- intersect(ref_sig, r$de_results$gene)
  sig_ref_dir <- sign(setNames(ref_de$logFC, ref_de$gene)[testable])
  sub_lfc <- setNames(r$de_results$logFC, r$de_results$gene)[testable]
  sub_fdr <- setNames(r$de_results$adj.P.Val, r$de_results$gene)[testable]
  rows[[run_id]] <- data.frame(
    run = run_id, dropped = ds,
    n_samples = nrow(sub$phenodata),
    n_first = sum(sub$phenodata[[grp_col]] == config$phenotype$baseline),
    n_second = sum(sub$phenodata[[grp_col]] == config$phenotype$contrast),
    combat = if (config$normalization$ref_batch %in% names(sub$exprs_list))
               paste("ref", config$normalization$ref_batch) else "plain (ref dropped)",
    n_genes_tested = m$n_genes_tested,
    n_deg = m$n_deg,
    jaccard = m$jaccard_vs_full,
    retention = m$overlap_vs_full,
    ref_degs_testable = length(testable),
    retention_testable = length(intersect(r$sig_genes, testable)) / length(testable),
    ref_degs_same_direction = mean(sign(sub_lfc) == sig_ref_dir, na.rm = TRUE),
    ref_degs_fdr05 = mean(sub_fdr < 0.05, na.rm = TRUE),
    ccc_logfc_all = m$logfc_ccc_vs_full,
    new_degs = length(setdiff(r$sig_genes, ref_sig)),
    stringsAsFactors = FALSE)
  cat(sprintf("drop %-24s %3d samples (%d/%d) | %4d DEGs | Jaccard %.3f | retention %.3f (testable %.3f) | CCC %.3f | %.0fs\n",
              run_id, rows[[run_id]]$n_samples, rows[[run_id]]$n_first, rows[[run_id]]$n_second, m$n_deg,
              m$jaccard_vs_full, m$overlap_vs_full, rows[[run_id]]$retention_testable,
              m$logfc_ccc_vs_full, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
config <- cfg0
tab <- do.call(rbind, rows)
write.table(tab, file.path(OUT_DIR, "lodo_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(tab, row.names = FALSE)
cat("Saved:", file.path(OUT_DIR, "lodo_summary.tsv"), "\n")

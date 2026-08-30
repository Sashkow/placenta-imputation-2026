#!/usr/bin/env Rscript
##
## Full-data reference run per imputer (arm F of the unified holdout).
##
## For imputer X: merged matrix (nothing hidden) -> X -> ComBat-ref with the
## production covariates -> limma with weight 0 on originally-missing cells.
## The resulting DE table is the ground truth against which the masked arm
## for imputer X is scored, so masking effects are not confounded with how X
## differs from softImpute. The softImpute reference is checked against the
## shipped production DE table. A side table compares the five reference DEG
## lists with softImpute's ("does imputer choice change the published result?").
##
## Usage:
##   Rscript scripts/holdout_reference_runs.R [--methods=a,b] [--force]

source(file.path({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "."
}, "holdout_common.R"))

cat("=== Unified holdout: reference runs ===\n")
cat("Companion:", COMPANION, "\nOutput:", OUT_DIR, "\nMethods:",
    paste(METHODS, collapse = ", "), "\n\n")

D <- prepare_data(config)
tools <- make_de_tools(D, config)

for (m in METHODS) {
  p <- reference_paths(m)
  if (file.exists(p$de) && file.exists(p$corrected) && !FORCE) {
    log_msg("reference ", m, ": exists, skipping (use --force to redo)")
    next
  }
  log_msg("reference ", m, ": imputing")
  t0 <- Sys.time()
  imp <- run_imputer(m, D$incomplete, method_cfg(m))
  t_imp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  log_msg("reference ", m, ": ComBat-ref")
  corrected <- tryCatch(tools$combat_correct(imp$matrix), error = function(e) e)
  if (inherits(corrected, "error")) {
    log_msg("reference ", m, ": ComBat FAILED: ", conditionMessage(corrected))
    write.csv(data.frame(method = m, status = "combat_error",
                         message = conditionMessage(corrected), stringsAsFactors = FALSE),
              p$meta, row.names = FALSE)
    next
  }
  de <- tools$de_fit(corrected, D$orig_mask)
  de$sig <- sig_call(de$logFC, de$adj.P.Val)
  saveRDS(corrected, p$corrected)
  write.csv(de, p$de, row.names = FALSE)
  write.csv(data.frame(method = m, status = "ok", message = "",
                       n_genes = nrow(de), n_deg = sum(de$sig),
                       n_nonfinite = sum(!is.finite(de$logFC)),
                       imputation_seconds = round(t_imp, 1), stringsAsFactors = FALSE),
            p$meta, row.names = FALSE)
  log_msg("reference ", m, ": ", nrow(de), " genes, ", sum(de$sig), " DEGs, ",
          sum(!is.finite(de$logFC)), " non-finite logFC")
}

## ---------------------------------------------- softImpute vs shipped table
if ("softimpute" %in% METHODS && file.exists(reference_paths("softimpute")$de)) {
  ref <- read.csv(reference_paths("softimpute")$de, stringsAsFactors = FALSE)
  shipped <- read.delim(config$paths$reference_de, stringsAsFactors = FALSE)
  s_sig <- shipped$gene[sig_call(shipped$logFC, shipped$adj.P.Val)]
  r_sig <- ref$gene[ref$sig]
  common <- intersect(ref$gene, shipped$gene)
  lfc_r <- cor(ref$logFC[match(common, ref$gene)], shipped$logFC[match(common, shipped$gene)])
  only_ref <- setdiff(r_sig, s_sig); only_ship <- setdiff(s_sig, r_sig)
  chk <- data.frame(
    shipped_degs = length(s_sig), reference_degs = length(r_sig),
    shared_degs = length(intersect(s_sig, r_sig)),
    jaccard = length(intersect(s_sig, r_sig)) / length(union(s_sig, r_sig)),
    logfc_r_shared_genes = lfc_r, n_shared_genes = length(common),
    only_reference = paste(only_ref, collapse = ";"),
    only_shipped = paste(only_ship, collapse = ";"),
    stringsAsFactors = FALSE)
  write.csv(chk, file.path(OUT_DIR, "reference", "softimpute_vs_shipped.csv"), row.names = FALSE)
  cat(sprintf("\nsoftImpute reference vs shipped: %d vs %d DEGs, %d shared (Jaccard %.3f), logFC r = %.4f\n",
              chk$reference_degs, chk$shipped_degs, chk$shared_degs, chk$jaccard, lfc_r))
  if (length(only_ref) || length(only_ship)) {
    cat("  only in reference:", if (length(only_ref)) paste(only_ref, collapse = " ") else "-", "\n")
    cat("  only in shipped:  ", if (length(only_ship)) paste(only_ship, collapse = " ") else "-", "\n")
    d <- shipped[match(c(only_ref, only_ship), shipped$gene), c("gene", "logFC", "adj.P.Val")]
    cat("  shipped |logFC| of the discordant genes: ",
        paste(sprintf("%.3f", abs(d$logFC)), collapse = " "), "\n")
  }
}

## ---------------------------------------------- imputer comparison table
avail <- METHODS[sapply(METHODS, function(m) file.exists(reference_paths(m)$de))]
if (length(avail) && "softimpute" %in% avail) {
  base <- read.csv(reference_paths("softimpute")$de, stringsAsFactors = FALSE)
  b_sig <- base$gene[base$sig]
  rows <- lapply(avail, function(m) {
    de <- read.csv(reference_paths(m)$de, stringsAsFactors = FALSE)
    meta <- read.csv(reference_paths(m)$meta, stringsAsFactors = FALSE)
    m_sig <- de$gene[de$sig]
    common <- intersect(de$gene[is.finite(de$logFC)], base$gene[is.finite(base$logFC)])
    x <- base$logFC[match(common, base$gene)]; y <- de$logFC[match(common, de$gene)]
    data.frame(method = m, label = METHOD_LABELS[m],
               n_genes_tested = sum(is.finite(de$logFC)),
               n_nonfinite = sum(!is.finite(de$logFC)),
               n_deg = length(m_sig),
               shared_with_softimpute = length(intersect(m_sig, b_sig)),
               jaccard_vs_softimpute = length(intersect(m_sig, b_sig)) / length(union(m_sig, b_sig)),
               logfc_r_vs_softimpute = cor(x, y),
               logfc_ccc_vs_softimpute = lin_ccc(x, y),
               imputation_seconds = meta$imputation_seconds %||% NA,
               stringsAsFactors = FALSE, row.names = NULL)
  })
  summ <- do.call(rbind, rows)
  write.csv(summ, file.path(OUT_DIR, "reference_runs_summary.csv"), row.names = FALSE)
  cat("\n=== Reference runs: DEGs per imputer ===\n")
  print(summ[, c("label", "n_genes_tested", "n_nonfinite", "n_deg", "shared_with_softimpute",
                 "jaccard_vs_softimpute", "logfc_ccc_vs_softimpute")], row.names = FALSE, digits = 3)
}
cat("\nWritten to:", file.path(OUT_DIR, "reference"), "\n")

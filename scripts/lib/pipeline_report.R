# Diagnostic outputs of a pipeline run: PCA and missingness plots, the
# comparison of method combinations, and the text summary.

# PCA of the samples: colour = comparison group, symbol = dataset.
plot_pca <- function(exprs, pd, config, file, title) {
  pc <- prcomp(t(exprs), center = TRUE, scale. = FALSE)
  var_pct <- 100 * summary(pc)$importance[2, ]
  group <- pd[[config$phenotype$group_column]]; dataset <- as.character(pd[[DATASET_ID]])
  groups <- unique(group); datasets <- unique(dataset)
  palette <- c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
  col <- setNames(palette[seq_along(groups)], groups)
  pch <- setNames(c(0:4, 6:8, 15:18)[seq_along(datasets)], datasets)
  png(file, width = 1400, height = 900, res = 120)
  par(mar = c(5, 5, 3, 2))
  plot(pc$x[, 1], pc$x[, 2], col = col[group], pch = pch[dataset], cex = 1.3,
       xlab = sprintf("PC1 (%.1f%%)", var_pct[1]), ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
       main = sprintf("PCA — %s (%d genes, %d samples)", title, nrow(exprs), ncol(exprs)))
  legend("bottomright",
         legend = c(sprintf("%s (n=%d)", groups, table(factor(group, levels = groups))), "",
                    sprintf("%s (n=%d)", datasets, table(factor(dataset, levels = datasets)))),
         col = c(col[groups], NA, rep("grey30", length(datasets))),
         pch = c(rep(15, length(groups)), NA, pch[datasets]),
         pt.cex = c(rep(1.5, length(groups)), NA, rep(1.3, length(datasets))),
         cex = 0.6, bty = "n")
  dev.off()
}

# One row per imputation x batch-correction combination, and the number of
# DEGs each pair of combinations shares.
write_method_comparison <- function(results, output_dir) {
  comparison <- do.call(rbind, lapply(names(results), function(key) with(results[[key]],
    data.frame(method = key, imputation = imputation, normalization = normalization,
               n_genes = nrow(de), n_significant = nrow(significant),
               n_up = sum(significant$logFC > 0), n_down = sum(significant$logFC < 0),
               stringsAsFactors = FALSE))))
  write.csv(comparison, file.path(output_dir, "method_comparison.csv"), row.names = FALSE)
  if (length(results) > 1) {
    degs <- lapply(results, function(r) r$significant$gene)
    overlap <- sapply(degs, function(a) sapply(degs, function(b) length(intersect(a, b))))
    write.csv(overlap, file.path(output_dir, "deg_overlap_matrix.csv"))
  }
  print(comparison)
  comparison
}

# Missingness "staircase" (genes x samples, sorted by coverage) with each
# combination's FDR as the row colour. `incomplete_matrix` should still hold
# the genes removed by the safety filter, so they show up in the plot.
plot_staircases <- function(incomplete_matrix, exprs_list, phenodata, results, config,
                            coverage_threshold, output_dir) {
  sample_dataset <- setNames(rep(names(exprs_list), sapply(exprs_list, ncol)),
                             unlist(lapply(exprs_list, colnames)))
  group <- phenodata[[config$phenotype$group_column]][match(names(sample_dataset), phenodata[[SAMPLE_ID]])]
  title <- paste0(config$phenotype$baseline, " vs ", config$phenotype$contrast)
  for (key in names(results)) {
    fdr <- setNames(results[[key]]$de$adj.P.Val, as.character(results[[key]]$de$gene))
    staircase <- prepare_staircase(incomplete_matrix, group, fdr, config$phenotype$baseline,
                                   config$phenotype$contrast, coverage_threshold = coverage_threshold)
    png(file.path(output_dir, paste0("na_staircase_", key, ".png")), width = 2400, height = 1800, res = 200)
    par(mar = c(4, 5, 2, 1))
    render_staircase(staircase, paste0(title, " [FDR: ", key, "]"), sample_dataset)
    dev.off()
  }
}

# summary.txt: what was run and the headline counts.
write_summary <- function(output_dir, config_path, config, coverage, dropped_genes,
                          comparison, validation) {
  n_imputed <- nrow(coverage)
  lines <- c(
    "Pipeline run summary", "====================", "",
    paste("Date:", Sys.time()), paste("Config:", config_path), "",
    paste("Comparison:", config$phenotype$contrast, "vs", config$phenotype$baseline),
    paste("Datasets:", paste(config$files$datasets, collapse = ", ")),
    paste("Minimum datasets per gene:", config$coverage$min_datasets), "",
    "Method comparison:", capture.output(print(comparison)), "",
    "Imputed-gene coverage per comparison group:",
    sprintf("Genes needing imputation somewhere:    %d", n_imputed),
    sprintf("  In BOTH groups with >=1 real value: %d / %d", sum(coverage$in_both_groups), n_imputed),
    sprintf("Dropping %d imputed genes without real values in both groups.", length(dropped_genes)), "",
    if (is.null(validation)) "Imputation validation: not run" else
      c("Imputation validation:", capture.output(print(aggregate(correlation ~ method, validation, mean)))))
  writeLines(lines, file.path(output_dir, "summary.txt"))
}

# Copy the previous contents of an output folder to <dir>/archive/<timestamp>/
# and clear it, so a re-run never silently overwrites results.
archive_outputs <- function(output_dir) {
  old <- list.files(output_dir, full.names = TRUE)
  old <- old[!grepl("archive", old)]
  if (!length(old)) return(invisible())
  archive_dir <- file.path(output_dir, "archive", format(Sys.time(), "%Y-%m-%d_%H%M%S"))
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(old, archive_dir, recursive = TRUE)
  unlink(old, recursive = TRUE)
  cat("Archived previous results to:", archive_dir, "\n")
}

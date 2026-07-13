#!/usr/bin/env Rscript
# Figure 5: Validation plots — convergence, retention, CCC, split-half (full textwidth)
# Run from repo root: Rscript scripts/fig_validation.R

source("scripts/_common.R")

source("scripts/lib/subsampling_helpers.R")
val_config <- parse_config_arg("config/config_validation.yaml")

ref_de <- read.delim(val_config$paths$reference_de)
ref_n_deg <- sum(ref_de$adj.P.Val < 0.05 & abs(ref_de$logFC) > 1, na.rm = TRUE)
ref_n_deg_fdr <- sum(ref_de$adj.P.Val < 0.05, na.rm = TRUE)

bal_n_deg <- NA; bal_n_deg_fdr <- NA
bal_de_path <- val_config$paths$balanced_reference_de
if (!is.null(bal_de_path) && file.exists(bal_de_path)) {
  bal_de <- read.delim(bal_de_path)
  bal_n_deg <- sum(bal_de$adj.P.Val < 0.05 & abs(bal_de$logFC) > 1, na.rm = TRUE)
  bal_n_deg_fdr <- sum(bal_de$adj.P.Val < 0.05, na.rm = TRUE)
}

t1_file  <- file.path(val_dir, "test1_first_trim_subsample.tsv")
t1b_file <- file.path(val_dir, "test1b_vs_balanced.tsv")
t2_file  <- file.path(val_dir, "test2_balanced_subsample.tsv")
t3_file  <- file.path(val_dir, "test3_split_half.tsv")

has_t1  <- file.exists(t1_file)
has_t1b <- file.exists(t1b_file)
has_t2  <- file.exists(t2_file)
has_t3  <- file.exists(t3_file)

if (has_t1)  { t1  <- read.delim(t1_file,  stringsAsFactors = FALSE); t1_ok  <- t1[t1$status == "ok", ] }
if (has_t1b) { t1b <- read.delim(t1b_file, stringsAsFactors = FALSE); t1b_ok <- t1b[t1b$status == "ok", ] }
if (has_t2)  { t2  <- read.delim(t2_file,  stringsAsFactors = FALSE); t2_ok  <- t2[t2$status == "ok", ] }
if (has_t3)  { t3  <- read.delim(t3_file,  stringsAsFactors = FALSE); t3_ok  <- t3[t3$status_a == "ok" & t3$status_b == "ok", ] }

# --- 5a. Convergence (2-panel) ---
if (has_t1 && has_t1b) {
  out_file <- file.path(fig_dir, "fig_validation_convergence.png")
  png(out_file, width = TW, height = TW * 0.35, units = "in", res = DPI, pointsize = PT)
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), mgp = c(2.8, 0.8, 0))

  t1_med <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok, FUN = median)
  t1_q25 <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok,
                       FUN = function(x) quantile(x, 0.25))
  t1_q75 <- aggregate(jaccard_vs_full ~ N_1st, data = t1_ok,
                       FUN = function(x) quantile(x, 0.75))
  plot(t1_med$N_1st, t1_med$jaccard_vs_full,
       type = "b", pch = 19, ylim = c(0, 1), cex = 0.8,
       xlab = "N 1st-trim samples", ylab = "Jaccard vs full",
       main = "Jaccard convergence (IQR)",
       cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
  arrows(t1_med$N_1st, t1_q25$jaccard_vs_full,
         t1_med$N_1st, t1_q75$jaccard_vs_full,
         angle = 90, code = 3, length = 0.03, col = "grey40")

  t1_ret <- aggregate(overlap_vs_full ~ N_1st, data = t1_ok, FUN = median)
  t1b_ret <- aggregate(overlap_vs_balanced ~ N_1st, data = t1b_ok, FUN = median)
  plot(t1_ret$N_1st, t1_ret$overlap_vs_full,
       type = "b", pch = 19, col = "blue", ylim = c(0, 1), cex = 0.8,
       xlab = "N 1st-trim samples", ylab = "Median retention",
       main = "Retention convergence",
       cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
  points(t1b_ret$N_1st, t1b_ret$overlap_vs_balanced,
         type = "b", pch = 17, col = "red", cex = 0.8)
  legend("bottomright",
         legend = c(paste0("Full DEGs (", ref_n_deg, ")"),
                    paste0("Balanced DEGs (", bal_n_deg, ")")),
         col = c("blue", "red"), pch = c(19, 17), bty = "n", cex = 0.8)

  dev.off()
  cat("Saved:", out_file, "\n")
}

# --- 5b. Retention (4-panel) ---
if (has_t1 && has_t1b) {
  out_file <- file.path(fig_dir, "fig_validation_retention.png")
  png(out_file, width = TW, height = TW * 0.35, units = "in", res = DPI, pointsize = PT)
  par(mfrow = c(1, 4), mar = c(5, 4.5, 3, 0.5), mgp = c(2.5, 0.8, 0))

  boxplot(overlap_vs_full ~ N_1st, data = t1_ok,
          col = "lightblue", xlab = "N 1st-trim", ylab = "Retention",
          main = paste0("Full DEGs (", ref_n_deg, ")"), ylim = c(0, 1),
          cex.main = 0.85, cex.axis = 0.9, cex.lab = 1.0)

  boxplot(overlap_vs_balanced ~ N_1st, data = t1b_ok,
          col = "lightsalmon", xlab = "N 1st-trim", ylab = "Retention",
          main = paste0("Balanced DEGs (", bal_n_deg, ")"), ylim = c(0, 1),
          cex.main = 0.85, cex.axis = 0.9, cex.lab = 1.0)

  boxplot(overlap_fdr_only ~ N_1st, data = t1_ok,
          col = "lightblue", xlab = "N 1st-trim", ylab = "Retention",
          main = paste0("Full FDR-only (", ref_n_deg_fdr, ")"), ylim = c(0, 1),
          cex.main = 0.85, cex.axis = 0.9, cex.lab = 1.0)

  if ("overlap_fdr_vs_balanced" %in% colnames(t1b_ok)) {
    boxplot(overlap_fdr_vs_balanced ~ N_1st, data = t1b_ok,
            col = "lightsalmon", xlab = "N 1st-trim", ylab = "Retention",
            main = paste0("Bal. FDR-only (", bal_n_deg_fdr, ")"),
            ylim = c(0, 1),
            cex.main = 0.85, cex.axis = 0.9, cex.lab = 1.0)
  } else {
    plot.new(); text(0.5, 0.5, "N/A", cex = 0.9)
  }

  dev.off()
  cat("Saved:", out_file, "\n")
}

# --- 5c. CCC (2-panel) ---
if (has_t1b) {
  ccc_full_col <- if ("logfc_ccc_vs_full" %in% colnames(t1b_ok))
    "logfc_ccc_vs_full" else "logfc_r_vs_full"
  ccc_bal_col <- if ("logfc_ccc_vs_balanced" %in% colnames(t1b_ok))
    "logfc_ccc_vs_balanced" else "logfc_r_vs_balanced"
  ccc_label <- if (grepl("ccc", ccc_full_col)) "Lin's CCC" else "Pearson r"

  out_file <- file.path(fig_dir, "fig_validation_ccc.png")
  png(out_file, width = TW, height = TW * 0.4, units = "in", res = DPI, pointsize = PT)
  par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), mgp = c(2.8, 0.8, 0))

  boxplot(t1b_ok[[ccc_full_col]] ~ t1b_ok$N_1st,
          col = "lightblue", xlab = "N first-trimester samples",
          ylab = paste(ccc_label, "(logFC)"),
          main = "logFC agreement vs full", ylim = c(0, 1),
          cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)

  boxplot(t1b_ok[[ccc_bal_col]] ~ t1b_ok$N_1st,
          col = "lightsalmon", xlab = "N first-trimester samples",
          ylab = paste(ccc_label, "(logFC)"),
          main = "logFC agreement vs balanced", ylim = c(0, 1),
          cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)

  dev.off()
  cat("Saved:", out_file, "\n")
}

# --- 5d. Split-half (2x2) ---
if (has_t3) {
  out_file <- file.path(fig_dir, "fig_validation_split_half.png")
  png(out_file, width = TW, height = TW * 0.7, units = "in", res = DPI, pointsize = PT)
  par(mfrow = c(2, 2), mar = c(5, 5, 3, 1), mgp = c(2.8, 0.8, 0))

  hist(t3_ok$jaccard_a_vs_b, breaks = 20, col = "plum",
       xlab = "Jaccard (half A vs half B)",
       main = "DEG overlap between halves", xlim = c(0, 1),
       cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
  abline(v = median(t3_ok$jaccard_a_vs_b, na.rm = TRUE),
         lty = 2, col = "red", lwd = 1.5)

  if ("logfc_ccc_a_vs_b" %in% colnames(t3_ok)) {
    hist(t3_ok$logfc_ccc_a_vs_b, breaks = 20, col = "plum",
         xlab = "Lin's CCC (half A vs half B)",
         main = "logFC agreement between halves", xlim = c(0, 1),
         cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
    abline(v = median(t3_ok$logfc_ccc_a_vs_b, na.rm = TRUE),
           lty = 2, col = "red", lwd = 1.5)
  } else {
    hist(t3_ok$logfc_pearson_a_vs_b, breaks = 20, col = "plum",
         xlab = "Pearson r (half A vs half B)",
         main = "logFC correlation between halves", xlim = c(0, 1),
         cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
    abline(v = median(t3_ok$logfc_pearson_a_vs_b, na.rm = TRUE),
           lty = 2, col = "red", lwd = 1.5)
  }

  if ("overlap_a_vs_full" %in% colnames(t3_ok)) {
    boxplot(list("Half A" = t3_ok$overlap_a_vs_full,
                 "Half B" = t3_ok$overlap_b_vs_full,
                 "Both" = c(t3_ok$overlap_a_vs_full, t3_ok$overlap_b_vs_full)),
            col = c("plum", "plum", "lightyellow"),
            ylab = "Retention of full-run DEGs",
            main = "FDR+logFC retention per half", ylim = c(0, 1),
            cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
  } else {
    all_j <- c(t3_ok$jaccard_a_vs_full, t3_ok$jaccard_b_vs_full)
    hist(all_j, breaks = 20, col = "plum",
         xlab = "Jaccard (half vs full)",
         main = "Each half vs full run", xlim = c(0, 1),
         cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
    abline(v = median(all_j, na.rm = TRUE), lty = 2, col = "red", lwd = 1.5)
  }

  if ("overlap_fdr_a_vs_full" %in% colnames(t3_ok)) {
    boxplot(list("Half A" = t3_ok$overlap_fdr_a_vs_full,
                 "Half B" = t3_ok$overlap_fdr_b_vs_full,
                 "Both" = c(t3_ok$overlap_fdr_a_vs_full, t3_ok$overlap_fdr_b_vs_full)),
            col = c("plum", "plum", "lightyellow"),
            ylab = "Retention of full-run FDR-only DEGs",
            main = "FDR-only retention per half", ylim = c(0, 1),
            cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
  } else {
    boxplot(list("FDR+logFC\nA vs B" = t3_ok$jaccard_a_vs_b,
                 "FDR only\nA vs B" = t3_ok$jaccard_fdr_a_vs_b),
            col = c("lightyellow", "plum"),
            main = "FDR+logFC vs FDR-only", ylim = c(0, 1),
            ylab = "Jaccard (A vs B)",
            cex.main = 1.0, cex.axis = 0.9, cex.lab = 1.0)
  }

  dev.off()
  cat("Saved:", out_file, "\n")
}

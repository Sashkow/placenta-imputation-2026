#!/usr/bin/env Rscript
# Figure 4: ComBat sensitivity — scatter + MAE histogram (0.49 textwidth each)
# Run from repo root: Rscript scripts/fig_combat_sensitivity.R

source("scripts/_common.R")

suppressPackageStartupMessages(library(org.Hs.eg.db))

exprs_inter_post <- as.matrix(read.delim(
  file.path(combat_dir, "exprs_none_combat_ref.tsv"),
  row.names = 1, check.names = FALSE))
exprs_soft_post <- as.matrix(read.delim(
  file.path(combat_dir, "exprs_softimpute_combat_ref.tsv"),
  row.names = 1, check.names = FALSE))

shared_genes <- intersect(rownames(exprs_inter_post), rownames(exprs_soft_post))
common_samp  <- intersect(colnames(exprs_inter_post), colnames(exprs_soft_post))
inter_post <- exprs_inter_post[shared_genes, common_samp]
soft_post  <- exprs_soft_post[shared_genes, common_samp]
gene_mae <- rowMeans(abs(inter_post - soft_post))
overall_r <- cor(as.vector(inter_post), as.vector(soft_post))

sym_map <- mapIds(org.Hs.eg.db, keys = shared_genes,
                  column = "SYMBOL", keytype = "ENTREZID",
                  multiVals = "first")
outlier_symbols <- c("SPP1", "F13A1", "HBG2", "HBE1", "CCK")
outlier_entrez  <- names(sym_map)[sym_map %in% outlier_symbols]

# --- Scatter plot ---
w <- TW * 0.49
set.seed(42)
idx <- sample.int(length(inter_post), min(50000, length(inter_post)))

out_file <- file.path(fig_dir, "fig_combat_post_scatter.png")
png(out_file, width = w, height = w, units = "in", res = DPI, pointsize = PT)
par(mar = c(4.5, 4.5, 3, 0.5), mgp = c(2.5, 0.8, 0))
plot(as.vector(inter_post)[idx], as.vector(soft_post)[idx],
     pch = 16, cex = 0.15, col = rgb(0, 0, 0, 0.1),
     xlab = "Corrected expr. (intersection)",
     ylab = "Corrected expr. (softImpute)",
     main = sprintf("Post-ComBat-ref\n(%d genes, r = %.4f)",
                    length(shared_genes), overall_r),
     cex.main = 0.9, cex.axis = 0.9, cex.lab = 1.0)
abline(0, 1, col = "red", lwd = 1.5)
dev.off()
cat("Saved:", out_file, "\n")

# --- MAE histogram ---
out_file <- file.path(fig_dir, "fig_combat_gene_mae.png")
png(out_file, width = w, height = w * 0.7, units = "in", res = DPI, pointsize = PT)
par(mar = c(4.5, 4.5, 3, 0.5), mgp = c(2.5, 0.8, 0))
hist(gene_mae, breaks = 50, col = "steelblue", border = "white",
     xlab = expression(paste("Per-gene MAE (", log[2], " units)")),
     ylab = "Number of genes",
     main = sprintf("Per-gene correction difference\n(n = %d)", length(shared_genes)),
     cex.main = 0.9, cex.axis = 0.9, cex.lab = 1.0)
abline(v = median(gene_mae), col = "red", lwd = 1.5, lty = 2)
text(median(gene_mae), par("usr")[4] * 0.6,
     sprintf("median\n%.3f", median(gene_mae)),
     col = "red", cex = 0.55, pos = 4)
outlier_maes <- sort(gene_mae[outlier_entrez])
outlier_eids <- names(outlier_maes)
y_top <- par("usr")[4]
for (i in seq_along(outlier_eids)) {
  eid <- outlier_eids[i]
  abline(v = gene_mae[eid], col = "darkorange3", lwd = 1, lty = 3)
}
y_positions <- y_top * c(0.92, 0.55, 0.92, 0.55, 0.92)
for (i in seq_along(outlier_eids)) {
  eid <- outlier_eids[i]
  text(gene_mae[eid], y_positions[i],
       labels = sym_map[eid], srt = 90, adj = c(1, -0.3),
       cex = 0.55, col = "black", font = 3)
}
dev.off()
cat("Saved:", out_file, "\n")

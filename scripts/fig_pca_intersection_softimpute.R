#!/usr/bin/env Rscript
# Figure 3: Combined PCA — intersection vs softImpute with shared legend (full textwidth)
# Run from repo root: Rscript scripts/fig_pca_intersection_softimpute.R

source("scripts/_common.R")

pca_data <- function(exprs_file) {
  exprs <- read.delim(file.path(base_dir, exprs_file), row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
  pd <- pheno[pheno$sample_id %in% colnames(exprs), ]
  exprs <- exprs[, pd$sample_id, drop = FALSE]
  pc <- prcomp(t(exprs), center = TRUE, scale. = FALSE)
  var_pct <- 100 * summary(pc)$importance[2, ]
  list(pc = pc, var_pct = var_pct, group = pd$trimester,
       batch = pd$dataset_id, ngenes = nrow(exprs), nsamples = ncol(exprs))
}

pca1 <- pca_data("exprs_none_combat_ref.tsv")
pca2 <- pca_data("exprs_softimpute_combat_ref.tsv")

group_lvls <- unique(pca1$group)
grp_pal <- c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F")
grp_col <- setNames(grp_pal[seq_along(group_lvls)], group_lvls)
ds_lvls <- unique(pca1$batch)
ds_pch <- setNames(c(0:4, 6:8, 15:18)[seq_along(ds_lvls)], ds_lvls)

out_png <- file.path(fig_dir, "fig_pca_intersection_softimpute.png")
png(out_png, width = TW, height = TW * 0.45, units = "in", res = DPI, pointsize = PT)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(5, 1.5))

plot_pca_panel <- function(pcd, title_label) {
  par(mar = c(4.5, 4.5, 3, 0.5), mgp = c(2.5, 0.8, 0))
  plot(pcd$pc$x[, 1], pcd$pc$x[, 2],
       col = grp_col[pcd$group], pch = ds_pch[pcd$batch], cex = 0.9,
       xlab = sprintf("PC1 (%.1f%%)", pcd$var_pct[1]),
       ylab = sprintf("PC2 (%.1f%%)", pcd$var_pct[2]),
       main = sprintf("PCA — %s (%d genes)", title_label, pcd$ngenes),
       cex.main = 0.9, cex.axis = 0.9, cex.lab = 1.0)
}

plot_pca_panel(pca1, "intersection + ComBat-ref")
plot_pca_panel(pca2, "softImpute + ComBat-ref")

par(mar = c(0, 0, 0, 0))
plot.new()
grp_n <- table(factor(pca1$group, levels = group_lvls))
ds_n <- table(factor(pca1$batch, levels = ds_lvls))
legend("top", legend = sprintf("%s (n=%d)", group_lvls, grp_n),
       col = grp_col[group_lvls], pch = rep(15, length(group_lvls)),
       pt.cex = 1.2, cex = 0.6, bty = "n", ncol = length(group_lvls),
       xpd = TRUE)
legend("bottom", legend = sprintf("%s (n=%d)", ds_lvls, ds_n),
       col = rep("grey30", length(ds_lvls)), pch = ds_pch[ds_lvls],
       pt.cex = 1.0, cex = 0.6, bty = "n", ncol = length(ds_lvls),
       xpd = TRUE)

dev.off()
cat("Saved:", out_png, "\n")

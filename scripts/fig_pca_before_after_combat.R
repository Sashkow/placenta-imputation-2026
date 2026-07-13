#!/usr/bin/env Rscript
# Figure: PCA before/after ComBat (0.75 textwidth, ggplot2 2x2)
# Run from repo root: Rscript scripts/fig_pca_before_after_combat.R

source("scripts/_common.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

make_pca_plot <- function(exprs_file, pd, color_by, title, legend_title) {
  mat <- as.matrix(read.delim(file.path(base_dir, exprs_file),
                              row.names = 1, check.names = FALSE))
  common <- intersect(colnames(mat), pd$sample_id)
  mat <- mat[, common]
  pdd <- pd[match(common, pd$sample_id), ]
  mat <- mat[apply(mat, 1, var, na.rm = TRUE) > 0 & complete.cases(mat), ]
  pc <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_pct <- 100 * summary(pc)$importance[2, 1:2]
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                   group = as.factor(pdd[[color_by]]))
  ggplot(df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 1.2, alpha = 0.7) +
    stat_ellipse(type = "norm", level = 0.95, linetype = 2) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", var_pct[1]),
         y = sprintf("PC2 (%.1f%%)", var_pct[2]),
         color = legend_title) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.key.size = unit(2.5, "mm"),
          legend.spacing.x = unit(1, "mm"),
          plot.margin = margin(2, 4, 2, 2, "pt")) +
    guides(color = guide_legend(nrow = 3))
}

p1 <- make_pca_plot("exprs_imputed_softimpute.tsv", pheno,
                    "dataset_id", "Before ComBat (dataset)", "Dataset")
p2 <- make_pca_plot("exprs_softimpute_combat_ref.tsv", pheno,
                    "dataset_id", "After ComBat-ref (dataset)", "Dataset")
p3 <- make_pca_plot("exprs_imputed_softimpute.tsv", pheno,
                    "trimester", "Before ComBat (trimester)", "Trimester")
p4 <- make_pca_plot("exprs_softimpute_combat_ref.tsv", pheno,
                    "trimester", "After ComBat-ref (trimester)", "Trimester")

w <- TW * 0.75
out_file <- file.path(fig_dir, "fig_pca_before_after_combat.png")
png(out_file, width = w, height = w * 1.0, units = "in", res = DPI)
grid.arrange(p1, p2, p3, p4, ncol = 2,
             top = textGrob("PCA: Before/After ComBat-ref (softImpute, 17,531 genes)",
                            gp = gpar(fontsize = 11, fontface = "bold")))
dev.off()
cat("Saved:", out_file, "\n")

#!/usr/bin/env Rscript
# Figure 7: RNA-seq concordance with Prater 2021 (0.70 textwidth, ggplot2)
# Run from repo root: Rscript scripts/fig_rnaseq_concordance.R

source("scripts/_common.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
})

lin_ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  mx <- mean(x); my <- mean(y)
  sx <- var(x); sy <- var(y)
  sxy <- cov(x, y)
  2 * sxy / (sx + sy + (mx - my)^2)
}

prater <- read.xlsx(loadWorkbook("data/references/prater_2021_supp_tables.xlsx"),
                    sheet = "T1 DEGs_results_table_l2fc1")
prater$entrez <- as.character(prater$entrezgene_id)
prater <- prater[!is.na(prater$entrez) & prater$entrez != "NA", ]
prater_logfc <- setNames(as.numeric(prater$log2FoldChange), prater$entrez)
prater_padj  <- setNames(as.numeric(prater$padj), prater$entrez)
prater_sig   <- prater$entrez[!is.na(prater_padj) & prater_padj < 0.05 &
                               abs(prater_logfc[prater$entrez]) > 1]

full_de <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                      stringsAsFactors = FALSE)
sig_de  <- read.delim(file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
                      stringsAsFactors = FALSE)

ours_logfc   <- setNames(full_de$logFC, full_de$gene)
ours_sig_set <- sig_de$gene[sig_de$adj.P.Val < 0.05 & abs(sig_de$logFC) > 1]

shared <- intersect(names(ours_logfc), names(prater_logfc))
df <- data.frame(gene = shared,
                 ours_logfc = ours_logfc[shared],
                 prater_logfc = prater_logfc[shared],
                 stringsAsFactors = FALSE)
df$ours_sig   <- df$gene %in% ours_sig_set
df$prater_sig <- df$gene %in% prater_sig

r_val   <- cor(df$ours_logfc, df$prater_logfc, use = "complete.obs")
ccc_val <- lin_ccc(df$ours_logfc, df$prater_logfc)

df$category <- "Not significant in either"
df$category[df$ours_sig & df$prater_sig]  <- "Both significant"
df$category[df$ours_sig & !df$prater_sig] <- "This study only"
df$category[!df$ours_sig & df$prater_sig] <- "Prater 2021 only"
df$category <- factor(df$category,
  levels = c("Both significant", "This study only",
             "Prater 2021 only", "Not significant in either"))

cat_colors <- c("Both significant" = "#D32F2F",
                "This study only" = "#1976D2",
                "Prater 2021 only" = "#F57C00",
                "Not significant in either" = "#BDBDBD")

w <- TW * 0.70
p <- ggplot(df, aes(x = ours_logfc, y = prater_logfc, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted",
             color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 0.6) +
  scale_color_manual(values = cat_colors,
                     name = expression("Significance (FDR < 0.05, |logFC| > 1)")) +
  annotate("text", x = -3.5, y = 4.5,
           label = sprintf("r = %.3f\nCCC = %.3f\nn = %s",
                           r_val, ccc_val, format(length(shared), big.mark = ",")),
           hjust = 0, vjust = 1, size = 3.5) +
  labs(x = "logFC (this study, microarray, 117 samples)",
       y = expression("log"[2]*"FC (Prater 2021, RNA-seq, 14 samples)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(3, "mm")) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 1.5, alpha = 0.8))) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-5, 5))

out_file <- file.path(fig_dir, "fig_rnaseq_concordance_prater.png")
png(out_file, width = w, height = w * 0.95, units = "in", res = DPI, type = "cairo")
print(p)
invisible(dev.off())
cat("Saved:", out_file, "\n")

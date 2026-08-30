#!/usr/bin/env Rscript
# Figure 7: RNA-seq concordance with Prater 2021 (0.70 textwidth, ggplot2)
# Run from repo root: Rscript scripts/fig_rnaseq_concordance.R

source("scripts/_common.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
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

## ---------------------------------------------------------------------------
## Fate of every DEG called by this study.
##
## Prater's supplementary table publishes only their DEGs -- every row has
## padj < 0.05 -- so genes we call that they never listed have no y coordinate
## and cannot appear in the scatter at all. Their absence is not evidence of
## non-replication: it conflates "tested and null" with "never tested", and the
## published data cannot separate the two. Panel B keeps them visible.
## ---------------------------------------------------------------------------
n_tested     <- nrow(full_de)
fate_repl    <- sum(ours_sig_set %in% prater_sig)
fate_disagree <- sum(ours_sig_set %in% shared & !ours_sig_set %in% prater_sig)
fate_absent  <- sum(!ours_sig_set %in% shared)

stopifnot(fate_repl + fate_disagree + fate_absent == length(ours_sig_set))

expected <- c(shared = 2560, replicated = 367, disagreeing = 14, absent = 149)
observed <- c(shared = length(shared), replicated = fate_repl,
              disagreeing = fate_disagree, absent = fate_absent)
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("Figure counts diverge from the values stated in the manuscript caption.\n",
       "  expected: ", paste(names(expected), expected, sep = "=", collapse = ", "), "\n",
       "  observed: ", paste(names(observed), observed, sep = "=", collapse = ", "), "\n",
       "Update the caption in main.tex and the expectations here together, ",
       "so figure and prose cannot drift apart.")
}

## Every plotted gene satisfies FDR < 0.05 in Prater, so the residual category
## must name the fold-change criterion rather than claim non-significance.
LAB_BOTH  <- "Significant in both"
LAB_OURS  <- "This study only"
LAB_PRAT  <- "Prater 2021 only"
LAB_BELOW <- "|logFC| ≤ 1 in both"

df$category <- LAB_BELOW
df$category[df$ours_sig & df$prater_sig]  <- LAB_BOTH
df$category[df$ours_sig & !df$prater_sig] <- LAB_OURS
df$category[!df$ours_sig & df$prater_sig] <- LAB_PRAT
df$category <- factor(df$category,
  levels = c(LAB_BOTH, LAB_OURS, LAB_PRAT, LAB_BELOW))

cat_colors <- setNames(c("#D32F2F", "#1976D2", "#F57C00", "#BDBDBD"),
                       c(LAB_BOTH, LAB_OURS, LAB_PRAT, LAB_BELOW))

w <- TW * 0.70
p_scatter <- ggplot(df, aes(x = ours_logfc, y = prater_logfc, color = category)) +
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
  annotate("text", x = -3.9, y = 4.7,
           label = sprintf("r = %.3f\nCCC = %.3f\n%s of %s genes tested",
                           r_val, ccc_val,
                           format(length(shared), big.mark = ","),
                           format(n_tested, big.mark = ",")),
           hjust = 0, vjust = 1, size = 3.1) +
  labs(x = "logFC (this study, microarray)",
       y = expression("log"[2]*"FC (Prater 2021, RNA-seq)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8.5),
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = 0, b = 0),
        axis.title = element_text(size = 9.5)) +
  guides(color = guide_legend(nrow = 2, title.position = "top",
                              override.aes = list(size = 1.5, alpha = 0.8))) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = -5:5) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-5, 5))

## ---------------------------------------------------------------------------
## Panel B -- fate of all 538 DEGs, including those the scatter cannot show.
## ---------------------------------------------------------------------------
FATE_REPL   <- "Replicated in Prater"
FATE_DIS    <- "Plotted, below Prater's threshold"
FATE_ABSENT <- "Absent from Prater's published table"

fate <- data.frame(
  fate = factor(c(FATE_REPL, FATE_DIS, FATE_ABSENT),
                levels = c(FATE_ABSENT, FATE_DIS, FATE_REPL)),
  n    = c(fate_repl, fate_disagree, fate_absent),
  stringsAsFactors = FALSE
)
fate$frac <- fate$n / sum(fate$n)

## Red and blue carry the same meaning in both panels. Grey does not: in panel A
## it marks genes that ARE in Prater's table but fall below the fold-change cut,
## whereas here it would mark genes absent from that table altogether. Reusing
## the swatch for opposite conditions is misleading, so the absent band is drawn
## unfilled and hatched -- a "no data" texture rather than a measured category.
fate_colors <- setNames(c("#D32F2F", "#1976D2", "#FAFAFA"),
                        c(FATE_REPL, FATE_DIS, FATE_ABSENT))

## Diagonal hatching for the absent band, clipped to its rectangle. Hand-rolled
## because ggpattern is not a dependency of this repo. `run` is the x-distance
## spanning the band height; it is tuned so the lines read as ~45 degrees at the
## final rendered aspect, not in data units.
hatch_rect <- function(x0, x1, y0, y1, run = 32, spacing = 11) {
  m  <- (y1 - y0) / run
  cs <- seq(x0 - run, x1, by = spacing)
  out <- lapply(cs, function(cc) {
    xa <- max(cc, x0); xb <- min(cc + run, x1)
    if (xa >= xb) return(NULL)
    data.frame(x = xa, xend = xb,
               y = y0 + (xa - cc) * m, yend = y0 + (xb - cc) * m)
  })
  do.call(rbind, out)
}

BAR_LO <- 0.75
BAR_HI <- 1.25
abs_x0 <- fate_repl + fate_disagree
abs_x1 <- sum(fate$n)
hatch  <- hatch_rect(abs_x0, abs_x1, BAR_LO, BAR_HI)

## Segments are labelled in place rather than through a legend: the bar is a
## single annotation strip, and a legend would cost more height than the bar.
## The 16-gene band is ~3% of the width, too narrow to hold text, so it is
## called out above the bar with a leader line.
fate$cum      <- cumsum(fate$n) - fate$n / 2
fate$inside   <- fate$frac >= 0.10
fate$label    <- c("367 replicated", "14", "149 absent")
fate$text_col <- c("white", "grey20", "grey25")

p_fate <- ggplot(fate, aes(x = n, y = 1, fill = fate)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.4,
           orientation = "y", show.legend = FALSE) +
  geom_segment(data = hatch,
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, color = "grey55", linewidth = 0.28) +
  annotate("rect", xmin = abs_x0, xmax = abs_x1,
           ymin = BAR_LO, ymax = BAR_HI,
           fill = NA, color = "grey55", linewidth = 0.4) +
  geom_text(data = fate[fate$inside & fate$fate == FATE_REPL, ],
            aes(x = cum, label = label, colour = text_col), y = 1,
            size = 3.0, fontface = "bold") +
  ## the absent label sits over hatching, so it needs an opaque backing
  geom_label(data = fate[fate$fate == FATE_ABSENT, ],
             aes(x = cum, label = label, colour = text_col), y = 1,
             size = 3.0, fontface = "bold", fill = "#FAFAFA",
             linewidth = 0, label.padding = unit(0.6, "mm")) +
  scale_colour_identity() +
  geom_segment(data = fate[!fate$inside, ],
               aes(x = cum, xend = cum), y = 1.27, yend = 1.40,
               color = "grey50", linewidth = 0.3) +
  geom_text(data = fate[!fate$inside, ],
            aes(x = cum, label = label), y = 1.52,
            color = "grey20", size = 3.0) +
  scale_fill_manual(values = fate_colors) +
  scale_y_continuous(limits = c(0.65, 1.7), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = sprintf(paste0("Fate of all %d DEGs of this study.\n",
                          "\"Absent\" = not in Prater's published table."),
                   sum(fate$n))) +
  theme_bw(base_size = 11) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 8.5, color = "grey25"),
        panel.grid   = element_blank(),
        panel.border = element_blank(),
        plot.margin  = margin(t = 2, r = 4, b = 2, l = 4))

## ---------------------------------------------------------------------------
## Compose
## ---------------------------------------------------------------------------
p <- p_scatter / p_fate +
  plot_layout(heights = c(4.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 11, face = "bold"))

out_file <- file.path(fig_dir, "fig_rnaseq_concordance_prater.png")
png(out_file, width = w, height = w * 1.06, units = "in", res = DPI, type = "cairo")
print(p)
invisible(dev.off())
cat("Saved:", out_file, "\n")
cat(sprintf("Panel B: %d replicated / %d below threshold / %d absent = %d\n",
            fate_repl, fate_disagree, fate_absent, sum(fate$n)))

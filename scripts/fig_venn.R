#!/usr/bin/env Rscript
# Figure 6: Three-way Venn diagram (0.65 textwidth)
# Run from repo root: Rscript scripts/fig_venn.R

source("scripts/_common.R")

suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(grid))

sig_softimpute <- read.delim(
  file.path(base_dir, "difexp_significant_softimpute_combat_ref.tsv"),
  stringsAsFactors = FALSE)
sig_intersection <- read.delim(
  file.path(base_dir, "difexp_significant_none_combat_ref.tsv"),
  stringsAsFactors = FALSE)
lykhenko_sig <- read.csv(
  "data/references/lykhenko_2021_deg.csv",
  stringsAsFactors = FALSE)

set_softimpute  <- as.character(sig_softimpute$gene)
set_intersection <- as.character(sig_intersection$gene)
set_lykhenko    <- as.character(lykhenko_sig$ENTREZID)

venn_list <- setNames(
  list(set_softimpute, set_intersection, set_lykhenko),
  c(paste0("softImpute (", length(set_softimpute), ")"),
    paste0("Intersection (", length(set_intersection), ")"),
    paste0("Lykhenko 2021 (", length(set_lykhenko), ")"))
)

w <- TW * 0.65
out_file <- file.path(fig_dir, "fig_venn_three_way.png")
png(out_file, width = w, height = w * 0.85, units = "in", res = DPI)
venn.plot <- venn.diagram(
  x = venn_list, filename = NULL,
  fill = c("#2166AC", "#F4A582", "#4DAF4A"),
  alpha = 0.5,
  cex = 1.1,
  cat.cex = 0.85,
  cat.dist = 0.07,
  margin = 0.15,
  main = "DEG overlap across three analyses",
  main.cex = 1.0
)
grid.draw(venn.plot)
dev.off()
cat("Saved:", out_file, "\n")

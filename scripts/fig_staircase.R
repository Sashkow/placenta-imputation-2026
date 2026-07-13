#!/usr/bin/env Rscript
# Figure 1: NA staircase diagram (full textwidth)
# Run from repo root: Rscript scripts/fig_staircase.R

source("scripts/_common.R")

suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(org.Hs.eg.db))

Sys.setenv(STAIRCASE_COLORS = "config/staircase_colors.yaml")
source("scripts/lib/plot_na_staircase.R")

cfg <- yaml::read_yaml("config/config_pipeline.yaml")

matrices <- list()
for (ds in cfg$files$datasets) {
  fn <- if (!is.null(cfg$files$file_map[[ds]])) cfg$files$file_map[[ds]] else
    paste0(ds, cfg$files$suffix, ".tsv")
  fp <- file.path(cfg$paths$mapped_data, fn)
  m <- read.delim(fp, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  ds_samples <- pheno$sample_id[pheno$dataset_id == ds]
  ds_samples <- ds_samples[ds_samples %in% colnames(m)]
  matrices[[ds]] <- m[, ds_samples, drop = FALSE]
}

all_genes <- unique(unlist(lapply(matrices, rownames)))
all_samples <- unlist(lapply(matrices, colnames))
merged <- matrix(NA_real_, nrow = length(all_genes), ncol = length(all_samples),
                 dimnames = list(all_genes, all_samples))
for (ds in names(matrices)) {
  m <- matrices[[ds]]
  merged[rownames(m), colnames(m)] <- as.matrix(m)
}

pc_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
merged <- merged[rownames(merged) %in% pc_genes, , drop = FALSE]

sample_ds <- pheno$dataset_id[match(colnames(merged), pheno$sample_id)]
sample_group <- pheno$trimester[match(colnames(merged), pheno$sample_id)]

de_si <- read.delim(file.path(base_dir, "difexp_softimpute_combat_ref.tsv"),
                    stringsAsFactors = FALSE)
gene_fdr <- setNames(de_si$adj.P.Val, as.character(de_si$gene))

staircase <- prepare_staircase(merged, sample_group, gene_fdr,
                               cfg$phenotype$baseline, cfg$phenotype$contrast,
                               coverage_threshold = NULL)

out_png <- file.path(fig_dir, "fig_na_staircase.png")
png(out_png, width = TW, height = TW * 0.55, units = "in", res = DPI, pointsize = PT)
par(mar = c(5, 5, 3, 1))
render_staircase(staircase, "", sample_ds,
                 cex_main = 1.0, cex_legend = 0.8,
                 cex_axis = 0.9, cex_lab = 1.0, cex_ds = 0.55,
                 line_lab = 3.5, line_ds = 1.8)
dev.off()
cat("Saved:", out_png, "\n")

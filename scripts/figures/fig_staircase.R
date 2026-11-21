#!/usr/bin/env Rscript
# Figure 1: the merged gene-by-sample matrix (missingness staircase).
#
# Columns are the 117 samples grouped by dataset (first-trimester samples, then
# second-trimester ones, marked by the strip on top); rows are the 17,829
# protein-coding genes. Bottom block: genes measured in all six datasets (the
# intersection); above it: genes missing from at least one dataset, grouped by
# the set of datasets that measured them; top: genes dropped because one
# trimester has no real value. Counts on the y axis and the right-hand brackets
# are computed from the data, never typed in.
#
# Run from repo root: Rscript scripts/figures/fig_staircase.R

source("scripts/figures/_common.R")
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(org.Hs.eg.db))

cfg <- yaml::read_yaml("config/main.yaml")
ds_order <- c("GSE100051", "GSE122214", "GSE28551", "GSE37901", "GSE93520", "GSE9984")
g1 <- cfg$phenotype$baseline; g2 <- cfg$phenotype$contrast

## ---- merged matrix (protein-coding genes x samples) ------------------------
matrices <- list()
for (ds in ds_order) {
  fn <- if (!is.null(cfg$files$file_map[[ds]])) cfg$files$file_map[[ds]] else
    paste0(ds, cfg$files$suffix, ".tsv")
  m <- read.delim(file.path(cfg$paths$mapped_data, fn), row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
  s <- pheno$sample_id[pheno$dataset_id == ds]
  s <- s[s %in% colnames(m)]
  s <- s[order(pheno$trimester[match(s, pheno$sample_id)] != g1)]   # 1st trimester first
  matrices[[ds]] <- m[, s, drop = FALSE]
}
pc_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
all_genes <- intersect(unique(unlist(lapply(matrices, rownames))), pc_genes)
samples <- unlist(lapply(matrices, colnames), use.names = FALSE)
obs <- matrix(FALSE, length(all_genes), length(samples), dimnames = list(all_genes, samples))
for (ds in ds_order) {
  m <- matrices[[ds]]; g <- intersect(rownames(m), all_genes)
  obs[g, colnames(m)] <- !is.na(as.matrix(m[g, , drop = FALSE]))
}
s_ds  <- pheno$dataset_id[match(samples, pheno$sample_id)]
s_grp <- pheno$trimester[match(samples, pheno$sample_id)]

## ---- gene classes and order ------------------------------------------------
ds_obs <- sapply(ds_order, function(d) rowSums(obs[, s_ds == d, drop = FALSE]) > 0)
n_ds <- rowSums(ds_obs)
in_both <- rowSums(obs[, s_grp == g1, drop = FALSE]) > 0 & rowSums(obs[, s_grp == g2, drop = FALSE]) > 0
full <- n_ds == length(ds_order)
cls <- ifelse(full, 0L, ifelse(in_both, 1L, 2L))      # 0 intersection, 1 imputed, 2 dropped
pattern <- apply(ds_obs, 1, function(x) paste(as.integer(x), collapse = ""))
## within a class: more datasets first (sits lower), then by pattern
ord <- order(cls, -n_ds, pattern)
obs <- obs[ord, ]; cls <- cls[ord]
n_int <- sum(cls == 0); n_imp <- sum(cls == 1); n_drop <- sum(cls == 2)
n_kept <- n_int + n_imp; n_all <- length(cls)
cat(sprintf("genes: %d total, %d intersection, %d imputed, %d dropped, %d kept\n",
            n_all, n_int, n_imp, n_drop, n_kept))

## code: 1 measured, 2 missing & imputed, 3 missing & gene dropped
code <- ifelse(obs, 1L, ifelse(cls == 2, 3L, 2L))
dim(code) <- dim(obs)
col_meas <- "#9DB2C7"; col_imp <- "#FDD44F"; col_drop <- "#D32F2F"
col_t1 <- "#66C2A5"; col_t2 <- "#FC8D62"; col_hi <- "#C62828"

## ---- draw ------------------------------------------------------------------
out_png <- file.path(fig_dir, "fig_na_staircase.png")
png(out_png, width = 8, height = 6, units = "in", res = 400, pointsize = 10)
layout(matrix(1:2, 2), heights = c(0.035, 0.965))
## trimester strip
par(mar = c(0, 6, 0.4, 9.5), xaxs = "i", yaxs = "i")
image(seq_along(samples), 1, matrix(ifelse(s_grp == g1, 1, 2), ncol = 1),
      col = c(col_t1, col_t2), axes = FALSE, xlab = "", ylab = "")
ds_end <- cumsum(table(factor(s_ds, levels = ds_order)))
abline(v = ds_end[-length(ds_end)] + 0.5, col = "white", lwd = 2)
## matrix
par(mar = c(7.2, 6, 0.6, 9.5), xaxs = "i", yaxs = "i")
image(seq_along(samples), seq_len(n_all), t(code), col = c(col_meas, col_imp, col_drop),
      breaks = c(0.5, 1.5, 2.5, 3.5), axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
abline(v = ds_end[-length(ds_end)] + 0.5, col = "white", lwd = 2)
abline(h = n_int + 0.5, lty = 2, lwd = 1.2)
abline(h = n_kept + 0.5, lty = 2, lwd = 0.8, col = col_hi)
box(lwd = 0.8)
## y axis: round ticks plus the two counts in red
ticks <- seq(0, n_all, by = 5000)
axis(2, at = ticks, labels = format(ticks, big.mark = ","), las = 1, cex.axis = 0.9)
axis(2, at = c(n_int, n_kept), labels = format(c(n_int, n_kept), big.mark = ","),
     las = 1, col.axis = col_hi, font.axis = 2, cex.axis = 0.9)
mtext(sprintf("Protein-coding genes (%s)", format(n_all, big.mark = ",")), side = 2, line = 4.3)
## dataset labels, staggered
mid <- (c(0, ds_end[-length(ds_end)]) + ds_end) / 2 + 0.5
n_s <- table(factor(s_ds, levels = ds_order))
for (i in seq_along(ds_order))
  mtext(sprintf("%s (n=%d)", ds_order[i], n_s[i]), side = 1, at = mid[i],
        line = if (i %% 2) 0.8 else 2.2, cex = 0.8)
mtext(sprintf("Samples, grouped by dataset (%d in total)", length(samples)), side = 1, line = 3.7)
## right-hand brackets
x0 <- length(samples) + 3.5; xt <- length(samples) + 5
bracket <- function(lo, hi, col = "black") {
  segments(x0, lo, x0, hi, col = col, lwd = 1.2, xpd = NA)
  segments(x0 - 1.2, c(lo, hi), x0, c(lo, hi), col = col, lwd = 1.2, xpd = NA)
}
bracket(0.5, n_int + 0.5)
bracket(n_int + 0.5, n_kept + 0.5)
bracket(n_kept + 0.5, n_all + 0.5, col_hi)
text(xt, n_int / 2, sprintf("%s genes\nmeasured in all\nsix datasets\n(intersection)",
                            format(n_int, big.mark = ",")), adj = 0, xpd = NA, cex = 0.85)
text(xt, n_int + n_imp / 2, sprintf("%s genes\nmissing from at\nleast one dataset,\nimputed",
                                    format(n_imp, big.mark = ",")), adj = 0, xpd = NA, cex = 0.85)
text(xt, n_all + 0.5, sprintf("%d dropped", n_drop), adj = c(0, 0.5), xpd = NA, cex = 0.85, col = col_hi)
## legend, in its own full-width overlay below the axis title
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
leg_lab <- c("measured", "missing (imputed)", "missing (gene dropped)", "1st trimester", "2nd trimester")
legend("bottom", inset = 0.01, horiz = TRUE, bty = "n", cex = 0.8, pch = 15, pt.cex = 1.4,
       col = c(col_meas, col_imp, col_drop, col_t1, col_t2), legend = leg_lab,
       text.width = strwidth(leg_lab, cex = 0.8) + 0.02)
dev.off()
cat("Saved:", out_png, "\n")

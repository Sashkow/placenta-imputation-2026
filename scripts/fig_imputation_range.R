#!/usr/bin/env Rscript
## Supplementary Figure S1: imputation accuracy as a range plot.
##
## Values are read from data/pipeline/holdout/main/table_s1.csv produced by
## holdout_tables.R (the unified holdout run); no transcribed table is involved. One horizontal segment per imputer spans its best-to-worst
## RMSE across the three masking schemes; markers show where each scheme lands.
##
## The two dashed verticals are no-skill anchors computed from the merged matrix
## (observed cells only, after the production gene drop):
##   grand mean          RMSE of predicting one overall mean for every cell
##   per-dataset average RMSE of predicting each dataset's own mean
##
## Usage (from the repository root): Rscript scripts/fig_imputation_range.R

source(file.path({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "."
}, "holdout_common.R"))

options(bitmapType = "cairo")
TW <- 6.27; DPI <- 300; PT <- 11
fig_dir <- "article/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tbl <- read.csv(file.path(OUT_DIR, "table_s1.csv"), stringsAsFactors = FALSE)
tbl <- tbl[is.finite(tbl$rmse_mean), ]
tbl$rmse <- tbl$rmse_mean

## ---- anchors from the merged matrix --------------------------------------
D <- prepare_data(config)
X <- D$X
obs <- !is.na(X)
anchor_grand <- sd(X[obs])
res <- unlist(lapply(D$datasets, function(d) {
  v <- X[, D$ds_of_sample == d]; v <- v[!is.na(v)]; v - mean(v)
}))
anchor_dataset <- sqrt(mean(res^2))
cat(sprintf("Anchors: grand mean RMSE = %.3f, per-dataset average RMSE = %.3f\n",
            anchor_grand, anchor_dataset))
write.csv(data.frame(anchor = c("grand_mean", "per_dataset_average"),
                     rmse = c(anchor_grand, anchor_dataset), n_genes = nrow(X)),
          file.path(OUT_DIR, "fig_s1_anchors.csv"), row.names = FALSE)

## ---- layout ---------------------------------------------------------------
worst <- tbl[tbl$scheme == "progressive_tax_block", ]
block <- tbl[tbl$scheme == "gene_dataset_block", ]
if (nrow(worst) && nrow(block)) {
  methods <- worst$method[order(worst$rmse, block$rmse[match(worst$method, block$method)], worst$method)]
} else {
  methods <- unique(tbl$method[order(tbl$rmse, decreasing = FALSE)])
}
methods <- unique(c(methods, unique(tbl$method)))

mask_pch <- c(random_cells = 16, gene_dataset_block = 18, progressive_tax_block = 17)
mask_lab <- c(random_cells = "random cells",
              gene_dataset_block = "gene–dataset block",
              progressive_tax_block = "worst-case block")

col_seg  <- "grey55"; col_mark <- "black"; col_anch <- "grey40"

out_file <- file.path(fig_dir, "fig_imputation_cv_range.png")
png(out_file, width = TW, height = TW * 0.52, units = "in", res = DPI, pointsize = PT)
par(mar = c(3.6, 6.2, 2.6, 0.6), mgp = c(2.2, 0.6, 0), xaxs = "i")

n <- length(methods)
ys <- rev(seq_len(n))
xmax <- max(4, ceiling(max(c(tbl$rmse, anchor_grand)) * 1.05))
plot(NULL, xlim = c(0, xmax), ylim = c(0.5, n + 2.0), axes = FALSE, xlab = "", ylab = "")
mtext(expression("RMSE (log"[2]*" units) — lower is better"), side = 1, line = 2.2, cex = 1)
axis(1, at = 0:xmax)
for (i in seq_along(methods)) {
  axis(2, at = ys[i], labels = METHOD_LABELS[methods[i]], las = 1, tick = FALSE,
       line = -0.4, font.axis = if (methods[i] == "softimpute") 2 else 1)
}

abline(v = c(anchor_dataset, anchor_grand), lty = 2, col = col_anch)
mtext(sprintf("per-dataset\naverage (%.2f)", anchor_dataset), side = 3,
      at = anchor_dataset, line = 0.15, cex = 0.72, col = col_anch)
mtext(sprintf("grand\nmean (%.2f)", anchor_grand), side = 3,
      at = anchor_grand, line = 0.15, cex = 0.72, col = col_anch)

for (i in seq_along(methods)) {
  d <- tbl[tbl$method == methods[i], ]
  segments(min(d$rmse), ys[i], max(d$rmse), ys[i], lwd = 2.4, col = col_seg)
  points(d$rmse, rep(ys[i], nrow(d)), pch = mask_pch[d$scheme], cex = 1.25, col = col_mark)
}

legend(x = 0.24, y = n + 2.02, bty = "n", pch = mask_pch, pt.cex = 1.15, cex = 0.85,
       legend = mask_lab[names(mask_pch)])
dev.off()
cat("Written:", out_file, "\n")

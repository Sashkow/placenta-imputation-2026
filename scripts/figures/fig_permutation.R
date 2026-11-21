#!/usr/bin/env Rscript
# Supplementary figure: label-permutation null (200 shuffles of trimester + GA)
# versus the real run. Reads the output of scripts/validation/label_permutation.R.
# Run from repo root: Rscript scripts/figures/fig_permutation.R

source("scripts/figures/_common.R")

FIG <- file.path(fig_dir, "fig_validation_permutation.png")
a <- read.delim(file.path(val_dir, "permutation_null.tsv"))
a <- a[a$status == "ok", ]
ref_de <- read.delim(yaml::read_yaml("config/validation.yaml")$paths$reference_de)
OBS_DEG <- sum(ref_de$adj.P.Val < 0.05 & abs(ref_de$logFC) > 1, na.rm = TRUE)
OBS_FDR <- sum(ref_de$adj.P.Val < 0.05, na.rm = TRUE)

png(FIG, width = 2400, height = 1000, res = 300)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 1), mgp = c(2.6, 0.7, 0), las = 1, cex = 0.8)
panel <- function(x, obs, xlab, main) {
  xl <- log10(x + 1)
  on.exit(title(main, line = 1.2))
  br <- seq(0, ceiling(log10(obs + 1) * 4) / 4 + 0.25, by = 0.125)
  hist(xl, breaks = br, col = "grey75", border = "white", xaxt = "n",
       xlab = xlab, ylab = "Permutations", main = "", xlim = range(br))
  ticks <- c(0, 1, 10, 100, 1000, 10000); ticks <- ticks[log10(ticks + 1) <= max(br)]
  axis(1, at = log10(ticks + 1), labels = format(ticks, big.mark = ",", scientific = FALSE, trim = TRUE))
  abline(v = log10(obs + 1), col = "#C0392B", lwd = 2)
  mtext(sprintf("null median %s, max %s  |  real run %s",
                format(round(median(x)), big.mark = ","), format(max(x), big.mark = ","),
                format(obs, big.mark = ",")), side = 3, line = 0.2, cex = 0.75, col = "grey25")
}
panel(a$n_deg, OBS_DEG, "DEGs (FDR < 0.05 and |logFC| > 1)", "Combined threshold")
panel(a$n_fdr_only, OBS_FDR, "Genes at FDR < 0.05", "FDR alone")
dev.off()
cat("Saved:", FIG, "\n")

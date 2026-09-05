#!/usr/bin/env Rscript
#
# qPCR benchmark recovery.
#
# Counts how many of the 17 benchmark genes in
# data/references/qpcr_benchmark_genes.csv are recovered as DEGs
# (adj.P.Val < 0.05 & |logFC| > 1) by a given DE table, and reports
# directional concordance.
#
# NOTE ON THE BENCHMARK: fifteen genes come from Uuskula et al. (2012),
# where TaqMan qPCR confirmed first-vs-second-trimester differential
# expression (ANOVA and t-test, FDR P < 0.005; ten additionally peak at
# mid-gestation before declining at term). They are NOT independent of
# this cohort: Uuskula et al. is the source study of GSE37901, whose four
# second-trimester discovery samples are analysed here. IDO1 (Blaschitz
# et al. 2011, immunohistochemistry, first-trimester vs term) and GH2
# (Mannik et al. 2012, term-placenta study) are external to the cohort
# but derive from other contrasts. See the article supplement.
#
# Usage:
#   Rscript scripts/qpcr_recovery.R
#   Rscript scripts/qpcr_recovery.R --de=data/pipeline/covariate_poly2/difexp_softimpute_combat_ref.tsv

args <- commandArgs(TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}

de_path    <- get_arg("de", "data/pipeline/main/difexp_softimpute_combat_ref.tsv")
genes_path <- get_arg("genes", "data/references/qpcr_benchmark_genes.csv")
fdr        <- as.numeric(get_arg("fdr", 0.05))
lfc        <- as.numeric(get_arg("logfc", 1.0))
out_path   <- get_arg("out", "")

stopifnot(file.exists(de_path), file.exists(genes_path))

bench <- read.csv(genes_path, stringsAsFactors = FALSE)
bench$entrez_id <- as.character(bench$entrez_id)

de <- read.delim(de_path, stringsAsFactors = FALSE)
id_col <- if ("ENTREZID" %in% colnames(de)) "ENTREZID" else colnames(de)[1]
de$.id <- as.character(de[[id_col]])

res <- bench
idx <- match(res$entrez_id, de$.id)
res$tested    <- !is.na(idx)
res$logFC     <- de$logFC[idx]
res$adj.P.Val <- de$adj.P.Val[idx]
res$is_deg    <- !is.na(res$adj.P.Val) & res$adj.P.Val < fdr & abs(res$logFC) > lfc
res$direction <- ifelse(is.na(res$logFC), NA, ifelse(res$logFC > 0, "up", "down"))

n_total     <- nrow(res)
n_tested    <- sum(res$tested)
n_recovered <- sum(res$is_deg)
n_up        <- sum(res$is_deg & res$direction == "up", na.rm = TRUE)

cat("=== qPCR benchmark recovery ===\n")
cat("DE table:        ", de_path, "\n")
cat("Thresholds:      adj.P.Val <", fdr, " & |logFC| >", lfc, "\n")
cat("Benchmark genes: ", n_total, " (", n_tested, " testable in this DE table)\n", sep = "")
cat("Recovered as DEG:", sprintf("%d/%d", n_recovered, n_total), "\n")
cat("Direction:       ", n_up, "of", n_recovered, "up-regulated",
    sprintf("(%.0f%% concordant with the expected gestational increase)\n",
            100 * n_up / max(1, n_recovered)))
cat("\n")

show <- res[order(-res$is_deg, res$adj.P.Val),
            c("symbol", "entrez_id", "source_study", "logFC", "adj.P.Val", "is_deg")]
print(show, row.names = FALSE, digits = 3)

if (nzchar(out_path)) {
  write.csv(res, out_path, row.names = FALSE)
  cat("\nWritten:", out_path, "\n")
}

invisible(list(n_recovered = n_recovered, n_total = n_total))

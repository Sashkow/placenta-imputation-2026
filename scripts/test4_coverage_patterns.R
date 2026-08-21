#!/usr/bin/env Rscript
#
# Task 4.1 of the block-holdout experiment: measure the EMPIRICAL coverage
# patterns of genes that only exist thanks to imputation.
#
# For every gene in the merged matrix, record which of the six datasets
# actually measured it (a "presence pattern"). Patterns are platform-determined,
# so only a handful occur. These real patterns are what test4 copies when it
# hides values from intersection genes: masking to a realistic pattern is what
# makes the held-out gene a fair stand-in for a real gained gene.
#
# Outputs (to --out, default data/pipeline/test4/):
#   coverage_patterns_all.csv     one row per pattern over all 17,531 genes
#   coverage_patterns_gained.csv  the same restricted to the 262 gained DEGs
#   coverage_summary.txt          human-readable summary
#
# Usage: Rscript scripts/test4_coverage_patterns.R

args <- commandArgs(TRUE)
get_arg <- function(n, d) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", n, "="), "", h[1]) else d }

out_dir  <- get_arg("out", "data/pipeline/test4")
main_dir <- get_arg("main", "data/pipeline/main")
expr_dir <- get_arg("expr", "data/expression")
REF_BATCH <- get_arg("ref_batch", "GSE100051")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

imp <- as.matrix(read.delim(file.path(main_dir, "exprs_imputed_softimpute.tsv"),
                            row.names = 1, check.names = FALSE))
genes <- rownames(imp); samples <- colnames(imp)

present <- NULL; dsnames <- character(0)
for (f in sort(list.files(expr_dir, pattern = "\\.tsv$", full.names = TRUE))) {
  d <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  if (!length(intersect(colnames(d), samples))) next
  present <- cbind(present, genes %in% rownames(d))
  dsnames <- c(dsnames, sub("\\.tsv$", "", basename(f)))
}
colnames(present) <- dsnames; rownames(present) <- genes
n_ds <- length(dsnames)

pattern_of <- function(m) apply(m, 1, function(r) paste(dsnames[r], collapse = "+"))
tabulate_patterns <- function(m, label) {
  if (!nrow(m)) return(NULL)
  p <- pattern_of(m)
  tb <- as.data.frame(table(p), stringsAsFactors = FALSE)
  names(tb) <- c("pattern", "n_genes")
  tb$k <- sapply(strsplit(tb$pattern, "\\+"), length)
  tb$has_ref_batch <- grepl(REF_BATCH, tb$pattern, fixed = TRUE)
  tb$frac <- tb$n_genes / sum(tb$n_genes)
  tb$set <- label
  tb[order(-tb$n_genes), c("set", "pattern", "k", "has_ref_batch", "n_genes", "frac")]
}

k <- rowSums(present)
all_tb <- tabulate_patterns(present, "all_genes")
non_tb <- tabulate_patterns(present[k < n_ds, , drop = FALSE], "non_intersection")

id_of <- function(d) as.character(if ("ENTREZID" %in% colnames(d)) d$ENTREZID else d[[1]])
si <- read.delim(file.path(main_dir, "difexp_significant_softimpute_combat_ref.tsv"))
ni <- read.delim(file.path(main_dir, "difexp_significant_none_combat_ref.tsv"))
gained <- intersect(setdiff(id_of(si), id_of(ni)), genes)
gain_tb <- tabulate_patterns(present[gained, , drop = FALSE], "gained_degs")

write.csv(rbind(all_tb, non_tb), file.path(out_dir, "coverage_patterns_all.csv"), row.names = FALSE)
write.csv(gain_tb, file.path(out_dir, "coverage_patterns_gained.csv"), row.names = FALSE)

## ---- reference-batch enrichment test ---------------------------------------
non_idx  <- k < n_ds
gained_non <- intersect(gained, genes[non_idx])
n_pool     <- sum(non_idx)
n_pool_ref <- sum(present[non_idx, REF_BATCH])
n_g        <- length(gained_non)
n_g_ref    <- sum(present[gained_non, REF_BATCH])
ft <- fisher.test(matrix(c(n_g_ref, n_g - n_g_ref,
                           n_pool_ref - n_g_ref,
                           (n_pool - n_pool_ref) - (n_g - n_g_ref)), nrow = 2))

con <- file(file.path(out_dir, "coverage_summary.txt"), "w")
w <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
w("=== Coverage patterns in the merged matrix ===")
w("Datasets: ", paste(dsnames, collapse = ", "), "  (reference batch: ", REF_BATCH, ")")
w("Genes: ", length(genes), "   intersection (k=", n_ds, "): ", sum(k == n_ds),
  "   non-intersection: ", n_pool)
w("")
w("Genes per coverage tier k:")
for (kk in sort(unique(k))) w("  k=", kk, ": ", sum(k == kk))
w("")
w("Distinct presence patterns among non-intersection genes: ", nrow(non_tb))
w("Non-intersection genes WITHOUT the reference batch: ", n_pool - n_pool_ref,
  sprintf(" (%.1f%%)", 100 * (n_pool - n_pool_ref) / n_pool))
w("")
w("=== Gained DEGs ===")
w("Gained DEGs: ", length(gained), " (", n_g, " outside the intersection, ",
  length(gained) - n_g, " intersection genes crossing threshold)")
w("Gained DEGs carrying the reference batch: ", sum(present[gained, REF_BATCH]),
  "/", length(gained))
w("")
w("Reference-batch enrichment among non-intersection gained DEGs:")
w("  observed with ref batch: ", n_g_ref, "/", n_g)
w("  expected if independent: ", sprintf("%.1f", n_g * n_pool_ref / n_pool))
w("  Fisher exact p = ", format.pval(ft$p.value, digits = 3),
  "   odds ratio = ", sprintf("%.2f", ft$estimate))
w("")
w("Interpretation: a gene whose blocks must be anchored entirely on imputed")
w("values in the reference batch is far less likely to be called differentially")
w("expressed. This is the population test4 masks into, so replicate masks must")
w("reproduce the real pattern frequencies rather than sampling k uniformly.")
close(con)

cat("\nTop non-intersection patterns:\n"); print(head(non_tb, 12), row.names = FALSE)
cat("\nGained-DEG patterns:\n"); print(gain_tb, row.names = FALSE)
cat("\nWritten to: ", out_dir, "\n")

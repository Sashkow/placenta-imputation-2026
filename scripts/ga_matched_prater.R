#!/usr/bin/env Rscript
#
# GA-matched cross-technology concordance with Prater et al. (2021).
# Reproduces the numbers reported in Supplementary Section 5.
#
# The full pipeline compares 4-12 weeks against 14-19 weeks, while Prater
# compared 7-8 weeks against 13-14 weeks. This script summarises a re-run of
# the pipeline on a cohort trimmed to Prater's windows, isolating how much of
# the moderate genome-wide correlation is attributable to the GA-range
# mismatch rather than to the change of measurement technology.
#
# Generate the input first:
#   Rscript scripts/pipeline/run_phase2b.R --config=config/config_ga_matched_prater.yaml
#
# Usage: Rscript scripts/ga_matched_prater.R

suppressPackageStartupMessages(library(openxlsx))

args <- commandArgs(TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
de_dir   <- get_arg("dir", "data/pipeline/ga_matched_prater")
out_path <- get_arg("out", "data/pipeline/ga_matched_prater/prater_concordance.csv")

FDR <- 0.05; LFC <- 1.0

lin_ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}

de_path <- file.path(de_dir, "difexp_softimpute_combat_ref.tsv")
if (!file.exists(de_path)) {
  stop("missing GA-matched pipeline output: ", de_path,
       "\n  run: Rscript scripts/pipeline/run_phase2b.R --config=config/config_ga_matched_prater.yaml")
}

de <- read.delim(de_path, stringsAsFactors = FALSE)
id_col <- if ("ENTREZID" %in% colnames(de)) "ENTREZID" else
          if ("gene" %in% colnames(de)) "gene" else colnames(de)[1]
de$.id <- as.character(de[[id_col]])

prater <- read.xlsx(loadWorkbook("data/references/prater_2021_supp_tables.xlsx"),
                    sheet = "T1 DEGs_results_table_l2fc1")
prater$entrez <- as.character(prater$entrezgene_id)
prater <- prater[!is.na(prater$entrez) & prater$entrez != "NA", ]
prater_logfc <- setNames(as.numeric(prater$log2FoldChange), prater$entrez)
prater_padj  <- setNames(as.numeric(prater$padj), prater$entrez)
prater_sig   <- prater$entrez[!is.na(prater_padj) & prater_padj < FDR &
                                abs(prater_logfc[prater$entrez]) > LFC]

ours_logfc <- setNames(de$logFC, de$.id)
ours_sig   <- de$.id[which(de$adj.P.Val < FDR & abs(de$logFC) > LFC)]

shared <- intersect(names(ours_logfc), names(prater_logfc))
x <- ours_logfc[shared]; y <- prater_logfc[shared]

shared_sig  <- intersect(ours_sig, prater_sig)
same_dir    <- sum(sign(ours_logfc[shared_sig]) == sign(prater_logfc[shared_sig]))

res <- data.frame(
  genes_tested      = nrow(de),
  degs              = length(ours_sig),
  shared_genes      = length(shared),
  prater_r          = cor(x, y, use = "complete.obs"),
  prater_ccc        = lin_ccc(x, y),
  shared_sig        = length(shared_sig),
  same_direction    = same_dir,
  same_direction_pct = 100 * same_dir / max(1, length(shared_sig))
)

cat("=== GA-matched Prater concordance (Supplementary Section 5) ===\n")
cat("Genes tested:        ", res$genes_tested, "\n")
cat("DEGs:                ", res$degs, "\n")
cat("Shared genes:        ", res$shared_genes, "\n")
cat("Pearson r:           ", round(res$prater_r, 3), "\n")
cat("Lin's CCC:           ", round(res$prater_ccc, 3), "\n")
cat("Shared significant:  ", res$shared_sig, "\n")
cat("Directional agreement:", sprintf("%d/%d (%.1f%%)\n",
    res$same_direction, res$shared_sig, res$same_direction_pct))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_path, row.names = FALSE)
cat("\nWritten:", out_path, "\n")

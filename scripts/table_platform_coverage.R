#!/usr/bin/env Rscript
# Produces Table 6: per-dataset ENTREZID coverage
# (Platform protein-coding genes and Preprocessed gene counts)
#
# Run from repo root:
#   Rscript scripts/table_platform_coverage.R

source("scripts/_common.R")

# Platform protein-coding gene counts (deterministic, from external sources):
#   Affy HG-U133 Plus 2.0: Brainarray v25 CDF (hgu133plus2hsentrezg), 16,924 PC
#   Illumina HumanHT-12 V4.0: GPL10558 annotation, 18,420 PC
#   ABI Human Genome Survey v2: GPL2986 annotation, 16,005 PC
#   Agilent Whole Human Genome 4x44K: GPL6480 annotation, 17,504 PC
platform_pc <- c(
  GSE100051 = 18420L,
  GSE122214 = 16924L,
  GSE28551  = 16005L,
  GSE37901  = 16924L,
  GSE93520  = 17504L,
  GSE9984   = 16924L
)

platforms <- c(
  GSE100051 = "Illumina HumanHT-12 V4.0",
  GSE122214 = "Affy HG-U133 Plus 2.0",
  GSE28551  = "ABI Human Genome Survey v2",
  GSE37901  = "Affy HG-U133 Plus 2.0",
  GSE93520  = "Agilent Whole Human Genome 4x44K",
  GSE9984   = "Affy HG-U133 Plus 2.0"
)

expr_dir <- "data/expression"

cat(sprintf("%-12s  %-35s  %12s  %12s\n",
            "Dataset", "Platform", "Platform PC", "Preprocessed"))
cat(paste(rep("-", 77), collapse = ""), "\n")

for (ds in datasets_6ds) {
  expr_file <- file.path(expr_dir, paste0(ds, ".tsv"))
  preproc <- nrow(read.delim(expr_file, row.names = 1, check.names = FALSE))

  cat(sprintf("%-12s  %-35s  %12s  %12s\n",
              ds, platforms[ds],
              formatC(platform_pc[ds], format = "d", big.mark = ","),
              formatC(preproc, format = "d", big.mark = ",")))
}

cat("\nPlatform PC sources:\n")
cat("  Affy: Brainarray v25 CDF (hgu133plus2hsentrezg)\n")
cat("  Non-Affy: GPL annotation files (GPL10558, GPL2986, GPL6480)\n")

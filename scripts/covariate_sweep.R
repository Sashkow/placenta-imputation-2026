#!/usr/bin/env Rscript
#
# ComBat biological-covariate sweep -- reproduces Supplementary Tables S2 and S3.
#
# Four specifications of the biological term in the ComBat model matrix:
#   categorical  ~ trimester        + fetus_sex_estimate
#   linear       ~ ga_weeks         + fetus_sex_estimate   (primary analysis)
#   poly2        ~ poly(ga_weeks,2) + fetus_sex_estimate
#   ns3          ~ ns(ga_weeks,3)   + fetus_sex_estimate
#
# This script SUMMARISES already-computed pipeline outputs; it does not run the
# pipeline. To (re)generate the inputs:
#
#   Rscript scripts/pipeline/run_phase2b.R --config=config/config_covariate_categorical.yaml
#   Rscript scripts/pipeline/run_phase2b.R --config=config/config_pipeline.yaml            # linear (primary)
#   Rscript scripts/pipeline/run_phase2b.R --config=config/config_covariate_poly2.yaml
#   Rscript scripts/pipeline/run_phase2b.R --config=config/config_covariate_ns3.yaml
#
#   Rscript scripts/pipeline/test1_first_trim_subsample.R --config=config/config_validation_covariate_categorical.yaml
#   Rscript scripts/pipeline/test1_first_trim_subsample.R --config=config/config_validation.yaml   # linear
#   Rscript scripts/pipeline/test1_first_trim_subsample.R --config=config/config_validation_covariate_poly2.yaml
#   Rscript scripts/pipeline/test1_first_trim_subsample.R --config=config/config_validation_covariate_ns3.yaml
#
# Usage: Rscript scripts/covariate_sweep.R [--out=data/pipeline/covariate_sweep]

suppressPackageStartupMessages(library(openxlsx))

args <- commandArgs(TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
out_dir <- get_arg("out", "data/pipeline/covariate_sweep")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FDR <- 0.05; LFC <- 1.0

VARIANTS <- list(
  list(key = "categorical", label = "Categorical trimester",
       pipeline = "data/pipeline/covariate_categorical",
       validation = "data/pipeline/validation_covariate_categorical"),
  list(key = "linear",      label = "Linear GA",
       pipeline = "data/pipeline/main",
       validation = "data/pipeline/validation_covariate_linear"),
  list(key = "poly2",       label = "Poly(GA, 2)",
       pipeline = "data/pipeline/covariate_poly2",
       validation = "data/pipeline/validation_covariate_poly2"),
  list(key = "ns3",         label = "ns(GA, df=3)",
       pipeline = "data/pipeline/covariate_ns3",
       validation = "data/pipeline/validation_covariate_ns3")
)

lin_ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}

## ---- Prater reference (published significant table only) -------------------
prater <- read.xlsx(loadWorkbook("data/references/prater_2021_supp_tables.xlsx"),
                    sheet = "T1 DEGs_results_table_l2fc1")
prater$entrez <- as.character(prater$entrezgene_id)
prater <- prater[!is.na(prater$entrez) & prater$entrez != "NA", ]
prater_logfc <- setNames(as.numeric(prater$log2FoldChange), prater$entrez)
prater_padj  <- setNames(as.numeric(prater$padj), prater$entrez)
prater_sig   <- prater$entrez[!is.na(prater_padj) & prater_padj < FDR &
                                abs(prater_logfc[prater$entrez]) > LFC]

## ---- qPCR benchmark --------------------------------------------------------
bench <- read.csv("data/references/qpcr_benchmark_genes.csv", stringsAsFactors = FALSE)
bench$entrez_id <- as.character(bench$entrez_id)

id_of <- function(d) as.character(if ("ENTREZID" %in% colnames(d)) d$ENTREZID else
                                  if ("gene" %in% colnames(d)) d$gene else d[[1]])

summarise_variant <- function(v) {
  de_path <- file.path(v$pipeline, "difexp_softimpute_combat_ref.tsv")
  if (!file.exists(de_path)) {
    warning("missing pipeline output for '", v$key, "': ", de_path,
            " -- run the config listed at the top of this script")
    return(NULL)
  }
  de <- read.delim(de_path, stringsAsFactors = FALSE)
  de$.id <- id_of(de)

  is_deg <- de$adj.P.Val < FDR & abs(de$logFC) > LFC
  ours_sig <- de$.id[which(is_deg)]

  ours_logfc <- setNames(de$logFC, de$.id)
  shared <- intersect(names(ours_logfc), names(prater_logfc))

  qp <- de[match(bench$entrez_id, de$.id), ]
  qp_hit <- !is.na(qp$adj.P.Val) & qp$adj.P.Val < FDR & abs(qp$logFC) > LFC

  data.frame(
    variant        = v$key,
    label          = v$label,
    genes_tested   = nrow(de),
    fdr_only       = sum(de$adj.P.Val < FDR, na.rm = TRUE),
    degs           = sum(is_deg, na.rm = TRUE),
    prater_r       = cor(ours_logfc[shared], prater_logfc[shared], use = "complete.obs"),
    prater_ccc     = lin_ccc(ours_logfc[shared], prater_logfc[shared]),
    prater_shared_sig = sum(ours_sig %in% prater_sig),
    qpcr_recovered = sum(qp_hit),
    qpcr_total     = nrow(bench),
    stringsAsFactors = FALSE
  )
}

s2 <- do.call(rbind, lapply(VARIANTS, summarise_variant))

cat("\n=== Supplementary Table S2: ComBat covariate specifications ===\n")
print(data.frame(
  Covariate    = s2$label,
  GenesTested  = s2$genes_tested,
  FDRonly      = s2$fdr_only,
  DEGs         = s2$degs,
  PraterR      = round(s2$prater_r, 3),
  CCC          = round(s2$prater_ccc, 3),
  SharedSig    = s2$prater_shared_sig,
  qPCR         = paste0(s2$qpcr_recovered, "/", s2$qpcr_total)
), row.names = FALSE)

write.csv(s2, file.path(out_dir, "table_s2_covariate_comparison.csv"), row.names = FALSE)

## ---- Table S3: within-size pairwise Jaccard --------------------------------
jac <- NULL
for (v in VARIANTS) {
  f <- file.path(v$validation, "test1_within_size_jaccard.tsv")
  if (!file.exists(f)) {
    warning("missing within-size Jaccard for '", v$key, "': ", f)
    next
  }
  t <- read.delim(f, stringsAsFactors = FALSE)
  col <- data.frame(N_1st = t$N_1st, x = t$median_pairwise_jaccard)
  names(col)[2] <- v$key
  jac <- if (is.null(jac)) col else merge(jac, col, by = "N_1st", all = TRUE)
}

if (!is.null(jac)) {
  jac <- jac[order(jac$N_1st), ]
  cat("\n=== Supplementary Table S3: median pairwise Jaccard between",
      "independent subsamples ===\n")
  print(jac, row.names = FALSE, digits = 3)
  write.csv(jac, file.path(out_dir, "table_s3_covariate_stability.csv"), row.names = FALSE)
}

cat("\nWritten to:", out_dir, "\n")

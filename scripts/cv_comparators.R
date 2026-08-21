#!/usr/bin/env Rscript
#
# Imputation cross-validation comparator sweep -- reproduces Supplementary Table S1.
#
# Masks a fraction of the OBSERVED cells of the merged matrix, imputes them with
# each method, and compares the filled values against the withheld truth. Three
# masking strategies:
#
#   random_cells             uniformly scattered single cells (in-distribution;
#                            optimistic, because each masked cell keeps
#                            informative neighbours in its own row and column)
#   gene_dataset_block       an entire (gene, dataset) block at once -- the
#                            realistic cross-platform pattern
#   gene_dataset_block_hard  as above, but never leaving a gene below the
#                            minimum observation floor (worst case)
#
# Usage:
#   Rscript scripts/cv_comparators.R [--config=config/config_cv_comparators.yaml]
#   Rscript scripts/cv_comparators.R --methods=softimpute,knn --n_repeats=1

suppressPackageStartupMessages(library(yaml))

script_dir <- {
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "scripts"
}
source(file.path(script_dir, "pipeline", "imputation.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}

config <- read_yaml(get_arg("config", "config/config_cv_comparators.yaml"))
cv <- config$cv

methods <- if (!is.null(get_arg("methods"))) {
  strsplit(get_arg("methods"), ",")[[1]]
} else unlist(cv$methods)
mask_types <- if (!is.null(get_arg("mask_types"))) {
  strsplit(get_arg("mask_types"), ",")[[1]]
} else unlist(cv$mask_types)
n_repeats  <- as.integer(get_arg("n_repeats", cv$n_repeats %||% 3))

out_dir <- config$paths$output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- load expression matrices ---------------------------------------------
pheno <- read.delim(config$paths$phenodata, stringsAsFactors = FALSE)
cm <- config$paths$column_map
for (from in names(cm)) if (from %in% colnames(pheno)) {
  colnames(pheno)[colnames(pheno) == from] <- cm[[from]]
}
keep <- rep(TRUE, nrow(pheno))
for (col in names(config$sample_filter)) {
  keep <- keep & pheno[[col]] %in% unlist(config$sample_filter[[col]])
}
pheno <- pheno[keep, ]

exprs_list <- list()
for (ds in unlist(config$files$datasets)) {
  fn <- config$files$file_map[[ds]] %||% paste0(ds, config$files$suffix %||% "", ".tsv")
  path <- file.path(config$paths$mapped_data, fn)
  if (!file.exists(path)) { warning("missing expression file: ", path); next }
  m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  m <- m[, intersect(colnames(m), pheno$sample_id), drop = FALSE]
  if (ncol(m)) exprs_list[[ds]] <- m
}
cat("Loaded", length(exprs_list), "datasets;",
    sum(sapply(exprs_list, ncol)), "samples\n")

## ---- run the sweep ---------------------------------------------------------
all_res <- list()
for (mt in mask_types) {
  cat("\n########## mask type:", mt, "##########\n")
  res <- validate_imputation(
    exprs_list,
    min_datasets       = config$coverage$min_datasets %||% 1L,
    leave_out_fraction = cv$leave_out_fraction %||% 0.1,
    n_repeats          = n_repeats,
    methods            = methods,
    rank_max           = cv$rank_max %||% 30,
    k                  = cv$k %||% 10,
    mask_type          = mt
  )
  all_res[[mt]] <- res
}

raw <- do.call(rbind, all_res)
write.csv(raw, file.path(out_dir, "cv_comparators_raw.csv"), row.names = FALSE)

## ---- Table S1 --------------------------------------------------------------
msd <- function(x) sprintf("%.3f +/- %.3f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
tbl <- do.call(rbind, lapply(split(raw, list(raw$mask_type, raw$method), drop = TRUE),
  function(g) data.frame(
    mask_type    = g$mask_type[1],
    method       = g$method[1],
    cells_masked = round(mean(g$n_masked, na.rm = TRUE)),
    pearson_r    = msd(g$correlation),
    rmse         = msd(g$rmse),
    mae          = msd(g$mae),
    converged    = sprintf("%d/%d", sum(g$converged, na.rm = TRUE), nrow(g)),
    stringsAsFactors = FALSE)))

ord_mask <- match(tbl$mask_type, mask_types)
ord_meth <- match(tbl$method, methods)
tbl <- tbl[order(ord_mask, ord_meth), ]

cat("\n=== Supplementary Table S1: imputation accuracy by method and mask ===\n")
print(tbl, row.names = FALSE)
write.csv(tbl, file.path(out_dir, "table_s1_cv_comparators.csv"), row.names = FALSE)
cat("\nWritten to:", out_dir, "\n")

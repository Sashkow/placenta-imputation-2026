#!/usr/bin/env Rscript
#
# The integration pipeline of the article:
#
#   load datasets -> safety filter -> merge -> impute -> ComBat -> limma
#
# Each step is one function from scripts/lib/. The run is repeated for every
# enabled imputation method ("none" = only the genes all datasets share) and
# every enabled ComBat variant, so the article's four method combinations come
# out of one run.
#
# Usage (from the repository root):
#   Rscript scripts/pipeline.R --config=config/main.yaml
#   Rscript scripts/pipeline.R --config=config/main.yaml --output_dir=some/dir
#
# Previous outputs are copied to <output>/archive/<timestamp>/ first
# (--no_archive skips the copy).

library(yaml)
for (f in c("config", "load_data", "gene_coverage", "imputation", "masking", "batch_correction",
            "differential_expression", "staircase_plot", "pipeline_report"))
  source(file.path("scripts/lib", paste0(f, ".R")))

config_path <- get_arg("config", "config/main.yaml")
config      <- load_config(config_path)
output_dir  <- get_arg("output_dir", config$paths$output)
ph          <- config$phenotype                       # groups and covariates
covariates  <- ph$covariates %||% c()
fdr         <- config$thresholds$fdr %||% 0.05
logfc       <- config$thresholds$logfc %||% 1.0

if (!has_flag("no_archive")) archive_outputs(output_dir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
save_table <- function(x, name, ...) write.table(x, file.path(output_dir, name), sep = "\t", quote = FALSE, ...)

# ---- 1. Load ---------------------------------------------------------------
# One ENTREZID x sample matrix per dataset; healthy first- and second-trimester
# samples, protein-coding genes (see sample_filter / gene_filter in the config).
data       <- load_datasets(config)
exprs_list <- data$exprs_list
phenodata  <- data$phenodata

# ---- 2. Safety filter ------------------------------------------------------
# A gene missing from some datasets will be imputed there. If it has no
# measured sample at all in one trimester, that trimester's expression would
# be entirely synthetic, so the gene is removed.
coverage <- imputed_gene_group_coverage(exprs_list, phenodata, config)
if (nrow(coverage))
  write.csv(coverage, file.path(output_dir, "imputed_gene_group_coverage.csv"), row.names = FALSE)
unsafe_genes <- if (config$coverage$drop_imputed_genes_missing_in_group %||% TRUE)
  coverage$gene[!coverage$in_both_groups] else character(0)
exprs_all  <- exprs_list                              # kept for the missingness plot
exprs_list <- drop_genes(exprs_list, unsafe_genes)
cat("Safety filter removed", length(unsafe_genes), "genes\n")

# ---- 3. Merge --------------------------------------------------------------
# Union of genes (those measured in at least `min_datasets` datasets), with NA
# where a dataset does not measure a gene.
min_datasets <- config$coverage$min_datasets
incomplete   <- create_incomplete_matrix(exprs_list, min_datasets)

# Leave-out check of the imputers on this matrix (hide 10% of the measured
# cells, impute, compare) -- a quick accuracy report; the article's holdout
# validation is scripts/validation/holdout_*.R.
imputers   <- intersect(names(IMPUTERS), names(Filter(function(m) isTRUE(m$enabled), config$imputation)))
validation <- NULL
if ((config$validation$leave_out_fraction %||% 0) > 0 && length(imputers)) {
  validation <- validate_imputation(
    exprs_list, min_datasets = min_datasets,
    leave_out_fraction = config$validation$leave_out_fraction,
    n_repeats = config$validation$n_repeats %||% 5, methods = imputers,
    rank_max = config$imputation$softimpute$rank_max %||% 50,
    k = config$imputation$knn$k %||% 10,
    path_stop_step = config$imputation$softimpute$path_stop_step %||% 7L,
    mask_type = config$validation$mask_type %||% "random_cells",
    min_obs_per_gene = config$validation$min_obs_per_gene %||% 4L)
  write.csv(validation, file.path(output_dir, "imputation_validation.csv"), row.names = FALSE)
}

results <- list()
for (imputation in c(if (config$imputation$none$enabled %||% TRUE) "none", imputers)) {

  # ---- 4. Impute -----------------------------------------------------------
  if (imputation == "none") {
    merged <- merge_common_genes(exprs_list)          # genes shared by all datasets
    imputed_mask <- matrix(FALSE, nrow(merged), ncol(merged), dimnames = dimnames(merged))
  } else {
    fit <- run_imputer(imputation, incomplete, config$imputation[[imputation]])
    merged <- fit$matrix
    imputed_mask <- is.na(incomplete$matrix)          # TRUE = this cell was filled in
    if (!is.null(fit$lambda))                         # softImpute: record lambda and rank
      write.csv(data.frame(method = imputation, lambda = fit$lambda, rank = fit$rank,
                           rank_max = config$imputation[[imputation]]$rank_max %||% NA,
                           path_stop_step = config$imputation[[imputation]]$path_stop_step %||% NA),
                file.path(output_dir, paste0("imputation_fit_", imputation, ".csv")), row.names = FALSE)
  }
  save_table(merged, paste0("exprs_imputed_", imputation, ".tsv"))

  # Sample table in matrix order; samples with a missing covariate are left out
  pd   <- phenodata_for(phenodata, colnames(merged), covariates)
  keep <- if (length(covariates)) complete.cases(pd[, covariates, drop = FALSE]) else rep(TRUE, nrow(pd))
  pd <- pd[keep, ]; merged <- merged[, keep, drop = FALSE]; imputed_mask <- imputed_mask[, keep, drop = FALSE]

  for (correction in names(Filter(function(m) isTRUE(m$enabled), config$normalization))) {
    key <- paste0(imputation, "_", correction)
    cat("\n=== ", key, " ===\n", sep = "")

    # ---- 5. Batch correction ----------------------------------------------
    # ComBat removes dataset effects while keeping what the model matrix
    # describes: gestational age and fetal sex. "combat_ref" leaves the
    # reference dataset (the largest one with both trimesters) untouched.
    corrected <- combat_correct(
      merged, batch = pd[[DATASET_ID]], mod = combat_model_matrix(pd, config),
      ref_batch = if (correction == "combat_ref") config$normalization$combat_ref$ref_batch)
    save_table(corrected, paste0("exprs_", key, ".tsv"))
    plot_pca(corrected, pd, config, file.path(output_dir, paste0("pca_", key, ".png")), key)

    # ---- 6. Differential expression ---------------------------------------
    # limma, second vs first trimester, fetal sex as covariate. Imputed cells
    # get weight `de.imputed_cell_weight` (0 in the article: they do not
    # contribute to the test of the gene they belong to).
    in_comparison <- pd[[ph$group_column]] %in% c(ph$baseline, ph$contrast)
    de_exprs <- corrected[, in_comparison]
    de_mask  <- imputed_mask[rownames(de_exprs), in_comparison, drop = FALSE]
    de <- fit_limma(de_exprs, de_design(pd[in_comparison, ], config),
                    imputed_cell_weights(de_mask, config$de$imputed_cell_weight %||% 1.0))
    significant <- de[is_significant(de, fdr, logfc), ]
    save_table(de, paste0("difexp_", key, ".tsv"), row.names = FALSE)
    save_table(significant, paste0("difexp_significant_", key, ".tsv"), row.names = FALSE)
    cat(sprintf("%d genes tested, %d DEGs (%d up, %d down)\n", nrow(de), nrow(significant),
                sum(significant$logFC > 0), sum(significant$logFC < 0)))

    results[[key]] <- list(imputation = imputation, normalization = correction,
                           de = de, significant = significant)
  }
}

# ---- 7. Report -------------------------------------------------------------
comparison <- write_method_comparison(results, output_dir)
plot_staircases(create_incomplete_matrix(exprs_all, min_datasets)$matrix, exprs_list, phenodata,
                results, config, coverage_threshold = min_datasets / length(exprs_list), output_dir)
write_summary(output_dir, config_path, config, coverage, unsafe_genes, comparison, validation)
cat("\nDone. Outputs in", output_dir, "\n")

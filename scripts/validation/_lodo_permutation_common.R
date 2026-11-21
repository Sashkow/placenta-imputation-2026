# Shared setup for the leave-one-dataset-out (LODO) and label-permutation tests.
#
# Both tests reuse the validation machinery unchanged
# (scripts/lib/subsampling.R: run_lean_pipeline = drop genes with
# one all-missing trimester -> reuse the full-cohort softImpute matrix ->
# ComBat-ref with GA covariate -> limma with imputed cells at weight 0), the same
# code behind the subsampling tests.
#
# Run from repo root; sourced by leave_one_dataset_out.R and label_permutation.R.
#
# Gene universe: restricted to the 17,531 genes of the full-cohort imputed
# matrix, so every run can reuse that matrix (no re-imputation). This matches the
# "imputation run once" design of the subsampling tests and carries the same caveat.

suppressPackageStartupMessages({
  source("scripts/lib/subsampling.R")
  library(parallel)
})
config <- config_from_cli("config/validation.yaml")
OUT_DIR <- config$paths$output
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
DATA <- load_pipeline_data(config)

# restrict every dataset to the reference imputed matrix's genes
REF_GENES <- rownames(DATA$ref_imputed)
DATA$exprs_list <- lapply(DATA$exprs_list, function(m) m[rownames(m) %in% REF_GENES, , drop = FALSE])

# GA weeks fixed per sample (the ComBat biological covariate)
DATA$phenodata$ga_weeks <- compute_ga_weeks(DATA$phenodata)

run_on <- function(exprs_list, pdata) {
  run_lean_pipeline(exprs_list, pdata, config, ref_imputed = DATA$ref_imputed)
}

sample_col <- SAMPLE_ID
ds_col     <- DATASET_ID
grp_col    <- config$phenotype$group_column

# sanity: the lean pipeline on all samples must reproduce the reference list
check_full_run <- function() {
  r <- run_on(DATA$exprs_list, DATA$phenodata)
  m <- compute_subsample_metrics(r, DATA$ref_de, DATA$ref_sig_genes, DATA$ref_sig_genes_fdr_only)
  cat(sprintf("Full-cohort re-run: %d DEGs (reference %d), Jaccard %.3f, CCC %.3f\n",
              m$n_deg, length(DATA$ref_sig_genes), m$jaccard_vs_full, m$logfc_ccc_vs_full))
  invisible(m)
}

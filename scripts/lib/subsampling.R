# Re-running batch correction and differential expression on a subset of the
# samples: the engine behind the subsampling, split-half, label-permutation
# and leave-one-dataset-out tests.
#
# The steps are the main pipeline's (scripts/pipeline.R). Two things differ,
# both stated in the article: the imputed matrix of the full cohort is reused
# instead of re-imputing every subset, and a subset that lost the reference
# dataset falls back to ComBat without a reference.

library(yaml)
for (f in c("config", "load_data", "gene_coverage", "imputation", "batch_correction",
            "differential_expression", "metrics"))
  source(file.path("scripts/lib", paste0(f, ".R")))

# Datasets, sample table and the full-cohort reference results a validation
# config points to (paths.reference_de, paths.reference_imputed).
load_pipeline_data <- function(config) {
  data <- load_datasets(config)
  exprs_list <- data$exprs_list
  phenodata  <- data$phenodata
  covariates <- config$phenotype$covariates %||% c()

  # samples with a missing covariate value are left out from the start
  samples  <- unlist(lapply(exprs_list, colnames))
  pd       <- phenodata_for(phenodata, intersect(phenodata[[SAMPLE_ID]], samples), covariates)
  complete <- if (length(covariates)) complete.cases(pd[, covariates, drop = FALSE]) else rep(TRUE, nrow(pd))
  exprs_list <- keep_samples(exprs_list, pd[[SAMPLE_ID]][complete])
  phenodata  <- phenodata[phenodata[[SAMPLE_ID]] %in% unlist(lapply(exprs_list, colnames)), ]

  group <- phenodata[[config$phenotype$group_column]]
  ref_de <- read.delim(config$paths$reference_de, stringsAsFactors = FALSE)
  fdr   <- config$thresholds$fdr %||% 0.05
  logfc <- config$thresholds$logfc %||% 1.0
  ref_imputed <- NULL
  if (!is.null(config$paths$reference_imputed) && file.exists(config$paths$reference_imputed))
    ref_imputed <- as.matrix(read.delim(config$paths$reference_imputed,
                                        row.names = 1, check.names = FALSE))

  list(exprs_list = exprs_list,
       phenodata = phenodata,
       ref_de = ref_de,
       ref_sig_genes          = ref_de$gene[is_significant(ref_de, fdr, logfc)],
       ref_sig_genes_fdr_only = ref_de$gene[!is.na(ref_de$adj.P.Val) & ref_de$adj.P.Val < fdr],
       ref_sig_genes_lfc15    = ref_de$gene[is_significant(ref_de, fdr, 1.5)],
       ref_sig_genes_lfc20    = ref_de$gene[is_significant(ref_de, fdr, 2.0)],
       ref_imputed = ref_imputed,
       first_trim_samples  = phenodata[[SAMPLE_ID]][group == config$phenotype$baseline],
       second_trim_samples = phenodata[[SAMPLE_ID]][group == config$phenotype$contrast])
}

# The datasets and sample table restricted to the given samples
subsample_exprs_list <- function(exprs_list, phenodata, sample_ids) {
  list(exprs_list = keep_samples(exprs_list, sample_ids),
       phenodata  = phenodata[phenodata[[SAMPLE_ID]] %in% sample_ids, ])
}

# Safety filter -> merge -> impute (or reuse) -> ComBat -> limma on a subset.
# Returns the DE table and the DEG lists at several thresholds; a failed run
# returns status "error:<message>" instead of stopping the caller's loop.
run_lean_pipeline <- function(exprs_list, pdata, config, ref_imputed = NULL) {
  tryCatch(run_lean_pipeline_inner(exprs_list, pdata, config, ref_imputed),
           error = function(e) list(de_results = NULL, sig_genes = character(0),
                                    n_sig = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
                                    n_genes_tested = NA_integer_,
                                    status = paste0("error:", conditionMessage(e))))
}

run_lean_pipeline_inner <- function(exprs_list, pdata, config, ref_imputed = NULL) {
  ph <- config$phenotype

  # 1. Safety filter and merge
  if (isTRUE(config$coverage$drop_imputed_genes_missing_in_group)) {
    coverage   <- imputed_gene_group_coverage(exprs_list, pdata, config)
    exprs_list <- drop_genes(exprs_list, coverage$gene[!coverage$in_both_groups])
  }
  incomplete   <- create_incomplete_matrix(exprs_list, min_datasets = config$coverage$min_datasets %||% 1L)
  imputed_mask <- is.na(incomplete$matrix)

  # 2. Impute -- or take the cells from the full-cohort imputed matrix
  genes <- rownames(incomplete$matrix); samples <- colnames(incomplete$matrix)
  reuse <- !is.null(ref_imputed) &&
    all(samples %in% colnames(ref_imputed)) && all(genes %in% rownames(ref_imputed))
  merged <- if (reuse) ref_imputed[genes, samples]
            else run_imputer("softimpute", incomplete, config$imputation$softimpute)$matrix

  # 3. Batch handling: ComBat, or (diagnostic) dataset as a term of the limma model
  pd    <- phenodata_for(pdata, samples)
  batch <- as.factor(pd[[DATASET_ID]])
  covariates <- informative_covariates(pd, ph$covariates %||% c())
  mod   <- combat_model_matrix(pd, config, covariates)
  batch_in_limma <- identical(config$normalization$method, "batch_in_limma")
  ref_batch <- config$normalization$ref_batch
  corrected <-
    if (batch_in_limma) merged
    else if (!is.null(ref_batch) && ref_batch %in% levels(batch))
      tryCatch(combat_correct(merged, batch, mod, ref_batch = ref_batch),
               error = function(e) combat_correct(merged, batch, mod))
    else combat_correct(merged, batch, mod)

  # 4. limma, imputed cells down-weighted
  in_comparison <- pd[[ph$group_column]] %in% c(ph$baseline, ph$contrast)
  de_exprs <- corrected[, in_comparison]
  de_mask  <- imputed_mask[rownames(de_exprs), in_comparison, drop = FALSE]
  design   <- de_design(pd[in_comparison, ], config, covariates, batch_in_design = batch_in_limma)
  de <- fit_limma(de_exprs, design,
                  imputed_cell_weights(de_mask, config$de$imputed_cell_weight %||% 1.0))

  fdr   <- config$thresholds$fdr %||% 0.05
  logfc <- config$thresholds$logfc %||% 1.0
  sig   <- is_significant(de, fdr, logfc)
  fdr_only <- !is.na(de$adj.P.Val) & de$adj.P.Val < fdr
  list(de_results = de,
       sig_genes          = de$gene[sig],
       sig_genes_fdr_only = de$gene[fdr_only],
       sig_genes_lfc15    = de$gene[is_significant(de, fdr, 1.5)],
       sig_genes_lfc20    = de$gene[is_significant(de, fdr, 2.0)],
       n_sig = sum(sig), n_sig_fdr_only = sum(fdr_only),
       n_sig_lfc15 = sum(is_significant(de, fdr, 1.5)),
       n_sig_lfc20 = sum(is_significant(de, fdr, 2.0)),
       n_up = sum(de$logFC[sig] > 0), n_down = sum(de$logFC[sig] < 0),
       n_genes_tested = nrow(de), status = "ok")
}

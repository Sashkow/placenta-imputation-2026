## Shared setup for the holdout validation scripts (holdout_*.R) and the
## lambda / ncp selection scripts. Run everything from the repository root.
##
## Loads the config, prepares the merged matrix exactly as the main pipeline
## does, and wraps the shared ComBat and limma steps (scripts/lib/) for use on
## masked copies of that matrix.

source("scripts/lib/subsampling.R")   # loads the rest of scripts/lib/

## a sourcing script may set HOLDOUT_DEFAULT_CONFIG to change the config used
## when --config is not given
config <- config_from_cli(if (exists("HOLDOUT_DEFAULT_CONFIG")) HOLDOUT_DEFAULT_CONFIG
                          else "config/holdout.yaml")
OUT_DIR <- get_arg("out", config$paths$output %||% "output")
for (d in c("", "masks", "checkpoints", "reference", "values", "de"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

HCFG      <- config$holdout
REF_BATCH <- HCFG$ref_batch %||% "GSE100051"
FDR       <- config$thresholds$fdr   %||% 0.05
LFC       <- config$thresholds$logfc %||% 1.0
METHODS    <- split_arg("methods", unlist(HCFG$methods))
MASK_TYPES <- split_arg("schemes", unlist(HCFG$mask_types))
N_REPEATS  <- as.integer(get_arg("n_repeats", HCFG$n_repeats %||% 3))
LEAVE_OUT  <- as.numeric(get_arg("frac", HCFG$leave_out_fraction %||% 0.1))
FORCE      <- has_flag("force")

## R's yaml parser reads "1e-4" as a STRING (YAML 1.1 floats need a dot), and
## softImpute's `while (ratio > thresh)` then compares as strings and stops
## after one ALS iteration. Every numeric-looking string is coerced here so the
## imputers get the numbers the config intends.
coerce_numeric <- function(cfg) {
  lapply(cfg, function(v) {
    if (is.character(v) && length(v) == 1 && !is.na(suppressWarnings(as.numeric(v))))
      as.numeric(v) else v
  })
}
## --emulate_production_thresh (or holdout.emulate_production_thresh: true) leaves
## the strings uncoerced so softImpute behaves exactly as the shipped pipeline did.
EMULATE_PROD_THRESH <- has_flag("emulate_production_thresh") ||
  isTRUE(config$holdout$emulate_production_thresh)
method_cfg <- function(method) {
  cfg <- config$imputation[[method]] %||% list()
  if (EMULATE_PROD_THRESH) cfg else coerce_numeric(cfg)
}

sig_call <- function(l, p) !is.na(p) & !is.na(l) & p < FDR & abs(l) > LFC

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

## ------------------------------------------------------------------ data
## Merged matrix prepared as the production pipeline prepares it: sample
## filter, covariate completeness, min_datasets, and the drop of genes with no
## real measurement in one comparison group (applied BEFORE imputation, as
## scripts/pipeline.R does).
prepare_data <- function(config) {
  data <- load_pipeline_data(config)
  exprs_list <- data$exprs_list
  pdata <- data$phenodata

  incomplete <- create_incomplete_matrix(exprs_list,
                                         min_datasets = config$coverage$min_datasets %||% 1L)
  X <- incomplete$matrix
  samples <- colnames(X)
  orig_mask <- is.na(X)

  id_col <- SAMPLE_ID
  ds_of_sample  <- setNames(as.character(pdata$secondaryaccession), pdata[[id_col]])[samples]
  grp_col <- config$phenotype$group_column
  grp_of_sample <- setNames(as.character(pdata[[grp_col]]), pdata[[id_col]])[samples]

  dropped <- character(0)
  if (config$coverage$drop_imputed_genes_missing_in_group %||% TRUE) {
    real_base <- rowSums(!orig_mask[, grp_of_sample == config$phenotype$baseline, drop = FALSE])
    real_cont <- rowSums(!orig_mask[, grp_of_sample == config$phenotype$contrast, drop = FALSE])
    bad <- real_base == 0 | real_cont == 0
    if (any(bad)) {
      dropped <- rownames(X)[bad]
      log_msg("Dropping ", length(dropped),
              " genes with no real values in one group (as the production pipeline does)")
      X <- X[!bad, , drop = FALSE]
      orig_mask <- orig_mask[!bad, , drop = FALSE]
      incomplete$matrix <- X
      if (!is.null(incomplete$gene_info))
        incomplete$gene_info <- incomplete$gene_info[!bad, , drop = FALSE]
    }
  }
  genes <- rownames(X)

  pd <- pdata[match(samples, pdata[[id_col]]), ]
  rownames(pd) <- samples

  datasets <- names(exprs_list)
  ds_cols <- split(seq_len(ncol(X)), factor(ds_of_sample, levels = datasets))
  presence <- sapply(datasets, function(ds)
    rowSums(!is.na(X[, ds_cols[[ds]], drop = FALSE])) > 0)
  dimnames(presence) <- list(genes, datasets)

  log_msg("Merged matrix: ", nrow(X), " genes x ", ncol(X), " samples; ",
          sprintf("%.1f%% missing", 100 * mean(orig_mask)))

  list(X = X, orig_mask = orig_mask, incomplete = incomplete, pd = pd,
       genes = genes, samples = samples, datasets = datasets,
       ds_of_sample = ds_of_sample, grp_of_sample = grp_of_sample,
       ds_cols = ds_cols, presence = presence,
       n_ds_per_gene = rowSums(presence), dropped_genes = dropped)
}

## ------------------------------------------------------ per-gene visibility
## For a new mask (logical genes x samples), what each gene is left with:
## number of masked cells, visible datasets (>= 1 real cell after masking),
## whether the reference batch is among them, and whether a comparison group
## has no real cell left (production drops such genes before testing).
gene_visibility <- function(D, new_mask, config) {
  real_after <- !(D$orig_mask | new_mask)
  vis <- sapply(D$datasets, function(ds)
    rowSums(real_after[, D$ds_cols[[ds]], drop = FALSE]) > 0)
  dimnames(vis) <- list(D$genes, D$datasets)
  base_ok <- rowSums(real_after[, D$grp_of_sample == config$phenotype$baseline, drop = FALSE]) > 0
  cont_ok <- rowSums(real_after[, D$grp_of_sample == config$phenotype$contrast, drop = FALSE]) > 0
  data.frame(gene = D$genes,
             n_masked_cells = rowSums(new_mask),
             k_visible = rowSums(vis),
             ref_batch_visible = vis[, REF_BATCH],
             lost_group = !(base_ok & cont_ok),
             stringsAsFactors = FALSE, row.names = NULL)
}

## --------------------------------------------------------- ComBat + limma
## The shared steps, bound to this cohort. ComBat errors propagate (no silent
## fallback from ComBat-ref to plain ComBat) so the caller can record them.
##   combat_correct(mat)            ComBat-ref, or plain ComBat with
##                                  normalization.combat_variant: "plain"
##   de_fit(corrected, weight_mask) limma with the masked cells at
##                                  de.imputed_cell_weight
make_de_tools <- function(D, config) {
  ph <- config$phenotype
  covariates <- ph$covariates %||% character(0)
  pd <- D$pd
  keep_samples <- pd[[ph$group_column]] %in% c(ph$baseline, ph$contrast)

  correct <- function(mat) {
    pd_mat <- pd[colnames(mat), , drop = FALSE]
    mod <- combat_model_matrix(pd_mat, config, informative_covariates(pd_mat, covariates))
    plain <- identical(config$normalization$combat_variant %||% "ref", "plain")
    combat_correct(mat, pd_mat[[DATASET_ID]], mod,
                   ref_batch = if (!plain) config$normalization$combat_ref$ref_batch %||% REF_BATCH)
  }

  de_fit <- function(corrected, weight_mask) {
    E <- corrected[, keep_samples, drop = FALSE]
    pd_de <- pd[keep_samples, , drop = FALSE]
    covs <- informative_covariates(pd_de, covariates)
    ok <- stats::complete.cases(pd_de[, c(ph$group_column, covs), drop = FALSE])
    E <- E[, ok, drop = FALSE]
    W <- matrix(1.0, nrow(E), ncol(E), dimnames = dimnames(E))
    W[weight_mask[rownames(E), colnames(E), drop = FALSE]] <- config$de$imputed_cell_weight %||% 0.0
    de <- suppressWarnings(fit_limma(E, de_design(pd_de[ok, , drop = FALSE], config, covs), W))
    data.frame(gene = de$gene, logFC = de$logFC, adj.P.Val = de$adj.P.Val,
               stringsAsFactors = FALSE, row.names = NULL)
  }

  list(de_fit = de_fit, combat_correct = correct, keep_samples = keep_samples)
}

## ------------------------------------------------------------- reference
reference_paths <- function(method) list(
  de        = file.path(OUT_DIR, "reference", paste0(method, "_difexp.csv")),
  corrected = file.path(OUT_DIR, "reference", paste0(method, "_corrected.rds")),
  meta      = file.path(OUT_DIR, "reference", paste0(method, "_meta.csv"))
)

load_reference <- function(method) {
  p <- reference_paths(method)
  if (!file.exists(p$de) || !file.exists(p$corrected))
    stop("Reference run for '", method, "' not found in ", file.path(OUT_DIR, "reference"),
         ". Run scripts/validation/holdout_reference_runs.R first.")
  list(de = read.csv(p$de, stringsAsFactors = FALSE),
       corrected = readRDS(p$corrected))
}

## Pretty names used in tables and figures
METHOD_LABELS <- c(softimpute = "softImpute", missmda = "imputePCA", knn = "KNN",
                   knn_centred = "KNN (centered)", sample_knn = "sample KNN",
                   sample_knn_centred = "sample KNN (centered)",
                   biscale_zero = "biScale fill", gene_mean = "gene mean", batch_mean = "batch mean")
SCHEME_LABELS <- c(random_cells = "random cells",
                   gene_dataset_block = "gene--dataset block",
                   progressive_tax_block = "worst-case block")

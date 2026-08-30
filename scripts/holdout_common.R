## holdout_common.R -- shared setup for the unified holdout validation.
##
## Sourced by the other holdout_*.R scripts. Run everything from the repository root.
## Loads the config, resolves paths, sources the companion pipeline helpers
## (which hard-code paths relative to the companion root, so the working
## directory is switched there), loads and prepares the merged matrix exactly
## as the production pipeline does, and provides the ComBat / limma tools
## lifted from the companion's test4_block_holdout_de.R.

suppressPackageStartupMessages({
  library(yaml)
  library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

.args <- commandArgs(TRUE)
get_arg <- function(n, d = NULL) {
  h <- grep(paste0("^--", n, "="), .args, value = TRUE)
  if (length(h)) sub(paste0("^--", n, "="), "", h[1]) else d
}
has_flag <- function(n) any(.args == paste0("--", n))
split_arg <- function(n, d) {
  v <- get_arg(n)
  if (is.null(v)) d else strsplit(v, ",")[[1]]
}

script_dir <- {
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else getwd()
}
ALT_DIR <- normalizePath(file.path(script_dir, ".."))

config <- yaml::read_yaml(get_arg("config",
  file.path(ALT_DIR, "config", "config_holdout_unified.yaml")))

COMPANION <- normalizePath(get_arg("companion", config$paths$companion))
OUT_DIR <- {
  p <- get_arg("out", config$paths$output %||% "output")
  if (!grepl("^/", p)) p <- file.path(ALT_DIR, p)
  p
}
for (d in c("", "masks", "checkpoints", "reference", "values", "de"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

## The companion helpers source each other via "scripts/pipeline/..." and the
## config's data paths are relative to the companion root.
setwd(COMPANION)
source("scripts/pipeline/subsampling_helpers.R")  # also sources imputation.R, normalization.R

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
## run_phase2b.R does).
prepare_data <- function(config) {
  data <- load_pipeline_data(config)
  exprs_list <- data$exprs_list
  pdata <- data$phenodata

  incomplete <- create_incomplete_matrix(exprs_list,
                                         min_datasets = config$coverage$min_datasets %||% 1L)
  X <- incomplete$matrix
  samples <- colnames(X)
  orig_mask <- is.na(X)

  id_col <- "arraydatafile_exprscolumnnames"
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
## Lifted from the companion's test4_block_holdout_de.R. Differences: no
## silent fallback from ComBat-ref to plain ComBat (errors propagate so the
## caller can record them), and de_fit() also returns whether the estimate is
## finite.
make_de_tools <- function(D, config) {
  group_col  <- config$phenotype$group_column
  baseline   <- config$phenotype$baseline
  contrast   <- config$phenotype$contrast
  covariates <- config$phenotype$covariates %||% character(0)
  pd <- D$pd
  keep_samples <- pd[[group_col]] %in% c(baseline, contrast)

  build_design <- function(pdsub) {
    grp <- factor(pdsub[[group_col]], levels = c(baseline, contrast))
    dd <- data.frame(group = grp, stringsAsFactors = FALSE)
    terms <- "group"
    for (cov in covariates) {
      if (length(unique(na.omit(pdsub[[cov]]))) >= 2) {
        dd[[cov]] <- pdsub[[cov]]; terms <- c(terms, cov)
      }
    }
    ok <- stats::complete.cases(dd)
    list(design = model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))),
                               data = dd[ok, , drop = FALSE]),
         keep = ok)
  }

  de_fit <- function(corrected, weight_mask) {
    E <- corrected[, keep_samples, drop = FALSE]
    W_mask <- weight_mask[rownames(E), colnames(E), drop = FALSE]
    bd <- build_design(pd[keep_samples, , drop = FALSE])
    E <- E[, bd$keep, drop = FALSE]
    W_mask <- W_mask[, bd$keep, drop = FALSE]
    W <- matrix(1.0, nrow(E), ncol(E), dimnames = dimnames(E))
    W[W_mask] <- config$de$imputed_cell_weight %||% 0.0
    fit <- suppressWarnings(eBayes(lmFit(E, bd$design, weights = W)))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.frame(gene = rownames(tt), logFC = tt$logFC, adj.P.Val = tt$adj.P.Val,
               stringsAsFactors = FALSE, row.names = NULL)
  }

  combat_correct <- function(mat) {
    merged_pdata <- pd[colnames(mat), , drop = FALSE]
    batch <- as.factor(merged_pdata$secondaryaccession)
    cbc <- config$phenotype$combat_bio_covariate
    if (!is.null(cbc)) {
      if (!cbc %in% colnames(merged_pdata))
        merged_pdata[[cbc]] <- compute_ga_weeks(merged_pdata)
      md <- data.frame(ga_linear = merged_pdata[[cbc]])
      terms <- "ga_linear"
    } else {
      md <- data.frame(bio_group = merged_pdata[[group_col]], stringsAsFactors = FALSE)
      terms <- "bio_group"
    }
    for (cov in covariates)
      if (length(unique(na.omit(merged_pdata[[cov]]))) >= 2) {
        md[[cov]] <- merged_pdata[[cov]]; terms <- c(terms, cov)
      }
    mod <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = md)
    ## normalization.combat_variant: "ref" (reference-batch ComBat, production) or
    ## "plain" (ComBat without a reference batch) -- same model matrix either way
    variant <- config$normalization$combat_variant %||% "ref"
    if (identical(variant, "plain")) return(normalize_combat(mat, batch, mod = mod))
    rb <- config$normalization$combat_ref$ref_batch %||% REF_BATCH
    normalize_combat_ref(mat, batch, mod = mod, ref_batch = rb)
  }

  list(de_fit = de_fit, combat_correct = combat_correct, keep_samples = keep_samples)
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
         ". Run scripts/holdout_reference_runs.R first.")
  list(de = read.csv(p$de, stringsAsFactors = FALSE),
       corrected = readRDS(p$corrected))
}

## Pretty names used in tables and figures
METHOD_LABELS <- c(softimpute = "softImpute", missmda = "imputePCA", knn = "KNN",
                   gene_mean = "gene mean", batch_mean = "batch mean")
SCHEME_LABELS <- c(random_cells = "random cells",
                   gene_dataset_block = "gene--dataset block",
                   progressive_tax_block = "worst-case block")

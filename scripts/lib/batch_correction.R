# Batch correction with ComBat (sva): the model matrix that protects the
# biological signal, and the correction itself.

library(sva)

# Covariates that take at least two values in these samples
informative_covariates <- function(pdata, covariates) {
  Filter(function(cov) length(unique(na.omit(pdata[[cov]]))) >= 2, covariates)
}

# Model matrix of the variation ComBat must preserve.
#
# The biological term is gestational age in weeks when
# config$phenotype$combat_bio_covariate is set -- linear by default, a
# polynomial (combat_bio_degree) or natural spline (combat_bio_spline_df) for
# the covariate-sensitivity runs -- and the comparison group otherwise.
# `covariates` (e.g. fetal sex) are added as they are.
combat_model_matrix <- function(pdata, config, covariates = config$phenotype$covariates %||% c()) {
  ph <- config$phenotype
  if (is.null(ph$combat_bio_covariate)) {
    model_data <- data.frame(bio_group = pdata[[ph$group_column]], stringsAsFactors = FALSE)
  } else {
    ga <- pdata[[ph$combat_bio_covariate]] %||% compute_ga_weeks(pdata, ph$group_column)
    basis <- if (!is.null(ph$combat_bio_spline_df)) splines::ns(ga, df = ph$combat_bio_spline_df)
             else if ((ph$combat_bio_degree %||% 1) > 1) poly(ga, degree = ph$combat_bio_degree)
             else NULL
    prefix <- if (!is.null(ph$combat_bio_spline_df)) "ga_ns" else "ga_poly"
    if (is.null(basis)) {
      model_data <- data.frame(ga_linear = ga)
    } else {
      model_data <- as.data.frame(unclass(basis)[, seq_len(ncol(basis)), drop = FALSE])
      colnames(model_data) <- paste0(prefix, seq_len(ncol(basis)))
    }
  }
  for (cov in covariates) model_data[[cov]] <- pdata[[cov]]
  model.matrix(as.formula(paste("~", paste(colnames(model_data), collapse = " + "))),
               data = model_data)
}

# ComBat (parametric empirical Bayes). With `ref_batch`, that dataset is left
# untouched and the others are adjusted towards it (reference-batch ComBat);
# without it, all datasets are adjusted towards the grand mean.
# Genes with zero variance are removed first (ComBat cannot scale them).
combat_correct <- function(exprs, batch, mod = NULL, ref_batch = NULL) {
  batch <- as.factor(batch)
  if (!is.null(ref_batch) && !ref_batch %in% levels(batch))
    stop("ref_batch '", ref_batch, "' not found in batch levels: ",
         paste(levels(batch), collapse = ", "))
  gene_var <- apply(exprs, 1, var, na.rm = TRUE)
  exprs <- exprs[gene_var > 0 & !is.na(gene_var), ]
  ComBat(dat = as.matrix(exprs), batch = batch, mod = mod,
         par.prior = TRUE, prior.plots = FALSE, ref.batch = ref_batch)
}

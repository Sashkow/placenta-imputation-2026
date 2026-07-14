#' Normalization Module for Phase 2B
#'
#' Functions for cross-platform/batch normalization of gene expression data.
#' Implements ComBat (baseline), DWD, and other alternatives.
#'
#' @author Expression Integration Pipeline
#' @date 2026-01

library(sva)

#' Apply ComBat batch correction (baseline)
#'
#' Standard parametric empirical Bayes batch correction.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector (factor or character)
#' @param mod Model matrix for biological covariates (optional)
#' @param par_prior Use parametric priors (default: TRUE)
#' @return Batch-corrected expression matrix
#' @export
normalize_combat <- function(exprs, batch, mod = NULL, par_prior = TRUE) {

  cat("\n=== ComBat Batch Correction ===\n")

  batch <- as.factor(batch)
  cat("Batches:", nlevels(batch), "\n")
  cat("Batch sizes:", paste(table(batch), collapse = ", "), "\n")

  # Remove genes with zero variance
  gene_vars <- apply(exprs, 1, var, na.rm = TRUE)
  keep <- gene_vars > 0 & !is.na(gene_vars)
  if (sum(!keep) > 0) {
    cat("Removing", sum(!keep), "genes with zero variance\n")
    exprs <- exprs[keep, ]
  }

  # Apply ComBat
  corrected <- ComBat(
    dat = as.matrix(exprs),
    batch = batch,
    mod = mod,
    par.prior = par_prior,
    prior.plots = FALSE
  )

  cat("ComBat correction complete\n")

  corrected
}


#' Reference-batch ComBat
#'
#' Standard ComBat adjusts all batches toward the grand mean, which mixes
#' batch and biology when datasets are unbalanced. Reference-batch ComBat
#' leaves one dataset untouched and adjusts all others toward it.
#' Choosing a balanced dataset as the reference anchors the correction
#' to a batch where batch and biology are decoupled.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @param mod Model matrix for biological covariates
#' @param ref_batch Name of the reference batch (must be a level of batch)
#' @param par_prior Use parametric priors (default: TRUE)
#' @return Batch-corrected expression matrix
normalize_combat_ref <- function(exprs, batch, mod = NULL, ref_batch,
                                 par_prior = TRUE) {

  cat("\n=== Reference-Batch ComBat ===\n")

  batch <- as.factor(batch)
  cat("Batches:", nlevels(batch), "\n")
  cat("Reference batch:", ref_batch, "\n")

  if (!ref_batch %in% levels(batch))
    stop("ref_batch '", ref_batch, "' not found in batch levels: ",
         paste(levels(batch), collapse = ", "))

  ref_n <- sum(batch == ref_batch)
  cat(sprintf("Reference batch samples: %d / %d total\n", ref_n, length(batch)))
  cat("Batch sizes:", paste(names(table(batch)), table(batch),
                             sep = "=", collapse = ", "), "\n")

  gene_vars <- apply(exprs, 1, var, na.rm = TRUE)
  keep <- gene_vars > 0 & !is.na(gene_vars)
  if (sum(!keep) > 0) {
    cat("Removing", sum(!keep), "genes with zero variance\n")
    exprs <- exprs[keep, ]
  }

  corrected <- ComBat(
    dat = as.matrix(exprs),
    batch = batch,
    mod = mod,
    par.prior = par_prior,
    prior.plots = FALSE,
    ref.batch = ref_batch
  )

  cat("Reference-batch ComBat correction complete\n")

  corrected
}


#' Apply Distance Weighted Discrimination (DWD) normalization
#'
#' DWD is a margin-based method that works better than ComBat when:
#' - Group sizes are unequal
#' - Strong batch effects exist
#' - Data is high-dimensional relative to sample size
#'
#' This implementation uses a simplified DWD approach based on
#' mean-shift correction with distance weighting.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @param mod Model matrix for biological covariates (optional)
#' @return Batch-corrected expression matrix
#' @export
normalize_dwd <- function(exprs, batch, mod = NULL) {

  cat("\n=== DWD-style Batch Correction ===\n")

  batch <- as.factor(batch)
  batches <- levels(batch)
  n_batches <- length(batches)

  cat("Batches:", n_batches, "\n")
  cat("Batch sizes:", paste(table(batch), collapse = ", "), "\n")

  # Reference batch (largest)
  batch_sizes <- table(batch)
  ref_batch <- names(which.max(batch_sizes))
  cat("Reference batch:", ref_batch, "\n")

  exprs <- as.matrix(exprs)
  corrected <- exprs

  # Calculate reference batch statistics
  ref_idx <- which(batch == ref_batch)
  ref_mean <- rowMeans(exprs[, ref_idx, drop = FALSE], na.rm = TRUE)
  ref_sd <- apply(exprs[, ref_idx, drop = FALSE], 1, sd, na.rm = TRUE)
  ref_sd[ref_sd == 0] <- 1  # Avoid division by zero

  # Correct each non-reference batch
  for (b in setdiff(batches, ref_batch)) {
    batch_idx <- which(batch == b)
    cat("  Correcting batch", b, "(", length(batch_idx), "samples)...\n")

    # Calculate batch statistics
    batch_mean <- rowMeans(exprs[, batch_idx, drop = FALSE], na.rm = TRUE)
    batch_sd <- apply(exprs[, batch_idx, drop = FALSE], 1, sd, na.rm = TRUE)
    batch_sd[batch_sd == 0] <- 1

    # DWD-style correction: scale then shift
    # 1. Z-score within batch
    z_batch <- sweep(exprs[, batch_idx, drop = FALSE], 1, batch_mean, "-")
    z_batch <- sweep(z_batch, 1, batch_sd, "/")

    # 2. Rescale to reference distribution
    corrected[, batch_idx] <- sweep(z_batch, 1, ref_sd, "*")
    corrected[, batch_idx] <- sweep(corrected[, batch_idx], 1, ref_mean, "+")
  }

  # Verify correction
  for (b in batches) {
    batch_idx <- which(batch == b)
    new_mean <- mean(rowMeans(corrected[, batch_idx, drop = FALSE], na.rm = TRUE))
    cat("  Batch", b, "mean after correction:", round(new_mean, 2), "\n")
  }

  cat("DWD correction complete\n")

  corrected
}


#' Apply quantile normalization across batches
#'
#' Forces all samples to have the same distribution.
#' Aggressive but effective for strong batch effects.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector (used for reporting only)
#' @return Quantile-normalized expression matrix
#' @export
normalize_quantile <- function(exprs, batch = NULL) {

  cat("\n=== Quantile Normalization ===\n")

  if (!requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("Package 'preprocessCore' required. Install with: BiocManager::install('preprocessCore')")
  }

  exprs <- as.matrix(exprs)
  cat("Input:", nrow(exprs), "genes x", ncol(exprs), "samples\n")

  # Apply quantile normalization
  normalized <- preprocessCore::normalize.quantiles(exprs)
  rownames(normalized) <- rownames(exprs)
  colnames(normalized) <- colnames(exprs)

  cat("Quantile normalization complete\n")

  normalized
}


#' Apply batch mean-centering (simple baseline)
#'
#' Simplest batch correction: subtract batch mean from each gene.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @return Mean-centered expression matrix
#' @export
normalize_mean_center <- function(exprs, batch) {

  cat("\n=== Batch Mean Centering ===\n")

  batch <- as.factor(batch)
  batches <- levels(batch)

  cat("Batches:", length(batches), "\n")

  exprs <- as.matrix(exprs)
  corrected <- exprs

  # Calculate global mean
  global_mean <- rowMeans(exprs, na.rm = TRUE)

  # Center each batch to global mean

  for (b in batches) {
    batch_idx <- which(batch == b)
    batch_mean <- rowMeans(exprs[, batch_idx, drop = FALSE], na.rm = TRUE)

    # Shift to global mean
    shift <- global_mean - batch_mean
    corrected[, batch_idx] <- sweep(exprs[, batch_idx, drop = FALSE], 1, shift, "+")
  }

  cat("Mean centering complete\n")

  corrected
}


#' Apply RUV batch correction (RUVg via SVD on control genes)
#'
#' Estimates k unwanted factors W from negative control genes using SVD,
#' equivalent to RUVg (Gagnon-Bartsch & Speed 2012). W is returned
#' separately and should be included as covariates in the limma design
#' matrix. The expression matrix itself is returned unchanged.
#'
#' @param exprs Expression matrix (genes x samples), log2 scale
#' @param control_genes Character vector of control gene IDs (matching rownames)
#' @param k Number of unwanted factors to estimate (default: 2)
#' @return List with $exprs (unchanged matrix) and $W (samples x k matrix)
#' @export
normalize_ruv <- function(exprs, control_genes, k = 2) {

  cat("\n=== RUV Batch Correction (RUVg) ===\n")

  exprs <- as.matrix(exprs)
  control_idx <- which(rownames(exprs) %in% control_genes)
  cat(sprintf("Control genes in matrix: %d / %d provided\n",
              length(control_idx), length(control_genes)))

  if (length(control_idx) < k + 1)
    stop("Not enough control genes in matrix for k=", k,
         " (have ", length(control_idx), ")")

  max_k <- min(length(control_idx) - 1, ncol(exprs) - 1)
  if (k > max_k) {
    cat(sprintf("Warning: k=%d exceeds max=%d, capping\n", k, max_k))
    k <- max_k
  }

  Y_c <- t(exprs[control_idx, , drop = FALSE])
  Y_c <- scale(Y_c, center = TRUE, scale = FALSE)

  svd_c <- svd(Y_c, nu = k, nv = 0)
  W <- svd_c$u[, seq_len(k), drop = FALSE]
  colnames(W) <- paste0("W", seq_len(k))
  rownames(W) <- colnames(exprs)

  var_explained <- round(100 * svd_c$d[seq_len(k)]^2 / sum(svd_c$d^2), 1)
  cat(sprintf("Estimated %d unwanted factors\n", k))
  cat(sprintf("Variance explained by W: %s%%\n",
              paste(var_explained, collapse = "%, ")))

  list(exprs = exprs, W = W)
}


#' BRUV: Balanced RUV — estimate W from balanced datasets only
#'
#' Standard RUVg estimates W from all samples, which fails when batch and
#' biology are confounded (unbalanced datasets contribute only one group).
#' BRUV restricts the SVD to samples from "balanced" datasets — those with
#' samples in both baseline and contrast groups — where batch and biology
#' are decoupled.  The resulting W is then projected onto all samples
#' (including those from unbalanced datasets) for use in the limma design.
#'
#' @param exprs Expression matrix (genes x samples), log2 scale
#' @param control_genes Character vector of control gene IDs
#' @param k Number of unwanted factors (default: 2)
#' @param sample_dataset Factor/character vector mapping each sample (column)
#'   to its source dataset
#' @param sample_group Factor/character vector mapping each sample to its
#'   biological group (baseline or contrast)
#' @return List with $exprs (unchanged) and $W (samples x k matrix)
normalize_bruv <- function(exprs, control_genes, k = 2,
                           sample_dataset, sample_group) {

  cat("\n=== BRUV (Balanced RUV) Batch Correction ===\n")

  exprs <- as.matrix(exprs)
  n_samples <- ncol(exprs)

  # Identify balanced datasets (have samples in both groups)
  ds <- as.character(sample_dataset)
  grp <- as.character(sample_group)
  ds_groups <- split(grp, ds)
  balanced_ds <- names(which(sapply(ds_groups, function(g) length(unique(g)) > 1)))
  unbalanced_ds <- setdiff(unique(ds), balanced_ds)

  cat(sprintf("  Total datasets: %d\n", length(unique(ds))))
  cat(sprintf("  Balanced datasets (both groups): %d — %s\n",
              length(balanced_ds), paste(balanced_ds, collapse = ", ")))
  cat(sprintf("  Unbalanced datasets (one group): %d — %s\n",
              length(unbalanced_ds), paste(unbalanced_ds, collapse = ", ")))

  balanced_idx <- which(ds %in% balanced_ds)
  cat(sprintf("  Balanced samples: %d / %d\n", length(balanced_idx), n_samples))

  if (length(balanced_idx) < k + 2)
    stop("Not enough balanced samples for BRUV with k=", k)

  # Filter control genes
  control_idx <- which(rownames(exprs) %in% control_genes)
  cat(sprintf("  Control genes in matrix: %d / %d provided\n",
              length(control_idx), length(control_genes)))

  if (length(control_idx) < k + 1)
    stop("Not enough control genes for k=", k)

  max_k <- min(length(control_idx) - 1, length(balanced_idx) - 1)
  if (k > max_k) {
    cat(sprintf("  Warning: k=%d exceeds max=%d for balanced subset, capping\n", k, max_k))
    k <- max_k
  }

  # SVD on control genes from BALANCED samples only
  Y_c_balanced <- t(exprs[control_idx, balanced_idx, drop = FALSE])
  Y_c_balanced <- scale(Y_c_balanced, center = TRUE, scale = FALSE)

  svd_bal <- svd(Y_c_balanced, nu = k, nv = k)
  W_balanced <- svd_bal$u[, seq_len(k), drop = FALSE]

  var_explained <- round(100 * svd_bal$d[seq_len(k)]^2 / sum(svd_bal$d^2), 1)
  cat(sprintf("  Estimated %d unwanted factors from balanced subset\n", k))
  cat(sprintf("  Variance explained by W (balanced): %s%%\n",
              paste(var_explained, collapse = "%, ")))

  # Project W onto ALL samples using the right singular vectors (V)
  # For balanced samples: W_bal = U_k
  # For all samples: W_all = Y_c_all %*% V_k %*% diag(1/d_k)
  # This projects each sample (including unbalanced) into the same
  # factor space estimated from the balanced subset.
  Y_c_all <- t(exprs[control_idx, , drop = FALSE])
  bal_center <- colMeans(t(exprs[control_idx, balanced_idx, drop = FALSE]))
  Y_c_all <- sweep(Y_c_all, 2, bal_center, "-")

  V_k <- svd_bal$v[, seq_len(k), drop = FALSE]
  d_k <- svd_bal$d[seq_len(k)]
  W_all <- Y_c_all %*% V_k %*% diag(1 / d_k, nrow = k, ncol = k)

  colnames(W_all) <- paste0("W", seq_len(k))
  rownames(W_all) <- colnames(exprs)

  cat(sprintf("  Projected W onto all %d samples\n", n_samples))

  list(exprs = exprs, W = W_all)
}


normalize_ruvinv <- function(exprs, control_genes, group, covariates_df = NULL,
                              lambda = NULL) {

  cat("\n=== RUV Batch Correction (RUVinv) ===\n")
  cat("  No k parameter — integrates over all unwanted factors\n")

  exprs <- as.matrix(exprs)
  ctl <- rownames(exprs) %in% control_genes
  cat(sprintf("  Control genes in matrix: %d / %d provided\n",
              sum(ctl), length(control_genes)))

  if (sum(ctl) < 3)
    stop("Not enough control genes for RUVinv (have ", sum(ctl), ")")

  Y <- t(exprs)  # samples x genes (ruv convention)

  X <- as.numeric(group) - 1
  X <- matrix(X, ncol = 1)
  colnames(X) <- "group"

  Z <- NULL
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    Z <- model.matrix(~ ., data = covariates_df)[, -1, drop = FALSE]
    cat(sprintf("  Covariates in Z: %s\n", paste(colnames(Z), collapse = ", ")))
  }

  if (!is.null(lambda)) {
    cat(sprintf("  Using ridge parameter lambda = %g (RUVrinv)\n", lambda))
  }

  result <- ruv::RUVinv(Y = Y, X = X, ctl = ctl, Z = Z, lambda = lambda)

  betahat <- as.numeric(result$betahat)
  pvals   <- as.numeric(result$p)
  tvals   <- as.numeric(result$t)
  sigma2  <- as.numeric(result$sigma2)

  adj_pvals <- p.adjust(pvals, method = "BH")

  ave_expr <- colMeans(Y, na.rm = TRUE)

  de_results <- data.frame(
    gene      = colnames(Y),
    logFC     = betahat,
    AveExpr   = ave_expr,
    t         = tvals,
    P.Value   = pvals,
    adj.P.Val = adj_pvals,
    B         = NA_real_,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  RUVinv complete: %d genes tested\n", nrow(de_results)))
  cat(sprintf("  Effective df: %s\n", paste(round(result$df, 1), collapse = ", ")))

  list(exprs = exprs, de_results = de_results, W = result$W)
}


#' HarmonizR-style batch correction via matrix dissection
#'
#' Groups genes by batch membership pattern (which batches have real data),
#' runs ComBat once per group on all its batches, then reassembles.
#' Each gene is corrected exactly once. No imputation, no averaging.
#'
#' Based on: Voß et al. (2022) Nature Communications 13:6512
#'
#' @param incomplete_matrix Expression matrix with NAs (genes x samples)
#' @param batch Batch factor (length = ncol)
#' @param mod_data Data frame for ComBat model (bio_group + covariates)
#' @param mod_formula RHS of model formula as string (e.g. "bio_group + sex")
#' @param ref_batch Optional reference batch for ref-batch ComBat
#' @return Corrected matrix, same dims as input, NAs where gene is in only 1 batch
normalize_harmonizr <- function(incomplete_matrix, batch, mod_data = NULL,
                                mod_formula = NULL, ref_batch = NULL) {

  cat("\n=== HarmonizR (Matrix Dissection by Batch Membership) ===\n")

  incomplete_matrix <- as.matrix(incomplete_matrix)
  batch <- as.factor(batch)
  batches <- levels(batch)
  n_batches <- length(batches)
  n_genes <- nrow(incomplete_matrix)
  n_samples <- ncol(incomplete_matrix)

  cat(sprintf("  Input: %d genes x %d samples, %d batches\n",
              n_genes, n_samples, n_batches))
  cat(sprintf("  NA cells: %d / %d (%.1f%%)\n",
              sum(is.na(incomplete_matrix)), length(incomplete_matrix),
              100 * mean(is.na(incomplete_matrix))))

  result <- matrix(NA_real_, n_genes, n_samples,
                   dimnames = dimnames(incomplete_matrix))

  # Step 1: compute batch membership per gene (>=2 non-NA values per batch)
  batch_membership <- matrix(FALSE, n_genes, n_batches,
                             dimnames = list(rownames(incomplete_matrix), batches))
  for (b in batches) {
    b_cols <- which(batch == b)
    batch_membership[, b] <- rowSums(!is.na(incomplete_matrix[, b_cols, drop = FALSE])) >= 2
  }
  n_batches_per_gene <- rowSums(batch_membership)

  # Step 2: group genes by batch membership pattern
  membership_keys <- apply(batch_membership, 1, function(row) {
    paste(batches[row], collapse = "|")
  })

  single_batch <- n_batches_per_gene <= 1
  multi_batch  <- n_batches_per_gene >= 2

  cat(sprintf("  Single-batch genes (no correction): %d\n", sum(single_batch)))
  cat(sprintf("  Multi-batch genes (to correct): %d\n", sum(multi_batch)))

  # Copy single-batch genes unchanged
  for (g in which(single_batch)) {
    result[g, ] <- incomplete_matrix[g, ]
  }

  # Step 3: process each batch-membership group
  groups <- split(names(membership_keys[multi_batch]),
                  membership_keys[multi_batch])

  cat(sprintf("  Batch-membership groups: %d\n", length(groups)))

  n_corrected_genes <- 0
  n_fallback_genes  <- 0

  for (key in names(groups)) {
    gene_names <- groups[[key]]
    group_batches <- strsplit(key, "\\|")[[1]]

    # Samples from the relevant batches
    sub_cols <- which(batch %in% group_batches)
    sub_batch <- droplevels(batch[sub_cols])

    # Extract sub-matrix
    sub_mat <- incomplete_matrix[gene_names, sub_cols, drop = FALSE]

    # Drop genes with any remaining NA (ComBat requires complete data)
    complete <- complete.cases(sub_mat)
    if (sum(complete) < 2) {
      # Fall back: copy uncorrected
      for (g in gene_names) result[g, ] <- incomplete_matrix[g, ]
      n_fallback_genes <- n_fallback_genes + length(gene_names)
      next
    }
    incomplete_in_group <- gene_names[!complete]
    gene_names <- gene_names[complete]
    sub_mat <- sub_mat[complete, , drop = FALSE]

    # Copy incomplete genes uncorrected
    for (g in incomplete_in_group) result[g, ] <- incomplete_matrix[g, ]
    n_fallback_genes <- n_fallback_genes + length(incomplete_in_group)

    # Drop zero-variance genes
    gv <- apply(sub_mat, 1, var)
    zv_genes <- gene_names[gv <= 0]
    sub_mat <- sub_mat[gv > 0, , drop = FALSE]
    gene_names <- gene_names[gv > 0]
    for (g in zv_genes) result[g, ] <- incomplete_matrix[g, ]
    n_fallback_genes <- n_fallback_genes + length(zv_genes)

    if (nrow(sub_mat) < 2) {
      for (g in gene_names) result[g, ] <- incomplete_matrix[g, ]
      n_fallback_genes <- n_fallback_genes + length(gene_names)
      next
    }

    # Build mod matrix with fallback chain
    local_mod <- NULL
    if (!is.null(mod_data) && !is.null(mod_formula)) {
      local_md <- mod_data[sub_cols, , drop = FALSE]
      # Try 1: full formula (bio_group + covariates)
      tryCatch({
        m <- model.matrix(as.formula(paste("~", mod_formula)), data = local_md)
        if (qr(m)$rank >= ncol(m)) local_mod <- m
      }, error = function(e) NULL)
      # Try 2: bio_group only (drop covariates that may confound)
      if (is.null(local_mod) && "bio_group" %in% colnames(local_md)) {
        tryCatch({
          m <- model.matrix(~ bio_group, data = local_md)
          if (qr(m)$rank >= ncol(m)) local_mod <- m
        }, error = function(e) NULL)
      }
    }

    use_ref <- if (!is.null(ref_batch) && ref_batch %in% levels(sub_batch)) ref_batch else NULL

    # Run ComBat with fallback: with mod → without mod → uncorrected
    corrected_sub <- NULL

    # Attempt 1: ComBat with mod
    corrected_sub <- tryCatch({
      ComBat(dat = sub_mat, batch = sub_batch, mod = local_mod,
             par.prior = TRUE, prior.plots = FALSE, ref.batch = use_ref)
    }, error = function(e) NULL)

    # Attempt 2: ComBat without mod
    if (is.null(corrected_sub)) {
      corrected_sub <- tryCatch({
        ComBat(dat = sub_mat, batch = sub_batch, mod = NULL,
               par.prior = TRUE, prior.plots = FALSE, ref.batch = use_ref)
      }, error = function(e) NULL)
    }

    if (is.null(corrected_sub)) {
      cat(sprintf("    Group [%s]: %d genes — ComBat failed, keeping uncorrected\n",
                  key, length(gene_names)))
      for (g in gene_names) result[g, ] <- incomplete_matrix[g, ]
      n_fallback_genes <- n_fallback_genes + length(gene_names)
      next
    }

    # Store corrected values (only in columns belonging to this group's batches)
    for (g in rownames(corrected_sub)) {
      result[g, sub_cols] <- corrected_sub[g, ]
    }
    n_corrected_genes <- n_corrected_genes + nrow(corrected_sub)

    if (nrow(corrected_sub) >= 100) {
      cat(sprintf("    Group [%d batches]: %d genes corrected\n",
                  length(group_batches), nrow(corrected_sub)))
    }
  }

  # Summary
  n_na_cells <- sum(is.na(result))
  n_real_cells <- sum(!is.na(result))

  cat(sprintf("\n  Corrected genes: %d, uncorrected fallback: %d, single-batch: %d\n",
              n_corrected_genes, n_fallback_genes, sum(single_batch)))
  cat(sprintf("  Result cells: %d real, %d NA\n", n_real_cells, n_na_cells))
  cat("HarmonizR correction complete\n")

  result
}


#' Batch correction via the HarmonizR Bioconductor package
#'
#' Wraps HarmonizR::harmonizR() which performs matrix dissection + ComBat/limma.
#' See: Voß et al. (2022) Nature Communications 13:6512
#' GitHub: https://github.com/HSU-HPC/HarmonizR
#'
#' @param incomplete_matrix Expression matrix with NAs (genes x samples)
#' @param batch Batch factor (length = ncol)
#' @param algorithm "ComBat" or "limma"
#' @param ComBat_mode 1-4 (1=parametric prior, 2=non-parametric, 3=no mean adj, 4=no adj)
#' @param sort Sorting method: FALSE, "sparsity_sort", "seriation_sort", "jaccard_sort"
#' @param block How many batches to treat as one (NULL = off)
#' @return Corrected matrix, same dims as input, NAs remain where gene is single-batch
normalize_github_harmonizr <- function(incomplete_matrix, batch,
                                       algorithm = "ComBat", ComBat_mode = 1,
                                       sort = FALSE, block = NULL) {

  cat("\n=== GitHub HarmonizR (Bioconductor package) ===\n")

  if (!requireNamespace("HarmonizR", quietly = TRUE))
    stop("HarmonizR package not installed. Install with: ",
         'devtools::install_github("HSU-HPC/HarmonizR")')

  incomplete_matrix <- as.matrix(incomplete_matrix)
  batch <- as.factor(batch)

  cat(sprintf("  Input: %d genes x %d samples, %d batches\n",
              nrow(incomplete_matrix), ncol(incomplete_matrix), nlevels(batch)))
  cat(sprintf("  NA cells: %d / %d (%.1f%%)\n",
              sum(is.na(incomplete_matrix)), length(incomplete_matrix),
              100 * mean(is.na(incomplete_matrix))))
  cat(sprintf("  Algorithm: %s, ComBat_mode: %d, sort: %s\n",
              algorithm, ComBat_mode, as.character(sort)))

  data_df <- as.data.frame(incomplete_matrix)

  desc_df <- data.frame(
    ID    = colnames(incomplete_matrix),
    sample = seq_len(ncol(incomplete_matrix)),
    batch  = as.integer(batch),
    stringsAsFactors = FALSE
  )

  result_df <- HarmonizR::harmonizR(
    data_as_input        = data_df,
    description_as_input = desc_df,
    algorithm            = algorithm,
    ComBat_mode          = ComBat_mode,
    sort                 = sort,
    block                = block,
    output_file          = FALSE,
    verbosity            = 1,
    cores                = 1
  )

  result <- as.matrix(result_df)

  n_genes_out <- nrow(result)
  n_genes_lost <- nrow(incomplete_matrix) - n_genes_out
  cat(sprintf("  Output: %d genes (%d dropped by HarmonizR)\n",
              n_genes_out, n_genes_lost))

  if (n_genes_lost > 0) {
    full_result <- matrix(NA_real_, nrow(incomplete_matrix), ncol(incomplete_matrix),
                          dimnames = dimnames(incomplete_matrix))
    shared <- intersect(rownames(result), rownames(full_result))
    full_result[shared, ] <- result[shared, ]
    result <- full_result
  }

  cat("GitHub HarmonizR correction complete\n")
  result
}


#' Compare normalization methods
#'
#' Runs multiple normalization methods and compares their effects
#' using PCA and batch effect metrics.
#'
#' @param exprs Expression matrix (genes x samples)
#' @param batch Batch vector
#' @param biological_group Biological group of interest (for preservation)
#' @param methods Vector of methods to test
#' @return List with normalized matrices and comparison metrics
#' @export
compare_normalizations <- function(exprs,
                                    batch,
                                    biological_group = NULL,
                                    methods = c("none", "combat", "dwd", "mean_center")) {

  cat("\n=== Comparing Normalization Methods ===\n\n")

  results <- list()
  metrics <- list()

  for (method in methods) {
    cat("Method:", method, "\n")

    normalized <- switch(
      method,
      "none" = as.matrix(exprs),
      "combat" = normalize_combat(exprs, batch),
      "dwd" = normalize_dwd(exprs, batch),
      "mean_center" = normalize_mean_center(exprs, batch),
      "quantile" = normalize_quantile(exprs, batch),
      stop("Unknown method: ", method)
    )

    results[[method]] <- normalized

    # Calculate batch effect metrics
    metrics[[method]] <- calculate_batch_metrics(normalized, batch, biological_group)
    cat("\n")
  }

  # Summary comparison
  cat("\n=== Normalization Comparison Summary ===\n\n")
  metrics_df <- do.call(rbind, lapply(names(metrics), function(m) {
    data.frame(
      method = m,
      batch_variance_pct = metrics[[m]]$batch_variance_pct,
      bio_variance_pct = metrics[[m]]$bio_variance_pct,
      silhouette_batch = metrics[[m]]$silhouette_batch,
      stringsAsFactors = FALSE
    )
  }))
  print(metrics_df)

  list(
    normalized = results,
    metrics = metrics,
    summary = metrics_df
  )
}


#' Calculate batch effect metrics
#'
#' @param exprs Expression matrix
#' @param batch Batch vector
#' @param biological_group Biological group (optional)
#' @return List of metrics
calculate_batch_metrics <- function(exprs, batch, biological_group = NULL) {

  # PCA
  pca <- prcomp(t(exprs), scale. = TRUE, center = TRUE)

  # Variance explained by batch (using PC1-PC5)
  pc_scores <- pca$x[, 1:min(5, ncol(pca$x))]

  batch_variance <- 0
  bio_variance <- 0

  for (i in 1:ncol(pc_scores)) {
    # ANOVA for batch effect
    batch_anova <- summary(aov(pc_scores[, i] ~ batch))[[1]]
    batch_variance <- batch_variance + batch_anova$`Sum Sq`[1] / sum(batch_anova$`Sum Sq`)

    # ANOVA for biological effect
    if (!is.null(biological_group)) {
      bio_anova <- summary(aov(pc_scores[, i] ~ biological_group))[[1]]
      bio_variance <- bio_variance + bio_anova$`Sum Sq`[1] / sum(bio_anova$`Sum Sq`)
    }
  }

  batch_variance_pct <- round(100 * batch_variance / ncol(pc_scores), 1)
  bio_variance_pct <- round(100 * bio_variance / ncol(pc_scores), 1)

  # Silhouette score for batch clustering
  silhouette_batch <- NA
  if (requireNamespace("cluster", quietly = TRUE) && nlevels(as.factor(batch)) > 1) {
    dist_mat <- dist(pc_scores)
    sil <- cluster::silhouette(as.numeric(as.factor(batch)), dist_mat)
    silhouette_batch <- round(mean(sil[, 3]), 3)
  }

  cat("  Batch variance (PC1-5):", batch_variance_pct, "%\n")
  if (!is.null(biological_group)) {
    cat("  Biological variance (PC1-5):", bio_variance_pct, "%\n")
  }
  cat("  Batch silhouette:", silhouette_batch, "\n")

  list(
    batch_variance_pct = batch_variance_pct,
    bio_variance_pct = bio_variance_pct,
    silhouette_batch = silhouette_batch,
    pca = pca
  )
}

# Imputation methods.
#
# Every method takes the list returned by create_incomplete_matrix()
# (lib/gene_coverage.R) and returns a list whose `matrix` is the completed
# genes x samples matrix; observed cells are never changed.
# The article's method is softImpute; the others are the comparators of the
# holdout validation. run_imputer(name, incomplete, cfg) dispatches by name.

library(softImpute)

#' Number of non-zero singular values ("qualities") of a softImpute fit
#' @export
softimpute_rank <- function(fit) {
  if (is.null(fit$d)) 0L else sum(fit$d > 1e-6)
}


#' softImpute warm-start lambda path
#'
#' Fits the geometric grid lambda_k = lambda0 * min_ratio^(k / (n_steps - 1)),
#' k = 0, ..., stop_step, each fit initialised from the previous one
#' (Mazumder, Hastie & Tibshirani 2010). lambda0 is computed on the centred
#' matrix itself, so the grid is relative and transfers unchanged to masked
#' copies of the matrix. Shared by impute_softimpute() and
#' scripts/analysis/softimpute_lambda_selection.R, which chooses stop_step by hold-out.
#'
#' @param X_centered biScale()d incomplete matrix
#' @param stop_step Last grid index fitted (0-based); default the full path
#' @param on_step Optional callback(step, lambda, fit), e.g. to score held-out
#'   cells at every step
#' @return List with: fits (one per step), lambda (grid up to stop_step),
#'   rank (qualities per step), lambda0
#' @export
softimpute_path <- function(X_centered, rank_max, thresh, maxit, type = "als",
                            n_steps = 8L, min_ratio = 0.01,
                            stop_step = n_steps - 1L, on_step = NULL) {
  stop_step <- as.integer(stop_step)
  if (stop_step < 0L || stop_step > n_steps - 1L)
    stop("path_stop_step must lie in 0..", n_steps - 1L)
  lambda_max <- lambda0(X_centered)
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_max * min_ratio), length.out = n_steps))
  lambda_seq <- lambda_seq[seq_len(stop_step + 1L)]
  fits <- vector("list", length(lambda_seq))
  ranks <- integer(length(lambda_seq))
  fit <- NULL
  for (i in seq_along(lambda_seq)) {
    fit <- softImpute(
      X_centered,
      rank.max = rank_max,
      lambda = lambda_seq[i],
      type = type,
      thresh = thresh,
      maxit = maxit,
      warm.start = fit
    )
    fits[[i]] <- fit
    ranks[i] <- softimpute_rank(fit)
    cat(sprintf("  step %d: lambda = %.4f, rank = %d%s\n", i - 1L, lambda_seq[i],
                ranks[i], if (ranks[i] >= rank_max) "  ** RANK CAP REACHED **" else ""))
    if (!is.null(on_step)) on_step(i - 1L, lambda_seq[i], fit)
  }
  list(fits = fits, lambda = lambda_seq, rank = ranks, lambda0 = lambda_max)
}


#' Impute missing values using softImpute
#'
#' Uses nuclear norm regularization to complete the matrix,
#' exploiting low-rank structure typical of gene expression data.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param rank_max Maximum rank for approximation. Should stay above the rank
#'   of the selected fit; reaching it is reported as a warning.
#' @param lambda Regularization parameter (0 = warm-start path, see
#'   path_stop_step)
#' @param thresh Convergence threshold
#' @param maxit Maximum iterations
#' @param type Algorithm type: "als" or "svd"
#' @param seed RNG seed set before fitting. The ALS solver starts from a
#'   random matrix, so without a fixed seed each run converges to a
#'   slightly different completion.
#' @param path_stop_step With lambda = 0: index (0-based) of the 8-point
#'   warm-start grid at which the path stops; the fit there is used. Chosen
#'   by hold-out in scripts/analysis/softimpute_lambda_selection.R. Default 7 = the full
#'   path to 0.01 * lambda0 (the original setting).
#' @return List with: matrix (imputed), fit (softImpute object), lambda, rank
#' @export
impute_softimpute <- function(incomplete,
                               rank_max = 50,
                               lambda = 0,
                               thresh = 1e-5,
                               maxit = 100,
                               type = "als",
                               seed = 42,
                               path_stop_step = 7L) {

  cat("\n=== softImpute Matrix Completion ===\n")
  set.seed(seed)
  cat("RNG seed:", seed, "\n")

  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "/", length(X),
      "(", round(100 * sum(is.na(X)) / length(X), 1), "%)\n")

  # Center and scale for better convergence (reduced iterations for speed)
  X_centered <- biScale(X, maxit = 20, thresh = 1e-3)

  # lambda = 0: warm-start path over the fixed relative grid, stopped at
  # path_stop_step (not auto-selected here; the step is chosen by hold-out
  # in scripts/analysis/softimpute_lambda_selection.R and set in the config)
  if (lambda == 0) {
    cat("Warm-start lambda path, stopping at step", path_stop_step, "...\n")
    path <- softimpute_path(X_centered, rank_max = rank_max, thresh = thresh,
                            maxit = maxit, type = type,
                            stop_step = path_stop_step)
    final_fit <- path$fits[[length(path$fits)]]
    lambda_used <- path$lambda[length(path$lambda)]
    cat("lambda0 =", round(path$lambda0, 4), "; using lambda =",
        round(lambda_used, 4), "\n")

  } else {
    cat("Using specified lambda =", lambda, "\n")
    final_fit <- softImpute(
      X_centered,
      rank.max = rank_max,
      lambda = lambda,
      type = type,
      thresh = thresh,
      maxit = maxit
    )
    lambda_used <- lambda
  }

  # Complete the matrix. softImpute::complete() reads the biScale
  # attributes on X_centered and returns values on the ORIGINAL scale,
  # so no manual un-scaling is needed (and doing so double-corrects
  # and produces values off by row+column centers).
  X_imputed <- complete(X_centered, final_fit)

  # Preserve original values (only impute NAs)
  X_final <- incomplete$matrix
  X_final[is.na(X_final)] <- X_imputed[is.na(X_final)]

  # Report
  final_rank <- softimpute_rank(final_fit)
  cat("\nImputation complete:\n")
  cat("  Final rank:", final_rank, "( rank_max", rank_max, ")\n")
  if (final_rank >= rank_max)
    warning("softImpute final rank reached rank_max = ", rank_max,
            ": the cap, not lambda, set the number of qualities")
  cat("  Imputed values:", sum(is.na(incomplete$matrix)), "\n")

  list(
    matrix = X_final,
    fit = final_fit,
    lambda = lambda_used,
    rank = final_rank,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Impute missing values using KNN
#'
#' Simple k-nearest neighbors imputation based on gene similarity.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param k Number of neighbors
#' @return List with: matrix (imputed), validation
#' @export
impute_knn <- function(incomplete, k = 10) {

  cat("\n=== KNN Imputation ===\n")

  if (!requireNamespace("impute", quietly = TRUE)) {
    stop("Package 'impute' required. Install with: BiocManager::install('impute')")
  }

  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")

  # impute.knn expects genes as rows
  X_imputed <- impute::impute.knn(X, k = k)$data

  cat("Imputation complete\n")

  list(
    matrix = X_imputed,
    k = k,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Impute missing values using iterative regularized PCA (missMDA)
#'
#' Wraps \code{missMDA::imputePCA}, an EM-style algorithm that
#' alternates between a low-rank PCA fit and re-imputation of the
#' missing cells. Conceptually closest to \code{softImpute} --- both
#' exploit low-rank structure --- but uses a different regularisation
#' scheme and convergence criterion.
#'
#' missMDA expects observations in rows and variables in columns, so we
#' transpose internally (samples become rows, genes become columns).
#'
#' @param incomplete List from create_incomplete_matrix
#' @param ncp Number of principal components used for imputation
#' @param scale Whether to scale variables to unit variance before PCA
#' @param method Either "Regularized" or "EM"
#' @param input "raw" (genes centred, and scaled if \code{scale}, by
#'   imputePCA itself) or "biscale" (gene and sample offsets and scales
#'   removed by \code{biScale} first, as softImpute gets its input; the
#'   result is mapped back to the original scale)
#' @return List with: matrix (imputed), ncp, gene_info, sample_info
#' @export
impute_missmda <- function(incomplete, ncp = 5, scale = TRUE,
                           method = "Regularized", input = "raw") {

  cat("\n=== missMDA imputePCA ===\n")

  if (!requireNamespace("missMDA", quietly = TRUE)) {
    stop("Package 'missMDA' required. Install with: install.packages('missMDA')")
  }

  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")
  cat("ncp =", ncp, " method =", method, " input =", input, "\n")
  if (sum(is.na(X)) == 0L) {
    cat("Nothing to impute (no NAs). Returning input.\n")
    return(list(matrix = X, ncp = ncp,
                gene_info = incomplete$gene_info,
                sample_info = incomplete$sample_info))
  }
  cat("[", format(Sys.time(), "%H:%M:%S"),
      "] imputePCA EM loop starting (silent until convergence)...\n",
      sep = "")

  # imputePCA: observations = rows, variables = columns
  input <- match.arg(input, c("raw", "biscale"))
  if (input == "biscale") {
    Xc <- biScale(X, maxit = 20, thresh = 1e-3)
    X_samples_by_genes <- t(unclass(Xc))
    scale <- FALSE
  } else {
    X_samples_by_genes <- t(X)
  }
  t0 <- Sys.time()
  res <- missMDA::imputePCA(
    X_samples_by_genes,
    ncp = ncp,
    scale = scale,
    method = method
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  X_imputed <- t(res$completeObs)
  if (input == "biscale") {
    rw <- attr(Xc, "biScale:row"); cl <- attr(Xc, "biScale:column")
    X_imputed <- outer(rw$center, cl$center, "+") + outer(rw$scale, cl$scale) * X_imputed
  }
  # Preserve original cells exactly
  X_final <- incomplete$matrix
  X_final[is.na(X_final)] <- X_imputed[is.na(X_final)]

  cat(sprintf("[%s] Imputation complete (%.1fs)\n",
              format(Sys.time(), "%H:%M:%S"), elapsed))

  list(
    matrix = X_final,
    ncp = ncp,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}


#' Impute missing values using sample-space KNN
#'
#' A hand-rolled sample-space KNN imputer. Contrast with
#' \code{impute::impute.knn}, which searches in gene space. Here each
#' missing cell is filled by averaging k nearest SAMPLES (columns)
#' that have an observed value for that gene. Sample-sample distances
#' are computed once on the fully-observed gene subset, standardised
#' per gene so each contributes equally.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param k Number of neighbours
#' @return List with: matrix (imputed), k, gene_info, sample_info
#' @export
impute_sample_knn <- function(incomplete, k = 10) {

  cat("\n=== Sample-space KNN Imputation ===\n")

  X <- incomplete$matrix
  n_genes <- nrow(X)
  n_samples <- ncol(X)
  cat("Input matrix:", n_genes, "genes x", n_samples, "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")
  cat("k =", k, "\n")

  t0 <- Sys.time()

  na_per_gene <- rowSums(is.na(X))
  basis_idx <- which(na_per_gene == 0)
  cat(sprintf("[%s] Distance basis: %d fully-observed genes\n",
              format(Sys.time(), "%H:%M:%S"), length(basis_idx)))
  if (length(basis_idx) < 10) {
    # Scattered missingness (e.g. random-cell masks) can leave no gene
    # fully observed: fall back to every gene, missing cells filled with
    # the gene's observed mean, for the distances only.
    cat("  Too few fully-observed genes; distance basis = all genes, NAs at gene mean\n")
    basis <- X
    na_b <- which(is.na(basis), arr.ind = TRUE)
    basis[na_b] <- rowMeans(X, na.rm = TRUE)[na_b[, 1]]
  } else {
    basis <- X[basis_idx, , drop = FALSE]
  }
  row_mu <- rowMeans(basis)
  row_sd <- apply(basis, 1, sd)
  row_sd[row_sd == 0 | is.na(row_sd)] <- 1
  basis_std <- (basis - row_mu) / row_sd

  d_mat <- as.matrix(stats::dist(t(basis_std)))
  diag(d_mat) <- Inf
  cat(sprintf("[%s] Sample distance matrix: %dx%d\n",
              format(Sys.time(), "%H:%M:%S"),
              nrow(d_mat), ncol(d_mat)))

  na_genes <- which(na_per_gene > 0)
  n_to_impute <- length(na_genes)
  cat(sprintf("[%s] Genes needing imputation: %d / %d (%.1f%%)\n",
              format(Sys.time(), "%H:%M:%S"),
              n_to_impute, n_genes,
              100 * n_to_impute / n_genes))

  if (n_to_impute == 0) {
    cat("Nothing to impute (no NAs). Returning input.\n")
    return(list(
      matrix = X,
      k = k,
      gene_info = incomplete$gene_info,
      sample_info = incomplete$sample_info
    ))
  }

  X_imputed <- X
  progress_every <- max(1L, floor(n_to_impute / 20))
  for (idx in seq_along(na_genes)) {
    g <- na_genes[idx]
    row_vals <- X[g, ]
    missing_samples <- which(is.na(row_vals))
    observed_samples <- which(!is.na(row_vals))
    if (length(observed_samples) == 0L) next
    for (s in missing_samples) {
      dists <- d_mat[s, observed_samples]
      ord <- order(dists)
      nn <- observed_samples[ord[seq_len(min(k, length(ord)))]]
      X_imputed[g, s] <- mean(row_vals[nn])
    }
    if (idx %% progress_every == 0L || idx == n_to_impute) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      pct <- 100 * idx / n_to_impute
      rate <- idx / max(elapsed, 1e-6)
      eta <- (n_to_impute - idx) / max(rate, 1e-6)
      cat(sprintf("[%s] %d/%d genes (%.0f%%)  elapsed=%.1fs  ETA=%.1fs\n",
                  format(Sys.time(), "%H:%M:%S"),
                  idx, n_to_impute, pct, elapsed, eta))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] Imputation complete (%.1fs)\n",
              format(Sys.time(), "%H:%M:%S"), elapsed))

  list(
    matrix = X_imputed,
    k = k,
    gene_info = incomplete$gene_info,
    sample_info = incomplete$sample_info
  )
}

#' Impute missing values using per-gene mean
#'
#' For each gene, replaces NAs with the mean of observed values for that gene.
#'
#' @param incomplete List from create_incomplete_matrix
#' @return List with: matrix (imputed)
#' @export
impute_gene_mean <- function(incomplete) {
  cat("\n=== Per-gene Mean Imputation ===\n")
  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")
  row_means <- rowMeans(X, na.rm = TRUE)
  na_idx <- which(is.na(X), arr.ind = TRUE)
  X[na_idx] <- row_means[na_idx[, 1]]
  cat("Imputed values:", nrow(na_idx), "\n")
  list(matrix = X)
}


#' Impute missing values using per-gene-per-batch mean
#'
#' For each gene and batch (dataset), replaces NAs with the mean of observed
#' values for that gene within samples from the same batch. Falls back to the
#' global gene mean for batches where the gene is entirely missing.
#'
#' @param incomplete List from create_incomplete_matrix
#' @return List with: matrix (imputed)
#' @export
impute_batch_mean <- function(incomplete) {
  cat("\n=== Per-gene-per-batch Mean Imputation ===\n")
  X <- incomplete$matrix
  sample_info <- incomplete$sample_info
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")

  datasets <- unique(sample_info$dataset)
  ds_cols <- split(seq_len(ncol(X)), sample_info$dataset)
  row_means <- rowMeans(X, na.rm = TRUE)

  n_imputed <- 0L
  for (ds in datasets) {
    cols <- ds_cols[[ds]]
    sub <- X[, cols, drop = FALSE]
    na_mask <- is.na(sub)
    if (!any(na_mask)) next
    batch_means <- rowMeans(sub, na.rm = TRUE)
    fallback <- is.nan(batch_means)
    batch_means[fallback] <- row_means[fallback]
    na_idx <- which(na_mask, arr.ind = TRUE)
    sub[na_idx] <- batch_means[na_idx[, 1]]
    X[, cols] <- sub
    n_imputed <- n_imputed + nrow(na_idx)
  }
  cat("Imputed values:", n_imputed, "\n")
  list(matrix = X)
}


#' Impute missing values with the biScale centres only
#'
#' Fits softImpute's biScale model x_ij = a_i + b_j + s_i * t_j * z_ij on the
#' observed cells (production settings: maxit 20, thresh 1e-3) and fills every
#' gap with z = 0, i.e. x_ij = a_i + b_j (gene centre + sample offset). This is
#' softImpute with zero qualities: the same values complete() returns for a fit
#' at lambda >= lambda0.
#'
#' @param incomplete List from create_incomplete_matrix
#' @return List with: matrix (imputed)
#' @export
impute_biscale_zero <- function(incomplete) {
  cat("\n=== biScale-only (zero-quality) Imputation ===\n")
  X <- incomplete$matrix
  cat("Input matrix:", nrow(X), "genes x", ncol(X), "samples\n")
  cat("Missing values:", sum(is.na(X)), "\n")
  Xc <- biScale(X, maxit = 20, thresh = 1e-3)
  row_center <- attr(Xc, "biScale:row")$center
  col_center <- attr(Xc, "biScale:column")$center
  na_idx <- which(is.na(X), arr.ind = TRUE)
  X[na_idx] <- row_center[na_idx[, 1]] + col_center[na_idx[, 2]]
  cat("Imputed values:", nrow(na_idx), "\n")
  list(matrix = X)
}


#' Run an imputer on the biScale-centred matrix and map the result back
#'
#' biScale removes gene and sample offsets and scales (as softImpute gets
#' its input); \code{imputer} fills the centred matrix; the fill is then
#' un-centred. Observed cells are returned unchanged.
#'
#' @param incomplete List from create_incomplete_matrix
#' @param imputer Function taking an incomplete-style list and returning a
#'   list with \code{matrix}
#' @return List with: matrix (imputed)
#' @export
impute_centred <- function(incomplete, imputer) {
  X <- incomplete$matrix
  Xc <- biScale(X, maxit = 20, thresh = 1e-3)
  inner <- incomplete
  inner$matrix <- matrix(unclass(Xc), nrow(X), dimnames = dimnames(X))
  Z <- imputer(inner)$matrix
  rw <- attr(Xc, "biScale:row"); cl <- attr(Xc, "biScale:column")
  X_back <- outer(rw$center, cl$center, "+") + outer(rw$scale, cl$scale) * Z
  X[is.na(X)] <- X_back[is.na(X)]
  list(matrix = X, gene_info = incomplete$gene_info,
       sample_info = incomplete$sample_info)
}


#' Imputation method registry
#'
#' Maps an imputation method name (as used in configs and CLI args) to a
#' uniform adapter function. Each adapter takes \code{(incomplete, cfg)}
#' where \code{incomplete} is the list returned by
#' \code{create_incomplete_matrix()} and \code{cfg} is a named list of
#' method-specific options (usually \code{config$imputation[[method]]}).
#' Every adapter must return a list with at least a \code{matrix} element
#' containing the imputed gene x sample matrix. Adding a new method
#' means adding one entry here.
#'
#' @export
IMPUTERS <- list(
  softimpute = function(incomplete, cfg = list()) {
    ## YAML 1.1 (R's yaml package) reads exponent notation without a dot, such
    ## as `thresh: 1e-4`, as a STRING. softImpute's ALS loop then compares
    ## `ratio > thresh` as strings ("1" > "1e-4" is FALSE) and never iterates,
    ## leaving the biScale gene+sample means as the "imputation". Coerce every
    ## numeric option so the config value is honoured regardless of notation.
    num <- function(v, d) if (is.null(v)) d else as.numeric(v)
    impute_softimpute(
      incomplete,
      rank_max = num(cfg$rank_max, 50),
      lambda   = num(cfg$lambda,   0),
      thresh   = num(cfg$thresh,   1e-5),
      maxit    = num(cfg$maxit,    100),
      type     = cfg$type     %||% "als",
      seed     = num(cfg$seed,     42),
      path_stop_step = num(cfg$path_stop_step, 7L)
    )
  },
  knn = function(incomplete, cfg = list()) {
    impute_knn(incomplete, k = cfg$k %||% 10)
  },
  missmda = function(incomplete, cfg = list()) {
    impute_missmda(
      incomplete,
      ncp    = cfg$ncp    %||% 5,
      scale  = cfg$scale  %||% TRUE,
      method = cfg$method %||% "Regularized",
      input  = cfg$input  %||% "raw"
    )
  },
  sample_knn = function(incomplete, cfg = list()) {
    impute_sample_knn(incomplete, k = cfg$k %||% 10)
  },
  knn_centred = function(incomplete, cfg = list()) {
    impute_centred(incomplete, function(x) impute_knn(x, k = cfg$k %||% 10))
  },
  sample_knn_centred = function(incomplete, cfg = list()) {
    impute_centred(incomplete, function(x) impute_sample_knn(x, k = cfg$k %||% 10))
  },
  biscale_zero = function(incomplete, cfg = list()) {
    impute_biscale_zero(incomplete)
  },
  gene_mean = function(incomplete, cfg = list()) {
    impute_gene_mean(incomplete)
  },
  batch_mean = function(incomplete, cfg = list()) {
    impute_batch_mean(incomplete)
  }
)


#' Run a registered imputation method by name
#'
#' Thin wrapper around \code{IMPUTERS} that gives a clear error for
#' unknown methods and keeps dispatch sites one-liners.
#'
#' @param method Name of a registered imputation method (e.g. "softimpute").
#' @param incomplete List from \code{create_incomplete_matrix()}.
#' @param cfg Named list of method-specific options. Typically
#'   \code{config$imputation[[method]]}.
#' @return Whatever the adapter returns (at least a list with
#'   \code{matrix}).
#' @export
run_imputer <- function(method, incomplete, cfg = list()) {
  if (!method %in% names(IMPUTERS)) {
    stop("Unknown imputation method: '", method,
         "'. Registered methods: ",
         paste(names(IMPUTERS), collapse = ", "))
  }
  IMPUTERS[[method]](incomplete, cfg %||% list())
}

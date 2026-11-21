# Masking schemes for holdout validation, and the leave-out accuracy check
# the main pipeline runs on its own imputer.
# Each mask function returns list(mask_idx, n_blocks).


#' Mask random individual cells
mask_random_cells <- function(observed_idx, n_mask_target) {
  mask_idx <- sample(observed_idx, min(n_mask_target, length(observed_idx)))
  list(mask_idx = mask_idx, n_blocks = NA_integer_)
}


#' Mask gene-dataset blocks, stopping at n_mask_target cells.
#' @param obs_pairs Data frame with gene/dataset columns
#' @param X The expression matrix (for computing linear indices)
#' @param ds_cols_cache Named list: dataset -> column indices
#' @param row_idx_map Named vector: gene name -> row index
#' @param n_mask_target Target number of cells to mask
#' @param min_obs_per_gene Minimum observed cells to keep per gene row
mask_gene_dataset_block <- function(obs_pairs, X, ds_cols_cache, row_idx_map,
                                    n_mask_target, min_obs_per_gene) {
  shuffled <- sample.int(nrow(obs_pairs))
  mask_idx_chunks <- vector("list", length(shuffled))
  n_masked_so_far <- 0L
  n_blocks <- 0L
  n_blocks_skipped_floor <- 0L
  gene_obs_remaining <- rowSums(!is.na(X))
  for (p in shuffled) {
    g  <- obs_pairs$gene[p]
    ds <- obs_pairs$dataset[p]
    cols <- ds_cols_cache[[ds]]
    r <- row_idx_map[[g]]
    block_idx <- (cols - 1L) * nrow(X) + r
    block_idx <- block_idx[!is.na(X[block_idx])]
    if (length(block_idx) == 0) next
    if (gene_obs_remaining[r] - length(block_idx) < min_obs_per_gene) {
      n_blocks_skipped_floor <- n_blocks_skipped_floor + 1L
      next
    }
    n_blocks <- n_blocks + 1L
    mask_idx_chunks[[n_blocks]] <- block_idx
    gene_obs_remaining[r] <- gene_obs_remaining[r] - length(block_idx)
    n_masked_so_far <- n_masked_so_far + length(block_idx)
    if (n_masked_so_far >= n_mask_target) break
  }
  mask_idx <- unlist(mask_idx_chunks[seq_len(n_blocks)], use.names = FALSE)
  cat("  Blocks selected:", n_blocks,
      " (skipped to keep >=", min_obs_per_gene, "obs/gene:",
      n_blocks_skipped_floor, ")\n")
  list(mask_idx = mask_idx, n_blocks = n_blocks)
}


#' Progressive-tax block masking: mask each gene down to n_min datasets,
#' but only a random subset of genes so total masked ≈ n_mask_target.
#' @param X The expression matrix
#' @param presence Boolean presence matrix (genes x datasets)
#' @param n_ds_per_gene Integer vector of per-gene dataset counts
#' @param ds_cols_cache Named list: dataset -> column indices
#' @param n_min Minimum coverage (genes at this level are untouched)
#' @param n_mask_target Target number of cells to mask
mask_progressive_tax_block <- function(X, presence, n_ds_per_gene,
                                       ds_cols_cache, n_min, n_mask_target) {
  genes_above <- which(n_ds_per_gene > n_min)
  # Estimate cells per gene if fully masked down to n_min, then sample
  # enough genes to hit the target.
  cells_per_gene <- numeric(length(genes_above))
  for (i in seq_along(genes_above)) {
    gi <- genes_above[i]
    gene_datasets <- names(which(presence[gi, ]))
    n_mask_ds <- length(gene_datasets) - n_min
    # estimate: avg cells per dataset-block for this gene
    total_obs <- sum(!is.na(X[gi, ]))
    cells_per_gene[i] <- total_obs * (n_mask_ds / length(gene_datasets))
  }
  # Shuffle and accumulate until we hit the target
  gene_order <- sample.int(length(genes_above))
  mask_idx_chunks <- list()
  n_blocks <- 0L
  n_genes_masked <- 0L
  n_masked_so_far <- 0L
  for (oi in gene_order) {
    gi <- genes_above[oi]
    gene_datasets <- names(which(presence[gi, ]))
    keep_ds <- sample(gene_datasets, n_min)
    mask_ds <- setdiff(gene_datasets, keep_ds)
    for (ds in mask_ds) {
      cols <- ds_cols_cache[[ds]]
      block_idx <- (cols - 1L) * nrow(X) + gi
      block_idx <- block_idx[!is.na(X[block_idx])]
      if (length(block_idx) == 0) next
      n_blocks <- n_blocks + 1L
      mask_idx_chunks[[n_blocks]] <- block_idx
      n_masked_so_far <- n_masked_so_far + length(block_idx)
    }
    n_genes_masked <- n_genes_masked + 1L
    if (n_masked_so_far >= n_mask_target) break
  }
  mask_idx <- unlist(mask_idx_chunks, use.names = FALSE)
  cat("  Genes masked down to", n_min, "datasets:", n_genes_masked,
      "/", length(genes_above),
      " | Blocks masked:", n_blocks, "\n")
  list(mask_idx = mask_idx, n_blocks = n_blocks)
}


#' Validate imputation accuracy using leave-out cross-validation
#'
#' Masks a fraction of observed values, imputes them, and measures
#' correlation with true values. Two masking strategies are supported:
#'
#' \itemize{
#'   \item \code{"random_cells"} — uniform random draw of individual
#'     matrix cells (the classical in-distribution leave-out).
#'   \item \code{"gene_dataset_block"} — mask all cells for a
#'     randomly chosen (gene, dataset) pair, stopping once the total
#'     number of masked cells is at least
#'     \code{leave_out_fraction * n_observed}. This simulates the
#'     real missingness pattern where a gene is absent from one
#'     entire dataset because the chip did not probe it.
#' }
#'
#' @param exprs_list Named list of expression matrices
#' @param min_datasets Minimum number of datasets a gene must appear in
#' @param leave_out_fraction Fraction of values to mask (0-1)
#' @param n_repeats Number of CV repetitions
#' @param methods Character vector of methods to test: "softimpute", "knn"
#' @param rank_max Passed to impute_softimpute
#' @param k Passed to impute_knn
#' @param path_stop_step Passed to impute_softimpute
#' @param mask_type Either "random_cells", "gene_dataset_block", or
#'   "gene_dataset_block_hard"
#' @return Data frame with validation metrics
#' @export
validate_imputation <- function(exprs_list,
                                 min_datasets = 1L,
                                 leave_out_fraction = 0.1,
                                 n_repeats = 5,
                                 methods = c("softimpute"),
                                 rank_max = 50,
                                 k = 10,
                                 path_stop_step = 7L,
                                 mask_type = "random_cells",
                                 min_obs_per_gene = 4L) {

  mask_type <- match.arg(mask_type, c("random_cells", "gene_dataset_block",
                                      "gene_dataset_block_hard",
                                      "progressive_tax_block"))

  cat("\n=== Imputation Validation ===\n")
  cat("Mask type:", mask_type, "\n")
  cat("Leave-out fraction:", leave_out_fraction, "\n")
  cat("Repeats:", n_repeats, "\n")
  cat("Methods:", paste(methods, collapse = ", "), "\n\n")

  # Create the base incomplete matrix
  incomplete <- create_incomplete_matrix(exprs_list, min_datasets)
  X <- incomplete$matrix

  # Get indices of observed (non-NA) values
  observed_idx <- which(!is.na(X))
  n_observed <- length(observed_idx)
  n_mask_target <- round(n_observed * leave_out_fraction)

  cat("Observed values:", n_observed, "\n")
  cat("Target masked cells per repeat:", n_mask_target, "\n")

  # Pre-compute (gene, dataset) pair pool for block masking
  obs_pairs <- NULL
  ds_cols_cache <- NULL
  row_idx_map <- NULL
  presence <- NULL
  n_ds_per_gene <- NULL
  n_datasets_total <- NULL
  is_block_mask <- mask_type %in% c("gene_dataset_block",
                                     "gene_dataset_block_hard",
                                     "progressive_tax_block")
  if (is_block_mask) {
    sample_info <- incomplete$sample_info
    datasets <- unique(sample_info$dataset)
    ds_cols_cache <- split(seq_len(ncol(X)), sample_info$dataset)
    row_idx_map <- setNames(seq_len(nrow(X)), rownames(X))

    # Boolean presence matrix: genes x datasets (TRUE if gene has >=1 obs)
    presence <- sapply(datasets, function(ds) {
      rowSums(!is.na(X[, ds_cols_cache[[ds]], drop = FALSE])) > 0
    })
    colnames(presence) <- datasets
    n_ds_per_gene <- rowSums(presence)
    n_datasets_total <- length(datasets)
  }

  if (mask_type %in% c("gene_dataset_block", "gene_dataset_block_hard")) {
    pair_list <- list()
    for (ds in datasets) {
      if (mask_type == "gene_dataset_block") {
        # Only genes present in ALL datasets — safe for biScale convergence
        eligible <- presence[, ds] & (n_ds_per_gene == n_datasets_total)
      } else {
        # Hard mode: any gene present in this dataset is eligible
        eligible <- presence[, ds]
      }
      genes_in_ds <- rownames(X)[eligible]
      if (length(genes_in_ds) > 0) {
        pair_list[[ds]] <- data.frame(
          gene = genes_in_ds,
          dataset = ds,
          stringsAsFactors = FALSE
        )
      }
    }
    obs_pairs <- do.call(rbind, pair_list)
    rownames(obs_pairs) <- NULL
    if (mask_type == "gene_dataset_block") {
      cat("Observed (gene, dataset) pairs (gene in all",
          n_datasets_total, "datasets):", nrow(obs_pairs), "\n")
    } else {
      cat("Observed (gene, dataset) pairs (hard, any coverage):",
          nrow(obs_pairs), "\n")
      cat("  Gene coverage distribution:\n")
      for (k_tier in sort(unique(n_ds_per_gene))) {
        n_genes_tier <- sum(n_ds_per_gene == k_tier)
        cat("    ", k_tier, "/", n_datasets_total, " datasets: ",
            n_genes_tier, " genes\n", sep = "")
      }
    }
  }

  if (mask_type == "progressive_tax_block") {
    n_min <- min(n_ds_per_gene)
    n_maskable <- sum(n_ds_per_gene > n_min)
    cat("Progressive-tax block masking:\n")
    cat("  n_min (poorest gene coverage):", n_min, "/", n_datasets_total, "\n")
    cat("  Genes at n_min (untouched):", sum(n_ds_per_gene == n_min), "\n")
    cat("  Genes above n_min (will be masked down to", n_min, "):",
        n_maskable, "\n")
    cat("  Gene coverage distribution:\n")
    for (k_tier in sort(unique(n_ds_per_gene))) {
      n_genes_tier <- sum(n_ds_per_gene == k_tier)
      blocks_to_mask <- max(0L, k_tier - n_min)
      cat("    ", k_tier, "/", n_datasets_total, " datasets: ",
          n_genes_tier, " genes -> mask ", blocks_to_mask, " blocks each\n",
          sep = "")
    }
  }
  cat("\n")

  results <- list()

  for (rep in seq_len(n_repeats)) {
    cat(sprintf("[%s] Repeat %d / %d ...\n",
                format(Sys.time(), "%H:%M:%S"), rep, n_repeats))

    set.seed(rep * 123)  # Reproducible

    if (mask_type == "random_cells") {
      m <- mask_random_cells(observed_idx, n_mask_target)
    } else if (mask_type == "progressive_tax_block") {
      m <- mask_progressive_tax_block(X, presence, n_ds_per_gene,
                                      ds_cols_cache, n_min, n_mask_target)
    } else {
      m <- mask_gene_dataset_block(obs_pairs, X, ds_cols_cache, row_idx_map,
                                   n_mask_target, min_obs_per_gene)
    }
    mask_idx <- m$mask_idx
    n_blocks <- m$n_blocks
    n_mask <- length(mask_idx)
    cat("  Cells masked:", n_mask, "\n")

    # Create masked matrix
    X_masked <- X
    true_values <- X[mask_idx]
    X_masked[mask_idx] <- NA

    # Create new incomplete object with masked data
    incomplete_masked <- incomplete
    incomplete_masked$matrix <- X_masked

    # Test each method (both score against the SAME mask_idx / true_values)
    for (method in methods) {
      method_cfg <- list(rank_max = rank_max, k = k, path_stop_step = path_stop_step)
      err_msg <- NA_character_
      imputed <- tryCatch(
        run_imputer(method, incomplete_masked, method_cfg),
        error = function(e) {
          cat("  ", method, " failed: ", conditionMessage(e), "\n", sep = "")
          err_msg <<- conditionMessage(e)
          NULL
        }
      )

      if (is.null(imputed)) {
        results[[length(results) + 1]] <- data.frame(
          method = method,
          mask_type = mask_type,
          repeat_n = rep,
          correlation = NA_real_,
          rmse = NA_real_,
          mae = NA_real_,
          n_masked = n_mask,
          n_blocks = n_blocks,
          converged = FALSE,
          error_message = err_msg,
          stringsAsFactors = FALSE
        )
        next
      }

      # Get imputed values at masked positions
      imputed_values <- imputed$matrix[mask_idx]

      # Calculate metrics
      correlation <- cor(true_values, imputed_values, use = "complete.obs")
      rmse <- sqrt(mean((true_values - imputed_values)^2, na.rm = TRUE))
      mae <- mean(abs(true_values - imputed_values), na.rm = TRUE)

      results[[length(results) + 1]] <- data.frame(
        method = method,
        mask_type = mask_type,
        repeat_n = rep,
        correlation = correlation,
        rmse = rmse,
        mae = mae,
        n_masked = n_mask,
        n_blocks = n_blocks,
        converged = TRUE,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  results_df <- do.call(rbind, results)

  # Summary statistics
  cat("\n=== Validation Summary ===\n")
  summary_df <- aggregate(
    cbind(correlation, rmse, mae) ~ method,
    data = results_df,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )

  for (method in unique(results_df$method)) {
    method_results <- results_df[results_df$method == method, ]
    n_total <- nrow(method_results)
    n_converged <- sum(method_results$converged, na.rm = TRUE)
    n_failed <- n_total - n_converged
    cat("\n", method, ":\n", sep = "")
    cat("  Converged: ", n_converged, "/", n_total, sep = "")
    if (n_failed > 0) {
      cat("  (", round(100 * n_failed / n_total, 1), "% failed)", sep = "")
    }
    cat("\n")
    converged_results <- method_results[method_results$converged, ]
    if (nrow(converged_results) > 0) {
      cat("  Correlation: ", round(mean(converged_results$correlation, na.rm = TRUE), 3),
          " (+/- ", round(sd(converged_results$correlation, na.rm = TRUE), 3), ")\n", sep = "")
      cat("  RMSE: ", round(mean(converged_results$rmse, na.rm = TRUE), 3),
          " (+/- ", round(sd(converged_results$rmse, na.rm = TRUE), 3), ")\n", sep = "")
      cat("  MAE: ", round(mean(converged_results$mae, na.rm = TRUE), 3),
          " (+/- ", round(sd(converged_results$mae, na.rm = TRUE), 3), ")\n", sep = "")
    }
  }

  results_df
}

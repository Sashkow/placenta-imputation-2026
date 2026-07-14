#!/usr/bin/env Rscript

library(yaml)
library(limma)

`%||%` <- function(a, b) if (is.null(a)) b else a

phase2b_dir <- "scripts/pipeline"
source(file.path(phase2b_dir, "imputation.R"))
source(file.path(phase2b_dir, "normalization.R"))

# ============================================================
# Data loading
# ============================================================

load_pipeline_data <- function(config) {
  cat("=== Loading pipeline data ===\n")

  pheno_path <- config$paths$phenodata
  if (grepl("\\.tsv$", pheno_path)) {
    phenodata <- read.delim(pheno_path, stringsAsFactors = FALSE)
  } else {
    phenodata <- read.csv(pheno_path, stringsAsFactors = FALSE)
  }
  if (!is.null(config$paths$column_map)) {
    for (to_col in names(config$paths$column_map)) {
      from_col <- config$paths$column_map[[to_col]]
      if (from_col %in% colnames(phenodata)) {
        phenodata[[to_col]] <- phenodata[[from_col]]
      }
    }
  }
  cat("Loaded phenodata:", nrow(phenodata), "samples\n")

  datasets <- config$files$datasets
  file_suffix <- config$files$suffix %||% ""
  file_map <- config$files$file_map
  exprs_list <- list()

  for (ds in datasets) {
    if (!is.null(file_map) && !is.null(file_map[[ds]])) {
      fp <- file.path(config$paths$mapped_data, file_map[[ds]])
    } else {
      fp <- file.path(config$paths$mapped_data, paste0(ds, file_suffix, ".tsv"))
    }
    if (file.exists(fp)) {
      exprs_list[[ds]] <- read.delim(fp, row.names = 1, check.names = FALSE)
      cat("  Loaded", ds, ":", nrow(exprs_list[[ds]]), "genes x",
          ncol(exprs_list[[ds]]), "samples\n")
    } else {
      cat("  WARNING: File not found:", fp, "\n")
    }
  }

  if (!is.null(config$sample_filter) && length(config$sample_filter) > 0) {
    cat("\nApplying sample filter:\n")
    mask <- rep(TRUE, nrow(phenodata))
    for (col in names(config$sample_filter)) {
      allowed <- config$sample_filter[[col]]
      if (!col %in% colnames(phenodata)) stop("sample_filter column not in phenodata: ", col)
      mask <- mask & phenodata[[col]] %in% allowed
      cat("  ", col, ":", paste(allowed, collapse = ", "), "\n")
    }
    allowed_samples <- phenodata$arraydatafile_exprscolumnnames[mask]
    for (ds in names(exprs_list)) {
      keep_cols <- intersect(colnames(exprs_list[[ds]]), allowed_samples)
      exprs_list[[ds]] <- exprs_list[[ds]][, keep_cols, drop = FALSE]
    }
    exprs_list <- exprs_list[sapply(exprs_list, ncol) > 0]
    phenodata <- phenodata[mask, ]
  }

  covariates <- config$phenotype$covariates %||% c()
  if (length(covariates) > 0) {
    cat("\nDropping samples with missing covariates:", paste(covariates, collapse = ", "), "\n")
    all_merged_samples <- unlist(lapply(exprs_list, colnames))
    pdata_merged <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% all_merged_samples, ]
    for (cov in covariates) {
      v <- pdata_merged[[cov]]
      v[v %in% c("", "_", "NA")] <- NA
      pdata_merged[[cov]] <- v
    }
    complete <- complete.cases(pdata_merged[, covariates, drop = FALSE])
    if (any(!complete)) {
      drop_ids <- pdata_merged$arraydatafile_exprscolumnnames[!complete]
      cat("  Dropping", length(drop_ids), "samples with missing covariates\n")
      for (ds in names(exprs_list)) {
        keep <- setdiff(colnames(exprs_list[[ds]]), drop_ids)
        exprs_list[[ds]] <- exprs_list[[ds]][, keep, drop = FALSE]
      }
      exprs_list <- exprs_list[sapply(exprs_list, ncol) > 0]
      phenodata <- phenodata[!phenodata$arraydatafile_exprscolumnnames %in% drop_ids, ]
    }
  }

  if (isTRUE(config$gene_filter$protein_coding_only)) {
    cat("\nFiltering to protein-coding genes...\n")
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("org.Hs.eg.db package required for protein_coding_only filter")
    gene_types <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "ENTREZID"),
      columns = c("ENTREZID", "GENETYPE"),
      keytype = "ENTREZID"
    )
    protein_coding <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding"]
    for (ds in names(exprs_list)) {
      keep <- rownames(exprs_list[[ds]]) %in% protein_coding
      exprs_list[[ds]] <- exprs_list[[ds]][keep, , drop = FALSE]
    }
    cat("  Protein-coding filter applied\n")
  }

  all_samples <- unlist(lapply(exprs_list, colnames))
  phenodata <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% all_samples, ]

  group_col <- config$phenotype$group_column
  first_trim <- phenodata$arraydatafile_exprscolumnnames[
    phenodata[[group_col]] == config$phenotype$baseline
  ]
  second_trim <- phenodata$arraydatafile_exprscolumnnames[
    phenodata[[group_col]] == config$phenotype$contrast
  ]

  ref_de <- read.delim(config$paths$reference_de, stringsAsFactors = FALSE)
  fdr_thresh <- config$thresholds$fdr %||% 0.05
  logfc_thresh <- config$thresholds$logfc %||% 1.0
  sig_mask <- !is.na(ref_de$adj.P.Val) & !is.na(ref_de$logFC) &
    ref_de$adj.P.Val < fdr_thresh & abs(ref_de$logFC) > logfc_thresh
  ref_sig_genes <- ref_de$gene[sig_mask]

  fdr_only_mask <- !is.na(ref_de$adj.P.Val) & ref_de$adj.P.Val < fdr_thresh
  ref_sig_genes_fdr_only <- ref_de$gene[fdr_only_mask]

  ref_fdr_valid <- !is.na(ref_de$adj.P.Val) &
    !is.na(ref_de$logFC) & ref_de$adj.P.Val < fdr_thresh
  ref_sig_genes_lfc15 <- ref_de$gene[
    ref_fdr_valid & abs(ref_de$logFC) > 1.5
  ]
  ref_sig_genes_lfc20 <- ref_de$gene[
    ref_fdr_valid & abs(ref_de$logFC) > 2.0
  ]

  ref_imputed <- NULL
  ref_imp_path <- config$paths$reference_imputed
  if (!is.null(ref_imp_path) && file.exists(ref_imp_path)) {
    ref_imputed <- as.matrix(read.delim(
      ref_imp_path, row.names = 1, check.names = FALSE
    ))
    cat("Loaded reference imputed matrix:",
        nrow(ref_imputed), "x", ncol(ref_imputed), "\n")
  }

  cat("\nData loaded:", length(exprs_list), "datasets,",
      length(first_trim), "baseline,", length(second_trim), "contrast samples\n")
  cat("Reference DEGs:", length(ref_sig_genes),
      " FDR-only:", length(ref_sig_genes_fdr_only),
      " lfc1.5:", length(ref_sig_genes_lfc15),
      " lfc2.0:", length(ref_sig_genes_lfc20), "\n")

  list(
    exprs_list = exprs_list,
    phenodata = phenodata,
    ref_de = ref_de,
    ref_sig_genes = ref_sig_genes,
    ref_sig_genes_fdr_only = ref_sig_genes_fdr_only,
    ref_sig_genes_lfc15 = ref_sig_genes_lfc15,
    ref_sig_genes_lfc20 = ref_sig_genes_lfc20,
    ref_imputed = ref_imputed,
    first_trim_samples = first_trim,
    second_trim_samples = second_trim
  )
}

# ============================================================
# Subsample helper
# ============================================================

subsample_exprs_list <- function(exprs_list, phenodata, sample_ids) {
  sub_list <- list()
  for (ds in names(exprs_list)) {
    keep <- intersect(colnames(exprs_list[[ds]]), sample_ids)
    if (length(keep) > 0) {
      sub_list[[ds]] <- exprs_list[[ds]][, keep, drop = FALSE]
    }
  }
  sub_pdata <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% sample_ids, ]
  list(exprs_list = sub_list, phenodata = sub_pdata)
}

# ============================================================
# Lean pipeline
# ============================================================

run_lean_pipeline <- function(exprs_sub_list, pdata_sub, config,
                              ref_imputed = NULL) {
  tryCatch(
    run_lean_pipeline_inner(exprs_sub_list, pdata_sub, config, ref_imputed),
    error = function(e) {
      list(
        de_results = NULL, sig_genes = character(0),
        n_sig = NA_integer_, n_up = NA_integer_, n_down = NA_integer_,
        n_genes_tested = NA_integer_,
        status = paste0("error:", conditionMessage(e))
      )
    }
  )
}

run_lean_pipeline_inner <- function(exprs_sub_list, pdata_sub, config,
                                    ref_imputed = NULL) {
  group_col <- config$phenotype$group_column
  baseline <- config$phenotype$baseline
  contrast <- config$phenotype$contrast
  covariates <- config$phenotype$covariates %||% c()

  # --- Drop imputed genes missing in one group ---
  if (isTRUE(config$coverage$drop_imputed_genes_missing_in_group)) {
    gene_sets <- lapply(exprs_sub_list, rownames)
    common_genes <- Reduce(intersect, gene_sets)
    imputed_genes <- setdiff(Reduce(union, gene_sets), common_genes)

    if (length(imputed_genes) > 0) {
      dataset_samples <- lapply(exprs_sub_list, colnames)
      sample_group <- setNames(
        pdata_sub[[group_col]],
        pdata_sub$arraydatafile_exprscolumnnames
      )

      # Vectorized: build presence matrix (imputed_genes x datasets)
      presence <- sapply(gene_sets, function(s) imputed_genes %in% s)
      if (!is.matrix(presence)) presence <- matrix(presence, nrow = length(imputed_genes))

      # For each imputed gene, count real values per group
      n_baseline <- integer(length(imputed_genes))
      n_contrast <- integer(length(imputed_genes))
      for (di in seq_along(dataset_samples)) {
        ds_grps <- sample_group[dataset_samples[[di]]]
        nb <- sum(ds_grps == baseline, na.rm = TRUE)
        nc <- sum(ds_grps == contrast, na.rm = TRUE)
        has <- presence[, di]
        n_baseline[has] <- n_baseline[has] + nb
        n_contrast[has] <- n_contrast[has] + nc
      }
      drop_genes <- imputed_genes[n_baseline == 0 | n_contrast == 0]

      if (length(drop_genes) > 0) {
        for (ds in names(exprs_sub_list)) {
          keep <- !(rownames(exprs_sub_list[[ds]]) %in% drop_genes)
          exprs_sub_list[[ds]] <- exprs_sub_list[[ds]][keep, , drop = FALSE]
        }
      }
    }
  }

  # --- Determine min_datasets (fixed or auto-scan) ---
  n_ds <- length(exprs_sub_list)
  fixed_min_ds <- config$coverage$min_datasets

  if (!is.null(fixed_min_ds)) {
    fixed_min_ds <- as.integer(fixed_min_ds)
    if (fixed_min_ds < 1L || fixed_min_ds > n_ds)
      stop("coverage$min_datasets must be between 1 and ", n_ds,
           " (got ", fixed_min_ds, ")")
    incomplete <- create_incomplete_matrix(exprs_sub_list,
                                           min_datasets = fixed_min_ds)
  } else {
    max_imp <- config$coverage$max_imputation_allowed %||% 0.20
    incomplete <- NULL
    for (min_ds in seq_len(n_ds)) {
      incomplete <- create_incomplete_matrix(exprs_sub_list,
                                             min_datasets = min_ds)
      frac <- sum(is.na(incomplete$matrix)) / length(incomplete$matrix)
      if (frac <= max_imp) break
      incomplete <- NULL
    }
    if (is.null(incomplete))
      stop("No min_datasets (1..", n_ds, ") keeps missingness below ",
           round(100 * max_imp, 1), "%")
  }

  # --- Imputation ---
  imputed_mask <- is.na(incomplete$matrix)
  dimnames(imputed_mask) <- dimnames(incomplete$matrix)
  use_ref <- FALSE
  if (!is.null(ref_imputed)) {
    needed_genes <- rownames(incomplete$matrix)
    needed_samples <- colnames(incomplete$matrix)
    if (all(needed_samples %in% colnames(ref_imputed)) &&
        all(needed_genes %in% rownames(ref_imputed))) {
      merged_exprs <- ref_imputed[needed_genes, needed_samples]
      use_ref <- TRUE
    }
  }
  if (!use_ref) {
    imputed <- run_imputer("softimpute", incomplete,
                           config$imputation$softimpute)
    merged_exprs <- imputed$matrix
  }

  # --- Build phenodata for merged samples ---
  merged_samples <- colnames(merged_exprs)
  merged_pdata <- pdata_sub[pdata_sub$arraydatafile_exprscolumnnames %in% merged_samples, ]
  rownames(merged_pdata) <- merged_pdata$arraydatafile_exprscolumnnames
  merged_pdata <- merged_pdata[merged_samples, ]

  # --- ComBat ---
  batch <- as.factor(merged_pdata$secondaryaccession)
  bio_group <- merged_pdata[[group_col]]
  combat_mod_data <- data.frame(bio_group = bio_group, stringsAsFactors = FALSE)
  for (cov in covariates) combat_mod_data[[cov]] <- merged_pdata[[cov]]
  combat_mod_rhs <- paste(c("bio_group", covariates), collapse = " + ")

  # Drop covariates with single level
  for (cov in covariates) {
    if (length(unique(na.omit(combat_mod_data[[cov]]))) < 2) {
      combat_mod_data[[cov]] <- NULL
      combat_mod_rhs <- paste(setdiff(c("bio_group", covariates), cov), collapse = " + ")
      covariates <- setdiff(covariates, cov)
    }
  }
  combat_mod <- model.matrix(as.formula(paste("~", combat_mod_rhs)),
                             data = combat_mod_data)

  ref_batch <- config$normalization$ref_batch
  if (!is.null(ref_batch) && ref_batch %in% levels(batch)) {
    normalized <- tryCatch(
      normalize_combat_ref(merged_exprs, batch, mod = combat_mod,
                           ref_batch = ref_batch),
      error = function(e) {
        normalize_combat(merged_exprs, batch, mod = combat_mod)
      }
    )
  } else {
    normalized <- normalize_combat(merged_exprs, batch, mod = combat_mod)
  }

  # --- limma DE ---
  keep_samples <- merged_pdata[[group_col]] %in% c(baseline, contrast)
  de_exprs <- normalized[, keep_samples]
  de_pdata <- merged_pdata[keep_samples, ]
  de_mask <- imputed_mask[rownames(de_exprs), keep_samples, drop = FALSE]

  group <- factor(de_pdata[[group_col]], levels = c(baseline, contrast))
  design_data <- data.frame(group = group, stringsAsFactors = FALSE)
  design_terms <- c("group")
  for (cov in covariates) {
    design_data[[cov]] <- de_pdata[[cov]]
    design_terms <- c(design_terms, cov)
  }
  design <- model.matrix(
    as.formula(paste("~", paste(design_terms, collapse = " + "))),
    data = design_data
  )

  imp_w <- config$de$imputed_cell_weight %||% 1.0
  if (imp_w != 1.0 && any(de_mask)) {
    W <- matrix(1.0, nrow(de_exprs), ncol(de_exprs), dimnames = dimnames(de_exprs))
    W[de_mask] <- imp_w
    fit <- lmFit(de_exprs, design, weights = W)
  } else {
    fit <- lmFit(de_exprs, design)
  }
  fit <- eBayes(fit)

  de_results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  de_results$gene <- rownames(de_results)
  de_results <- de_results[, c("gene", "logFC", "AveExpr", "t",
                                "P.Value", "adj.P.Val", "B")]

  fdr_thresh <- config$thresholds$fdr %||% 0.05
  logfc_thresh <- config$thresholds$logfc %||% 1.0
  sig_mask <- !is.na(de_results$adj.P.Val) & !is.na(de_results$logFC) &
    de_results$adj.P.Val < fdr_thresh & abs(de_results$logFC) > logfc_thresh
  sig_genes <- de_results$gene[sig_mask]

  fdr_only_mask <- !is.na(de_results$adj.P.Val) &
    de_results$adj.P.Val < fdr_thresh
  sig_genes_fdr_only <- de_results$gene[fdr_only_mask]

  fdr_valid <- !is.na(de_results$adj.P.Val) &
    !is.na(de_results$logFC) &
    de_results$adj.P.Val < fdr_thresh
  sig_genes_lfc15 <- de_results$gene[
    fdr_valid & abs(de_results$logFC) > 1.5
  ]
  sig_genes_lfc20 <- de_results$gene[
    fdr_valid & abs(de_results$logFC) > 2.0
  ]

  list(
    de_results = de_results,
    sig_genes = sig_genes,
    sig_genes_fdr_only = sig_genes_fdr_only,
    sig_genes_lfc15 = sig_genes_lfc15,
    sig_genes_lfc20 = sig_genes_lfc20,
    n_sig = length(sig_genes),
    n_sig_fdr_only = length(sig_genes_fdr_only),
    n_sig_lfc15 = length(sig_genes_lfc15),
    n_sig_lfc20 = length(sig_genes_lfc20),
    n_up = sum(de_results$logFC[sig_mask] > 0),
    n_down = sum(de_results$logFC[sig_mask] < 0),
    n_genes_tested = nrow(de_results),
    status = "ok"
  )
}

# ============================================================
# Metrics
# ============================================================

jaccard <- function(a, b) {
  int <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  int / uni
}

lin_ccc <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  mx <- mean(x); my <- mean(y)
  sx <- var(x); sy <- var(y)
  sxy <- cov(x, y)
  2 * sxy / (sx + sy + (mx - my)^2)
}

compute_subsample_metrics <- function(result, ref_de, ref_sig_genes,
                                      ref_sig_genes_fdr_only = NULL) {
  na_row <- data.frame(
    n_deg = NA, n_up = NA, n_down = NA, n_genes_tested = NA,
    jaccard_vs_full = NA, overlap_vs_full = NA,
    logfc_pearson_vs_full = NA, logfc_ccc_vs_full = NA,
    logfc_spearman_vs_full = NA,
    same_direction_pct = NA,
    n_deg_fdr_only = NA, jaccard_fdr_only = NA,
    overlap_fdr_only = NA,
    stringsAsFactors = FALSE
  )
  if (result$status != "ok" || is.null(result$de_results))
    return(na_row)

  de <- result$de_results
  sig <- result$sig_genes
  j <- jaccard(sig, ref_sig_genes)
  overlap <- length(intersect(sig, ref_sig_genes)) /
    max(length(ref_sig_genes), 1)

  shared_genes <- intersect(de$gene, ref_de$gene)
  if (length(shared_genes) > 1) {
    de_sub <- de[match(shared_genes, de$gene), ]
    ref_sub <- ref_de[match(shared_genes, ref_de$gene), ]
    pearson_r <- cor(de_sub$logFC, ref_sub$logFC,
                     use = "complete.obs")
    ccc_r <- lin_ccc(de_sub$logFC, ref_sub$logFC)
    spearman_r <- cor(de_sub$logFC, ref_sub$logFC,
                      use = "complete.obs", method = "spearman")
    shared_sig <- intersect(sig, ref_sig_genes)
    if (length(shared_sig) > 0) {
      de_dir <- sign(de_sub$logFC[de_sub$gene %in% shared_sig])
      ref_dir <- sign(ref_sub$logFC[ref_sub$gene %in% shared_sig])
      same_dir <- sum(de_dir == ref_dir, na.rm = TRUE) /
        max(length(shared_sig), 1) * 100
    } else {
      same_dir <- NA_real_
    }
  } else {
    pearson_r <- NA_real_
    ccc_r <- NA_real_
    spearman_r <- NA_real_
    same_dir <- NA_real_
  }

  sig_fdr <- result$sig_genes_fdr_only
  n_fdr <- length(sig_fdr)
  j_fdr <- NA_real_
  ov_fdr <- NA_real_
  if (!is.null(ref_sig_genes_fdr_only) && length(sig_fdr) > 0) {
    j_fdr <- jaccard(sig_fdr, ref_sig_genes_fdr_only)
    ov_fdr <- length(intersect(sig_fdr, ref_sig_genes_fdr_only)) /
      max(length(ref_sig_genes_fdr_only), 1)
  }

  data.frame(
    n_deg = result$n_sig, n_up = result$n_up,
    n_down = result$n_down,
    n_genes_tested = result$n_genes_tested,
    jaccard_vs_full = j, overlap_vs_full = overlap,
    logfc_pearson_vs_full = pearson_r,
    logfc_ccc_vs_full = ccc_r,
    logfc_spearman_vs_full = spearman_r,
    same_direction_pct = same_dir,
    n_deg_fdr_only = n_fdr, jaccard_fdr_only = j_fdr,
    overlap_fdr_only = ov_fdr,
    stringsAsFactors = FALSE
  )
}

compute_pairwise_metrics <- function(result_a, result_b) {
  na_row <- data.frame(
    jaccard_between = NA, logfc_pearson_between = NA,
    logfc_ccc_between = NA,
    logfc_spearman_between = NA, same_direction_pct = NA,
    jaccard_fdr_only_between = NA,
    jaccard_lfc15_between = NA,
    jaccard_lfc20_between = NA,
    stringsAsFactors = FALSE
  )
  if (result_a$status != "ok" || result_b$status != "ok")
    return(na_row)

  j <- jaccard(result_a$sig_genes, result_b$sig_genes)
  j_fdr <- jaccard(result_a$sig_genes_fdr_only,
                    result_b$sig_genes_fdr_only)
  j_lfc15 <- jaccard(result_a$sig_genes_lfc15,
                      result_b$sig_genes_lfc15)
  j_lfc20 <- jaccard(result_a$sig_genes_lfc20,
                      result_b$sig_genes_lfc20)
  shared <- intersect(result_a$de_results$gene,
                      result_b$de_results$gene)
  if (length(shared) > 1) {
    a_sub <- result_a$de_results[
      match(shared, result_a$de_results$gene), ]
    b_sub <- result_b$de_results[
      match(shared, result_b$de_results$gene), ]
    pr <- cor(a_sub$logFC, b_sub$logFC,
              use = "complete.obs")
    cc <- lin_ccc(a_sub$logFC, b_sub$logFC)
    sr <- cor(a_sub$logFC, b_sub$logFC,
              use = "complete.obs", method = "spearman")
    shared_sig <- intersect(
      result_a$sig_genes, result_b$sig_genes)
    if (length(shared_sig) > 0) {
      da <- sign(a_sub$logFC[a_sub$gene %in% shared_sig])
      db <- sign(b_sub$logFC[b_sub$gene %in% shared_sig])
      same_dir <- sum(da == db, na.rm = TRUE) /
        max(length(shared_sig), 1) * 100
    } else {
      same_dir <- NA_real_
    }
  } else {
    pr <- NA_real_
    cc <- NA_real_
    sr <- NA_real_
    same_dir <- NA_real_
  }

  data.frame(
    jaccard_between = j, logfc_pearson_between = pr,
    logfc_ccc_between = cc,
    logfc_spearman_between = sr, same_direction_pct = same_dir,
    jaccard_fdr_only_between = j_fdr,
    jaccard_lfc15_between = j_lfc15,
    jaccard_lfc20_between = j_lfc20,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# CLI config loader
# ============================================================

parse_config_arg <- function(default_config = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- default_config
  for (arg in args) {
    if (grepl("^--config=", arg)) {
      config_path <- sub("^--config=", "", arg)
    }
  }
  if (is.null(config_path))
    stop("No config specified. Use --config=PATH")
  yaml::read_yaml(config_path)
}

parse_int_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  pat <- paste0("^--", name, "=")
  for (arg in args) {
    if (grepl(pat, arg)) return(as.integer(sub(pat, "", arg)))
  }
  default
}

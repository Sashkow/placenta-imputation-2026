#!/usr/bin/env Rscript

#' Run Phase 2B: Direct Merge Improvements
#'
#' This phase evaluates imputation and alternative normalization methods
#' to improve gene recovery when directly merging expression datasets.
#'
#' Usage:
#'   Rscript run_phase2b.R [options]
#'
#' Options:
#'   --config=PATH         Path to config YAML (default: config_phase2b.yaml)
#'   --imputation=METHOD   Run only specific imputation: softimpute, knn, none
#'   --normalization=METHOD Run only specific normalization: combat, dwd, mean_center
#'   --validate_only       Only run imputation validation (no DE analysis)
#'   --output_dir=PATH     Override output directory
#'   --no_archive          Don't archive previous results
#'   --help                Show this help message

library(yaml)
library(limma)

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Get script directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return("scripts/integrative_analysis/phase2b_direct_merge")
}

script_dir <- get_script_dir()

# Source utilities
utils_path <- "scripts/utils/logging_utils.R"
if (!file.exists(utils_path)) {
  utils_path <- file.path("..", "..", "utils", "logging_utils.R")
}
if (file.exists(utils_path)) {
  source(utils_path)
} else {
  # Fallback logging functions
  setup_logging <- function(dir, prefix = "log") {
    log_file <- file.path(dir, paste0(prefix, "_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".txt"))
    cat("Logging to:", log_file, "\n")
    log_file
  }
  close_logging <- function(log_file) invisible(NULL)
  archive_previous_results <- function(output_dir) {
    if (dir.exists(output_dir) && length(list.files(output_dir)) > 0) {
      archive_dir <- file.path(output_dir, "archive",
                               format(Sys.time(), "%Y-%m-%d_%H%M%S"))
      dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
      files <- list.files(output_dir, full.names = TRUE)
      files <- files[!grepl("archive", files)]
      if (length(files) > 0) {
        file.copy(files, archive_dir, recursive = TRUE)
        unlink(files, recursive = TRUE)
        cat("Archived previous results to:", archive_dir, "\n")
      }
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
}

# Source module files
source(file.path(script_dir, "imputation.R"))
source(file.path(script_dir, "normalization.R"))
source(file.path(script_dir, "plot_na_staircase.R"))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Defaults
config_path <- file.path(script_dir, "config_phase2b.yaml")
override_imputation <- NULL
override_normalization <- NULL
validate_only <- FALSE
override_output <- NULL
no_archive <- FALSE

# Parse arguments
for (arg in args) {
  if (arg == "--help" || arg == "-h") {
    cat("
Phase 2B: Direct Merge Improvements

Usage:
  Rscript run_phase2b.R [options]

Options:
  --config=PATH           Path to config YAML (default: config_phase2b.yaml)
  --imputation=METHOD     Run only specific imputation: softimpute, knn, none
  --normalization=METHOD  Run only specific normalization: combat, dwd, mean_center
  --validate_only         Only run imputation validation (no DE analysis)
  --output_dir=PATH       Override output directory
  --no_archive            Don't archive previous results
  --help                  Show this help message

Examples:
  Rscript run_phase2b.R
  Rscript run_phase2b.R --imputation=softimpute --normalization=dwd
  Rscript run_phase2b.R --validate_only
")
    quit(status = 0)
  } else if (grepl("^--config=", arg)) {
    config_path <- sub("^--config=", "", arg)
  } else if (grepl("^--imputation=", arg)) {
    override_imputation <- sub("^--imputation=", "", arg)
  } else if (grepl("^--normalization=", arg)) {
    override_normalization <- sub("^--normalization=", "", arg)
  } else if (arg == "--validate_only") {
    validate_only <- TRUE
  } else if (grepl("^--output_dir=", arg)) {
    override_output <- sub("^--output_dir=", "", arg)
  } else if (arg == "--no_archive") {
    no_archive <- TRUE
  }
}

# Load config
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}
config <- yaml::read_yaml(config_path)

# Compute numeric gestational age (weeks) from phenodata columns
compute_ga_weeks <- function(pdata) {
  ga <- numeric(nrow(pdata))
  for (i in seq_len(nrow(pdata))) {
    fw <- suppressWarnings(as.numeric(pdata$fetus_week[i]))
    if (!is.na(fw) && fw > 0) {
      ga[i] <- fw
    } else {
      rw <- as.character(pdata$fetus_range_week[i])
      if (!is.na(rw) && nchar(rw) > 0 && rw != "0") {
        weeks <- as.numeric(strsplit(rw, ",")[[1]])
        ga[i] <- mean(weeks, na.rm = TRUE)
      } else {
        tri <- as.character(pdata[[config$phenotype$group_column]][i])
        ga[i] <- switch(tri,
          "First Trimester" = 8, "First trimester" = 8,
          "Second Trimester" = 16, "Second trimester" = 16,
          "Third Trimester" = 39, "Third trimester" = 39,
          12)
      }
    }
  }
  ga
}

# Apply overrides
if (!is.null(override_output)) config$paths$output <- override_output

output_dir <- config$paths$output

cat("\n")
cat("============================================================\n")
cat("  Phase 2B: Direct Merge Improvements\n")
cat("============================================================\n\n")

cat("Configuration:\n")
cat("  Config file:     ", config_path, "\n")
cat("  Mapped path:     ", config$paths$mapped_data, "\n")
cat("  Phenodata:       ", config$paths$phenodata, "\n")
cat("  Output dir:      ", output_dir, "\n")
cat("  Comparison:      ", config$phenotype$contrast, "vs", config$phenotype$baseline, "\n")
if (!is.null(config$coverage$min_datasets)) {
  cat("  Coverage mode:    fixed min_datasets =", config$coverage$min_datasets, "\n")
} else {
  cat("  Coverage mode:    max_imputation_allowed =",
      config$coverage$max_imputation_allowed %||% 0.20, "\n")
}
cat("\n")

# Archive previous results
if (!no_archive) {
  archive_previous_results(output_dir)
} else {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Setup logging
log_file <- setup_logging(output_dir, prefix = "phase2b_log")

# Make sure any uncaught error lands in the log file before the script
# terminates. Without this, R's default error handler writes to stderr
# only, and the sink buffer may be lost on crash.
options(error = quote({
  tryCatch({
    if (sink.number() > 0) {
      cat("\n\n=== UNCAUGHT ERROR - script terminating ===\n")
      cat("Time:  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
      cat("Error: ", geterrmessage(), sep = "")
      cat("\nTraceback:\n")
      tb <- .traceback(2)
      if (length(tb)) {
        for (i in seq_along(tb)) {
          cat(sprintf("  %d: %s\n", i, paste(tb[[i]], collapse = " ")))
        }
      }
      while (sink.number() > 0) sink()
    }
  }, error = function(e) {})
  if (!interactive()) quit(save = "no", status = 1)
}))

# ============================================================
# Step 1: Load datasets
# ============================================================

cat("=== Step 1: Loading Datasets ===\n\n")

# Load phenodata
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
      cat("  Column mapped:", from_col, "->", to_col, "\n")
    }
  }
}
cat("Loaded phenodata:", nrow(phenodata), "samples\n")

# Load expression files
datasets <- config$files$datasets
exprs_list <- list()

file_suffix <- if (!is.null(config$files$suffix)) config$files$suffix else ""
file_map <- config$files$file_map
for (ds in datasets) {
  if (!is.null(file_map) && !is.null(file_map[[ds]])) {
    file_path <- file.path(config$paths$mapped_data, file_map[[ds]])
  } else {
    file_path <- file.path(config$paths$mapped_data, paste0(ds, file_suffix, ".tsv"))
  }
  if (file.exists(file_path)) {
    exprs <- read.delim(file_path, row.names = 1, check.names = FALSE)
    exprs_list[[ds]] <- exprs
    cat("  Loaded", ds, ":", nrow(exprs), "genes x", ncol(exprs), "samples\n")
  } else {
    cat("  WARNING: File not found:", file_path, "\n")
  }
}

cat("\nLoaded", length(exprs_list), "datasets\n")

# Apply sample filter (e.g. Diagnosis, Biological.Specimen, GA.Category).
# Restricts ComBat/imputation to samples matching phenodata criteria.
if (!is.null(config$sample_filter) && length(config$sample_filter) > 0) {
  cat("\nApplying sample filter:\n")
  mask <- rep(TRUE, nrow(phenodata))
  for (col in names(config$sample_filter)) {
    allowed <- config$sample_filter[[col]]
    cat("  ", col, ":", paste(allowed, collapse = ", "), "\n")
    if (!col %in% colnames(phenodata)) {
      stop("sample_filter column not in phenodata: ", col)
    }
    mask <- mask & phenodata[[col]] %in% allowed
  }
  allowed_samples <- phenodata$arraydatafile_exprscolumnnames[mask]
  cat("  ", length(allowed_samples), "samples pass filter globally\n")

  for (ds in names(exprs_list)) {
    before <- ncol(exprs_list[[ds]])
    keep_cols <- intersect(colnames(exprs_list[[ds]]), allowed_samples)
    exprs_list[[ds]] <- exprs_list[[ds]][, keep_cols, drop = FALSE]
    cat("  ", ds, ":", before, "->", ncol(exprs_list[[ds]]), "samples\n")
  }

  exprs_list <- exprs_list[sapply(exprs_list, ncol) > 0]
  cat("  ", length(exprs_list), "datasets retained after filter\n")
}

if (!is.null(config$per_dataset_filter)) {
  cat("\nApplying per-dataset filters:\n")
  for (ds in names(config$per_dataset_filter)) {
    if (!ds %in% names(exprs_list)) next
    ds_filter <- config$per_dataset_filter[[ds]]
    ds_samples <- colnames(exprs_list[[ds]])
    ds_pdata <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% ds_samples, ]
    ds_mask <- rep(TRUE, nrow(ds_pdata))
    for (col in names(ds_filter)) {
      vals <- ds_filter[[col]]
      ds_mask <- ds_mask & ds_pdata[[col]] %in% vals
      cat("  ", ds, col, ":", paste(vals, collapse = ", "), "\n")
    }
    keep <- ds_pdata$arraydatafile_exprscolumnnames[ds_mask]
    cat("  ", ds, ":", length(ds_samples), "->", length(keep), "samples\n")
    exprs_list[[ds]] <- exprs_list[[ds]][, keep, drop = FALSE]
  }
  exprs_list <- exprs_list[sapply(exprs_list, ncol) > 0]
}

# Filter to protein-coding genes (ENTREZIDs with GENETYPE == "protein-coding")
if (isTRUE(config$gene_filter$protein_coding_only)) {
  cat("\nFiltering to protein-coding genes (org.Hs.eg.db GENETYPE)...\n")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db package required for protein_coding_only filter")
  }
  gene_types <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "ENTREZID"),
    columns = c("ENTREZID", "GENETYPE"),
    keytype = "ENTREZID"
  )
  protein_coding <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding"]
  cat("  protein-coding reference set:", length(protein_coding), "genes\n")

  for (ds in names(exprs_list)) {
    before <- nrow(exprs_list[[ds]])
    keep <- rownames(exprs_list[[ds]]) %in% protein_coding
    exprs_list[[ds]] <- exprs_list[[ds]][keep, , drop = FALSE]
    cat("  ", ds, ":", before, "->", nrow(exprs_list[[ds]]), "genes\n")
  }
}

# ============================================================
# Step 2: Gene Coverage Analysis
# ============================================================

cat("\n=== Step 2: Gene Coverage Analysis ===\n")

# compare_gene_recovery() intentionally disabled: it ran softimpute on every
# threshold up-front (10x full imputations) just to build a reference CSV that
# isn't consumed downstream. The auto-selection scan below picks the threshold
# on its own; set gene_recovery to NULL so the summary block can skip it.
gene_recovery <- NULL

# ----------------------------------------------------------------
# For genes that will be imputed in the outer-join merged matrix,
# check whether they have non-imputed values in BOTH comparison
# groups. If a gene is missing from one group entirely, softimpute
# has to synthesize its expression under that condition, which
# could bias logFC estimates.
# ----------------------------------------------------------------

cat("\n=== Imputed-gene coverage per comparison group ===\n")

gene_sets <- lapply(exprs_list, rownames)
all_genes <- Reduce(union, gene_sets)
common_genes <- Reduce(intersect, gene_sets)
imputed_genes <- setdiff(all_genes, common_genes)

imputed_coverage <- list(
  n_common = length(common_genes),
  n_imputed = length(imputed_genes),
  n_both = NA_integer_,
  n_miss_baseline = NA_integer_,
  n_miss_contrast = NA_integer_,
  baseline_grp = config$phenotype$baseline,
  contrast_grp = config$phenotype$contrast,
  lines = character(0)
)

add_line <- function(s) {
  imputed_coverage$lines <<- c(imputed_coverage$lines, s)
  cat(s, "\n", sep = "")
}

add_line(sprintf("Genes in all datasets (no imputation): %d",
                 length(common_genes)))
add_line(sprintf("Genes needing imputation somewhere:    %d",
                 length(imputed_genes)))

if (length(imputed_genes) > 0) {
  baseline_grp <- config$phenotype$baseline
  contrast_grp <- config$phenotype$contrast

  # Per-dataset lookups: which samples belong to which group
  dataset_samples <- lapply(exprs_list, colnames)
  sample_group <- setNames(
    phenodata[[config$phenotype$group_column]],
    phenodata$arraydatafile_exprscolumnnames
  )

  non_imp_baseline <- integer(length(imputed_genes))
  non_imp_contrast <- integer(length(imputed_genes))

  for (i in seq_along(imputed_genes)) {
    g <- imputed_genes[i]
    has_gene <- vapply(gene_sets, function(s) g %in% s, logical(1))
    samples_with_value <- unlist(dataset_samples[has_gene], use.names = FALSE)
    grps <- sample_group[samples_with_value]
    non_imp_baseline[i] <- sum(grps == baseline_grp, na.rm = TRUE)
    non_imp_contrast[i] <- sum(grps == contrast_grp, na.rm = TRUE)
  }

  n_both <- sum(non_imp_baseline >= 1 & non_imp_contrast >= 1)
  n_miss_baseline <- sum(non_imp_baseline == 0)
  n_miss_contrast <- sum(non_imp_contrast == 0)

  imputed_coverage$n_both <- n_both
  imputed_coverage$n_miss_baseline <- n_miss_baseline
  imputed_coverage$n_miss_contrast <- n_miss_contrast

  add_line(sprintf("  In BOTH groups with >=1 real value: %d / %d",
                   n_both, length(imputed_genes)))
  add_line(sprintf("  Missing entirely from %-16s: %d",
                   baseline_grp, n_miss_baseline))
  add_line(sprintf("  Missing entirely from %-16s: %d",
                   contrast_grp, n_miss_contrast))

  coverage_df <- data.frame(
    gene = imputed_genes,
    n_datasets_with_gene = vapply(imputed_genes, function(g)
      sum(vapply(gene_sets, function(s) g %in% s, logical(1))),
      integer(1)),
    non_imputed_baseline = non_imp_baseline,
    non_imputed_contrast = non_imp_contrast,
    in_both_groups = (non_imp_baseline >= 1 & non_imp_contrast >= 1),
    stringsAsFactors = FALSE
  )
  colnames(coverage_df)[3:4] <- c(
    paste0("non_imputed_", make.names(baseline_grp)),
    paste0("non_imputed_", make.names(contrast_grp))
  )
  write.csv(coverage_df,
            file.path(output_dir, "imputed_gene_group_coverage.csv"),
            row.names = FALSE)
  cat("Saved: imputed_gene_group_coverage.csv\n")

  if (n_miss_baseline > 0 || n_miss_contrast > 0) {
    add_line(paste("WARNING: some imputed genes lack real values in one",
                   "comparison group - softimpute will synthesize their",
                   "expression under that condition."))
  }

  # Optionally drop imputed genes that are missing from one group entirely.
  # Flag: coverage.drop_imputed_genes_missing_in_group (default TRUE).
  drop_bad <- config$coverage$drop_imputed_genes_missing_in_group %||% TRUE
  bad_idx <- non_imp_baseline == 0 | non_imp_contrast == 0
  dropped_group_genes <- character(0)
  # Snapshot for staircase (before dropping genes)
  exprs_list_pre_drop <- lapply(exprs_list, function(e) e)
  if (drop_bad && any(bad_idx)) {
    dropped_group_genes <- imputed_genes[bad_idx]
    add_line(sprintf(
      "Dropping %d imputed genes without real values in both groups.",
      length(dropped_group_genes)))
    for (ds in names(exprs_list)) {
      keep <- !(rownames(exprs_list[[ds]]) %in% dropped_group_genes)
      exprs_list[[ds]] <- exprs_list[[ds]][keep, , drop = FALSE]
    }
  } else if (!drop_bad && any(bad_idx)) {
    add_line(sprintf(
      "Keeping %d imputed genes without real values in both groups (flag off).",
      sum(bad_idx)))
  }
}

# ----------------------------------------------------------------
# Determine min_datasets: the minimum number of datasets a gene
# must appear in to be included in the merged matrix.
#
# Two mutually exclusive config modes:
#   coverage$min_datasets  – use this value directly (no scan)
#   coverage$max_imputation_allowed – scan 1..N and pick the lowest
#       min_datasets whose missing fraction does not exceed the limit
#
# Runs after dropping group-imbalanced imputed genes so the missing
# fraction reflects what softimpute will actually see.
# ----------------------------------------------------------------

n_datasets_total <- length(exprs_list)
fixed_min_ds     <- config$coverage$min_datasets
max_imp_allowed  <- config$coverage$max_imputation_allowed

cat("\n=== Coverage threshold selection ===\n")
cat("Number of datasets:", n_datasets_total, "\n")

if (!is.null(fixed_min_ds)) {
  # --- Mode 1: fixed min_datasets from config ---
  fixed_min_ds <- as.integer(fixed_min_ds)
  if (fixed_min_ds < 1L || fixed_min_ds > n_datasets_total)
    stop("coverage$min_datasets must be between 1 and ", n_datasets_total,
         " (got ", fixed_min_ds, ")")
  selected_threshold <- fixed_min_ds
  inc <- create_incomplete_matrix(exprs_list, min_datasets = selected_threshold)
  frac <- sum(is.na(inc$matrix)) / length(inc$matrix)
  cat(sprintf("  Fixed min_datasets=%d/%d -> %d genes, %.2f%% missing\n",
              selected_threshold, n_datasets_total,
              nrow(inc$matrix), 100 * frac))
} else {
  # --- Mode 2: auto-select via max_imputation_allowed ---
  max_imp_allowed <- max_imp_allowed %||% 0.20
  cat("Max imputation allowed:", round(100 * max_imp_allowed, 1), "%\n")

  selected_threshold <- NA_integer_
  for (min_ds in seq_len(n_datasets_total)) {
    inc <- create_incomplete_matrix(exprs_list, min_datasets = min_ds)
    frac <- sum(is.na(inc$matrix)) / length(inc$matrix)
    cat(sprintf("  min_datasets=%d/%d -> %d genes, %.2f%% missing\n",
                min_ds, n_datasets_total, nrow(inc$matrix), 100 * frac))
    if (frac <= max_imp_allowed) {
      selected_threshold <- min_ds
      break
    }
  }

  if (is.na(selected_threshold)) {
    stop("No min_datasets value (1..", n_datasets_total,
         ") keeps missing fraction below ",
         round(100 * max_imp_allowed, 1), "%.")
  }
}

cat("Selected min_datasets:", selected_threshold, "of", n_datasets_total, "\n")
config$coverage$threshold <- selected_threshold

# ============================================================
# Step 3: Imputation Validation (if enabled)
# ============================================================

validation_results <- NULL

if (config$validation$leave_out_fraction %||% 0 > 0) {
  cat("\n=== Step 3: Imputation Validation ===\n")

  # Validate any enabled imputer that is in the IMPUTERS registry.
  methods_to_validate <- intersect(
    names(IMPUTERS),
    names(Filter(function(cfg) isTRUE(cfg$enabled),
                 config$imputation))
  )

  if (length(methods_to_validate) > 0) {
    validation_results <- validate_imputation(
      exprs_list,
      min_datasets = config$coverage$threshold,
      leave_out_fraction = config$validation$leave_out_fraction %||% 0.1,
      n_repeats = config$validation$n_repeats %||% 5,
      methods = methods_to_validate,
      rank_max = config$imputation$softimpute$rank_max %||% 50,
      k = config$imputation$knn$k %||% 10,
      mask_type = config$validation$mask_type %||% "random_cells",
      min_obs_per_gene = config$validation$min_obs_per_gene %||% 4L
    )

    # Save validation results
    write.csv(validation_results,
              file.path(output_dir, "imputation_validation.csv"),
              row.names = FALSE)
    cat("\nSaved imputation validation results\n")
  }
}

if (validate_only) {
  cat("\n--validate_only flag set, stopping here.\n")
  close_logging(log_file)
  quit(status = 0)
}

# ============================================================
# Step 4: Create Merged Matrices with Different Methods
# ============================================================

cat("\n=== Step 4: Creating Merged Expression Matrices ===\n")

# Methods to run: "none" plus every enabled IMPUTERS entry.
imputation_methods <- c()
if (config$imputation$none$enabled %||% TRUE) imputation_methods <- c(imputation_methods, "none")
enabled_imputers <- intersect(
  names(IMPUTERS),
  names(Filter(function(cfg) isTRUE(cfg$enabled), config$imputation))
)
imputation_methods <- c(imputation_methods, enabled_imputers)

if (!is.null(override_imputation)) {
  imputation_methods <- override_imputation
}

normalization_methods <- c()
if (config[["normalization"]][["combat"]][["enabled"]] %||% TRUE) normalization_methods <- c(normalization_methods, "combat")
if (config[["normalization"]][["batch_in_limma"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "batch_in_limma")
if (config[["normalization"]][["ruv"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "ruv")
if (config[["normalization"]][["ruvinv"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "ruvinv")
if (config[["normalization"]][["bruv"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "bruv")
if (config[["normalization"]][["combat_ref"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "combat_ref")
if (config[["normalization"]][["harmonizr"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "harmonizr")
if (config[["normalization"]][["github_harmonizr"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "github_harmonizr")
# DWD disabled: normalize_dwd() fails on single-sample batches (see 2_3 run).
# if (config[["normalization"]][["dwd"]][["enabled"]] %||% FALSE) normalization_methods <- c(normalization_methods, "dwd")

if (!is.null(override_normalization)) {
  normalization_methods <- override_normalization
}

ruv_control_genes <- NULL
if (any(c("ruv", "ruvinv", "bruv") %in% normalization_methods)) {
  ruv_cfg <- config[["normalization"]][["ruv"]]
  if (is.null(ruv_cfg)) ruv_cfg <- config[["normalization"]][["ruvinv"]]
  if (is.null(ruv_cfg)) ruv_cfg <- config[["normalization"]][["bruv"]]
  ruv_file <- ruv_cfg[["control_genes_file"]]
  if (!is.null(ruv_file) && file.exists(ruv_file)) {
    ruv_control_genes <- readLines(ruv_file)
    cat(sprintf("RUV control genes loaded: %d from %s\n",
                length(ruv_control_genes), ruv_file))
  } else {
    stop("RUV enabled but control_genes_file not found: ", ruv_file)
  }
}

# Create incomplete matrix
incomplete <- create_incomplete_matrix(
  exprs_list,
  min_datasets = config$coverage$threshold
)

# Store all results for comparison
all_results <- list()

for (imp_method in imputation_methods) {
  cat("\n--- Imputation:", imp_method, "---\n")

  # Apply imputation. "none" is a special case (inner join, no
  # imputation); any other method name is dispatched through the
  # IMPUTERS registry in imputation.R.
  if (imp_method == "none") {
    common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
    cat("Common genes across all datasets:", length(common_genes), "\n")

    # Use unname to prevent list names from being prepended to column names
    merged_exprs <- do.call(cbind, unname(lapply(exprs_list, function(e) e[common_genes, ])))
    imputed <- list(matrix = merged_exprs, method = "none")
    imputed_mask <- matrix(FALSE, nrow(merged_exprs), ncol(merged_exprs),
                           dimnames = dimnames(merged_exprs))
  } else {
    imputed <- tryCatch(
      run_imputer(
        imp_method,
        incomplete,
        config$imputation[[imp_method]]
      ),
      error = function(e) {
        cat("\n!!! Imputation method '", imp_method,
            "' failed: ", conditionMessage(e), "\n", sep = "")
        crash_file <- file.path(
          output_dir,
          paste0("crash_state_", imp_method, ".rds")
        )
        saveRDS(
          list(
            imp_method         = imp_method,
            error_message      = conditionMessage(e),
            error              = e,
            coverage_threshold = config$coverage$threshold,
            incomplete         = incomplete,
            exprs_list         = exprs_list,
            session_info       = sessionInfo(),
            time               = Sys.time()
          ),
          crash_file
        )
        cat("Saved crash snapshot: ", crash_file, "\n", sep = "")
        NULL
      }
    )
    if (is.null(imputed)) {
      cat("Skipping '", imp_method, "' for downstream steps.\n", sep = "")
      next
    }
    merged_exprs <- imputed$matrix
    # Imputed cell = cell that was NA in the pre-imputation incomplete matrix.
    imputed_mask <- is.na(incomplete$matrix)
    dimnames(imputed_mask) <- dimnames(incomplete$matrix)
  }

  cat(sprintf("Imputed cells in merged matrix: %d / %d (%.2f%%)\n",
              sum(imputed_mask), length(imputed_mask),
              100 * sum(imputed_mask) / length(imputed_mask)))

  # Save imputed matrix
  if (config$output$save_imputed_matrices %||% FALSE) {
    imputed_file <- file.path(output_dir, paste0("exprs_imputed_", imp_method, ".tsv"))
    write.table(merged_exprs, imputed_file, sep = "\t", quote = FALSE)
    cat("Saved imputed matrix:", imputed_file, "\n")
  }

  # Create phenodata for merged samples
  merged_samples <- colnames(merged_exprs)
  merged_pdata <- phenodata[phenodata$arraydatafile_exprscolumnnames %in% merged_samples, ]

  # Match order
  rownames(merged_pdata) <- merged_pdata$arraydatafile_exprscolumnnames
  merged_pdata <- merged_pdata[merged_samples, ]

  # Covariates (e.g. Combined.Fetus.Sex) — used for both ComBat mod and limma
  covariates <- config$phenotype$covariates %||% c()
  if (length(covariates) > 0) {
    for (cov in covariates) {
      if (!cov %in% colnames(merged_pdata)) {
        stop("phenotype.covariate not in phenodata: ", cov)
      }
      v <- merged_pdata[[cov]]
      v[v %in% c("", "_", "NA")] <- NA
      merged_pdata[[cov]] <- v
    }
    complete <- complete.cases(merged_pdata[, covariates, drop = FALSE])
    n_drop <- sum(!complete)
    if (n_drop > 0) {
      cat("Dropping", n_drop,
          "samples with missing covariate values:",
          paste(covariates, collapse = ", "), "\n")
      merged_pdata <- merged_pdata[complete, ]
      merged_exprs <- merged_exprs[, complete, drop = FALSE]
      imputed_mask <- imputed_mask[, complete, drop = FALSE]
      merged_samples <- colnames(merged_exprs)
    }
    cat("Using covariates:", paste(covariates, collapse = ", "), "\n")
  }

  # Create batch variable
  batch <- as.factor(merged_pdata$secondaryaccession)

  # Create biological group
  bio_group <- merged_pdata[[config$phenotype$group_column]]

  # Build ComBat mod formula (biological group + covariates)
  combat_bio_cov <- config$phenotype$combat_bio_covariate
  if (!is.null(combat_bio_cov)) {
    if (!combat_bio_cov %in% colnames(merged_pdata)) {
      merged_pdata[[combat_bio_cov]] <- compute_ga_weeks(merged_pdata)
    }
    ga_vals <- merged_pdata[[combat_bio_cov]]
    bio_degree <- config$phenotype$combat_bio_degree
    bio_spline_df <- config$phenotype$combat_bio_spline_df
    combat_mod_data <- data.frame(row.names = seq_along(ga_vals))
    if (!is.null(bio_spline_df)) {
      ns_mat <- splines::ns(ga_vals, df = bio_spline_df)
      for (j in seq_len(ncol(ns_mat)))
        combat_mod_data[[paste0("ga_ns", j)]] <- ns_mat[, j]
      bio_terms <- paste0("ga_ns", seq_len(ncol(ns_mat)))
      cat("  ComBat bio covariate: ns(", combat_bio_cov, ", df=", bio_spline_df, ")\n")
    } else if (!is.null(bio_degree) && bio_degree > 1) {
      poly_mat <- poly(ga_vals, degree = bio_degree)
      for (j in seq_len(ncol(poly_mat)))
        combat_mod_data[[paste0("ga_poly", j)]] <- poly_mat[, j]
      bio_terms <- paste0("ga_poly", seq_len(ncol(poly_mat)))
      cat("  ComBat bio covariate: poly(", combat_bio_cov, ",", bio_degree, ")\n")
    } else {
      combat_mod_data$ga_linear <- ga_vals
      bio_terms <- "ga_linear"
      cat("  ComBat bio covariate:", combat_bio_cov, "(linear)\n")
    }
    for (cov in covariates) combat_mod_data[[cov]] <- merged_pdata[[cov]]
    combat_mod_rhs <- paste(c(bio_terms, covariates), collapse = " + ")
  } else {
    combat_mod_data <- data.frame(bio_group = bio_group,
                                  stringsAsFactors = FALSE)
    for (cov in covariates) combat_mod_data[[cov]] <- merged_pdata[[cov]]
    combat_mod_rhs <- paste(c("bio_group", covariates), collapse = " + ")
  }
  combat_mod <- model.matrix(as.formula(paste("~", combat_mod_rhs)),
                             data = combat_mod_data)

  for (norm_method in normalization_methods) {
    cat("\n  Normalization:", norm_method, "\n")

    if (norm_method %in% c("harmonizr", "github_harmonizr") && imp_method != "none") {
      cat("  Skipping ", norm_method, ": only valid with imp_method='none'\n")
      next
    }

    result_key <- paste0(imp_method, "_", norm_method)

    # Apply normalization
    # "dwd" branch disabled: fails on single-sample batches.
    ruv_W <- NULL
    normalized <- switch(
      norm_method,
      "combat" = normalize_combat(merged_exprs, batch, mod = combat_mod),
      "combat_ref" = {
        ref_batch <- config[["normalization"]][["combat_ref"]][["ref_batch"]]
        if (is.null(ref_batch)) stop("combat_ref requires ref_batch in config")
        normalize_combat_ref(merged_exprs, batch, mod = combat_mod,
                             ref_batch = ref_batch)
      },
      # "dwd" = normalize_dwd(merged_exprs, batch),
      "mean_center" = normalize_mean_center(merged_exprs, batch),
      "batch_in_limma" = merged_exprs,
      "ruv" = {
        k <- config[["normalization"]][["ruv"]][["k"]] %||% 2
        ruv_result <- normalize_ruv(merged_exprs, ruv_control_genes, k = k)
        ruv_W <- ruv_result$W
        ruv_result$exprs
      },
      "ruvinv" = {
        merged_exprs
      },
      "bruv" = {
        k <- config[["normalization"]][["bruv"]][["k"]] %||% 2
        bruv_result <- normalize_bruv(
          merged_exprs, ruv_control_genes, k = k,
          sample_dataset = as.character(merged_pdata$secondaryaccession),
          sample_group = merged_pdata[[config$phenotype$group_column]]
        )
        ruv_W <- bruv_result$W
        bruv_result$exprs
      },
      "harmonizr" = {
        harmonizr_ref <- config[["normalization"]][["harmonizr"]][["ref_batch"]]
        inc_mat <- incomplete$matrix[, merged_samples, drop = FALSE]
        normalize_harmonizr(
          incomplete_matrix = inc_mat,
          batch = batch,
          mod_data = combat_mod_data,
          mod_formula = combat_mod_rhs,
          ref_batch = harmonizr_ref
        )
      },
      "github_harmonizr" = {
        gh_cfg <- config[["normalization"]][["github_harmonizr"]]
        inc_mat <- incomplete$matrix[, merged_samples, drop = FALSE]
        normalize_github_harmonizr(
          incomplete_matrix = inc_mat,
          batch = batch,
          algorithm    = gh_cfg[["algorithm"]] %||% "ComBat",
          ComBat_mode  = gh_cfg[["ComBat_mode"]] %||% 1,
          sort         = gh_cfg[["sort"]] %||% FALSE,
          block        = gh_cfg[["block"]]
        )
      },
      merged_exprs
    )

    # HarmonizR returns NAs (uncorrected cells) — set up weight mask
    if (norm_method %in% c("harmonizr", "github_harmonizr")) {
      imputed_mask <- is.na(normalized)
      dimnames(imputed_mask) <- dimnames(normalized)
      normalized[is.na(normalized)] <- 0
    }

    # Save normalized matrix
    if (config$output$save_normalized_matrices %||% FALSE) {
      norm_file <- file.path(output_dir, paste0("exprs_", result_key, ".tsv"))
      write.table(normalized, norm_file, sep = "\t", quote = FALSE)
      cat("  Saved normalized matrix:", norm_file, "\n")
    }

    # PCA plot of normalized data
    pc <- prcomp(t(normalized), center = TRUE, scale. = FALSE)
    var_pct <- 100 * summary(pc)$importance[2, ]
    pca_group <- merged_pdata[[config$phenotype$group_column]]
    pca_batch <- as.character(merged_pdata$secondaryaccession)

    group_lvls <- unique(pca_group)
    grp_pal <- c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F",
                 "#E5C494", "#B3B3B3")
    grp_col <- setNames(grp_pal[seq_along(group_lvls)], group_lvls)

    ds_lvls <- unique(pca_batch)
    ds_pch <- setNames(c(0:4, 6:8, 15:18)[seq_along(ds_lvls)], ds_lvls)

    pca_file <- file.path(output_dir, paste0("pca_", result_key, ".png"))
    png(pca_file, width = 1400, height = 900, res = 120)
    par(mar = c(5, 5, 3, 2))
    plot(pc$x[, 1], pc$x[, 2],
         col = grp_col[pca_group], pch = ds_pch[pca_batch], cex = 1.3,
         xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
         ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
         main = sprintf("PCA — %s (%d genes, %d samples)",
                        result_key, nrow(normalized), ncol(normalized)))
    grp_n <- table(factor(pca_group, levels = group_lvls))
    ds_n <- table(factor(pca_batch, levels = ds_lvls))
    legend("bottomright",
           legend = c(sprintf("%s (n=%d)", group_lvls, grp_n),
                      "", sprintf("%s (n=%d)", ds_lvls, ds_n)),
           col = c(grp_col[group_lvls], NA, rep("grey30", length(ds_lvls))),
           pch = c(rep(15, length(group_lvls)), NA, ds_pch[ds_lvls]),
           pt.cex = c(rep(1.5, length(group_lvls)), NA, rep(1.3, length(ds_lvls))),
           cex = 0.6, bty = "n")
    dev.off()
    cat("  PCA plot saved:", pca_file, "\n")

    # ============================================================
    # Step 5: Differential Expression Analysis
    # ============================================================

    cat("\n  Running differential expression...\n")

    # Filter to comparison groups
    keep_samples <- merged_pdata[[config$phenotype$group_column]] %in%
      c(config$phenotype$baseline, config$phenotype$contrast)

    de_exprs <- normalized[, keep_samples]
    de_pdata <- merged_pdata[keep_samples, ]

    # Align the imputed-cell mask to the DE matrix (ComBat may have
    # dropped zero-variance rows).
    de_mask <- imputed_mask[rownames(de_exprs), keep_samples, drop = FALSE]

    # Create design matrix (group + covariates, optionally + batch)
    group <- factor(de_pdata[[config$phenotype$group_column]],
                    levels = c(config$phenotype$baseline,
                               config$phenotype$contrast))

    if (norm_method == "ruvinv") {
      # RUVinv does its own DE testing — bypass limma
      cov_df <- NULL
      if (length(covariates) > 0) {
        cov_df <- de_pdata[, covariates, drop = FALSE]
      }
      lambda <- config[["normalization"]][["ruvinv"]][["lambda"]]
      ruvinv_result <- normalize_ruvinv(
        de_exprs, ruv_control_genes, group,
        covariates_df = cov_df, lambda = lambda
      )
      de_results <- ruvinv_result$de_results
    } else {
      design_data <- data.frame(group = group, stringsAsFactors = FALSE)
      for (cov in covariates) design_data[[cov]] <- de_pdata[[cov]]
      design_terms <- c("group", covariates)
      if (norm_method == "batch_in_limma") {
        de_batch <- droplevels(as.factor(de_pdata$secondaryaccession))
        if (nlevels(de_batch) > 1) {
          design_data$batch <- de_batch
          design_terms <- c("group", "batch", covariates)
          cat(sprintf("  Including batch in limma design (%d levels)\n", nlevels(de_batch)))
        } else {
          cat("  Warning: only 1 batch level after filtering, skipping batch term\n")
        }
      }
      if (norm_method %in% c("ruv", "bruv") && !is.null(ruv_W)) {
        de_W <- ruv_W[keep_samples, , drop = FALSE]
        for (j in seq_len(ncol(de_W))) {
          wname <- colnames(de_W)[j]
          design_data[[wname]] <- de_W[, j]
          design_terms <- c(design_terms, wname)
        }
        cat(sprintf("  Including %d RUV factors (W) in limma design\n", ncol(de_W)))
      }
      design_rhs <- paste(design_terms, collapse = " + ")
      design <- model.matrix(as.formula(paste("~", design_rhs)),
                             data = design_data)

      # Handle technical replicates via duplicateCorrelation (e.g. GSE55439).
      block_col <- de_pdata$technical_replicate_block
      if (!is.null(block_col) && any(nzchar(block_col))) {
        block_vec <- block_col
        block_vec[!nzchar(block_vec)] <- paste0("singleton_", seq_len(sum(!nzchar(block_vec))))
        block_vec <- as.factor(block_vec)
        cat(sprintf("  duplicateCorrelation: %d samples in %d blocks (%d singleton)\n",
                    length(block_vec), nlevels(block_vec), sum(grepl("^singleton_", block_vec))))
        corfit <- duplicateCorrelation(de_exprs, design, block = block_vec)
        cat(sprintf("  Consensus correlation: %.4f\n", corfit$consensus.correlation))
      } else {
        block_vec <- NULL
        corfit <- NULL
      }

      # Per-cell limma weights for imputed cells.
      imp_w <- config$de$imputed_cell_weight %||% 1.0
      if (imp_w != 1.0 && any(de_mask)) {
        W <- matrix(1.0, nrow(de_exprs), ncol(de_exprs),
                    dimnames = dimnames(de_exprs))
        W[de_mask] <- imp_w
        cat(sprintf("  Imputed-cell weight: %.3f  (downweighting %d cells)\n",
                    imp_w, sum(de_mask)))
        if (!is.null(corfit)) {
          fit <- lmFit(de_exprs, design, weights = W,
                       block = block_vec, correlation = corfit$consensus.correlation)
        } else {
          fit <- lmFit(de_exprs, design, weights = W)
        }
      } else {
        if (!is.null(corfit)) {
          fit <- lmFit(de_exprs, design,
                       block = block_vec, correlation = corfit$consensus.correlation)
        } else {
          fit <- lmFit(de_exprs, design)
        }
      }
      fit <- eBayes(fit)

      de_results <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
      de_results$gene <- rownames(de_results)
      de_results <- de_results[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    }

    # Apply thresholds
    fdr_thresh <- config$thresholds$fdr %||% 0.05
    logfc_thresh <- config$thresholds$logfc %||% 1.0

    sig_mask <- !is.na(de_results$adj.P.Val) & !is.na(de_results$logFC) &
      de_results$adj.P.Val < fdr_thresh & abs(de_results$logFC) > logfc_thresh
    de_significant <- de_results[sig_mask, ]

    cat("  Total genes:", nrow(de_results), "\n")
    cat("  Significant (FDR <", fdr_thresh, ", |logFC| >", logfc_thresh, "):",
        nrow(de_significant), "\n")
    cat("    Up-regulated:", sum(de_significant$logFC > 0), "\n")
    cat("    Down-regulated:", sum(de_significant$logFC < 0), "\n")

    # Save DE results
    de_file <- file.path(output_dir, paste0("difexp_", result_key, ".tsv"))
    write.table(de_results, de_file, sep = "\t", row.names = FALSE, quote = FALSE)

    sig_file <- file.path(output_dir, paste0("difexp_significant_", result_key, ".tsv"))
    write.table(de_significant, sig_file, sep = "\t", row.names = FALSE, quote = FALSE)

    # Store results for comparison
    all_results[[result_key]] <- list(
      imputation = imp_method,
      normalization = norm_method,
      n_genes = nrow(de_results),
      n_significant = nrow(de_significant),
      n_up = sum(de_significant$logFC > 0),
      n_down = sum(de_significant$logFC < 0),
      de_results = de_results,
      de_significant = de_significant
    )
  }
}

# ============================================================
# Step 6: Method Comparison
# ============================================================

cat("\n=== Step 6: Method Comparison ===\n\n")

# Summary table
comparison_df <- do.call(rbind, lapply(names(all_results), function(key) {
  r <- all_results[[key]]
  data.frame(
    method = key,
    imputation = r$imputation,
    normalization = r$normalization,
    n_genes = r$n_genes,
    n_significant = r$n_significant,
    n_up = r$n_up,
    n_down = r$n_down,
    stringsAsFactors = FALSE
  )
}))

print(comparison_df)

# Save comparison
write.csv(comparison_df, file.path(output_dir, "method_comparison.csv"),
          row.names = FALSE)

# Calculate overlaps between methods
if (length(all_results) > 1) {
  cat("\n=== DEG Overlap Between Methods ===\n\n")

  method_names <- names(all_results)
  overlap_matrix <- matrix(NA, nrow = length(method_names), ncol = length(method_names),
                           dimnames = list(method_names, method_names))

  for (i in seq_along(method_names)) {
    for (j in seq_along(method_names)) {
      genes_i <- all_results[[method_names[i]]]$de_significant$gene
      genes_j <- all_results[[method_names[j]]]$de_significant$gene
      overlap_matrix[i, j] <- length(intersect(genes_i, genes_j))
    }
  }

  cat("Overlap matrix (significant genes):\n")
  print(overlap_matrix)

  write.csv(overlap_matrix, file.path(output_dir, "deg_overlap_matrix.csv"))
}

# ============================================================
# Step 7: Visualization (NA staircases with FDR overlay)
# ============================================================

cat("\n=== Step 7: Visualization ===\n")

gene_fdr_list <- list()
for (key in names(all_results)) {
  gene_fdr_list[[key]] <- setNames(
    all_results[[key]]$de_results$adj.P.Val,
    as.character(all_results[[key]]$de_results$gene)
  )
}

staircase_sample_ds <- rep(names(exprs_list), sapply(exprs_list, ncol))
names(staircase_sample_ds) <- unlist(lapply(exprs_list, colnames))
staircase_title <- paste0(config$phenotype$baseline, " vs ", config$phenotype$contrast)
staircase_group <- phenodata[[config$phenotype$group_column]][
  match(names(staircase_sample_ds), phenodata$arraydatafile_exprscolumnnames)
]

# Build full pre-drop matrix so dropped genes appear as "1 group all-NA"
if (length(dropped_group_genes) > 0) {
  staircase_inc <- create_incomplete_matrix(
    exprs_list_pre_drop, min_datasets = config$coverage$threshold
  )
  staircase_matrix <- staircase_inc$matrix
} else {
  staircase_matrix <- incomplete$matrix
}

plot_na_staircase(staircase_matrix, staircase_sample_ds,
                  staircase_title, file.path(output_dir, "na_staircase.png"),
                  sample_group = staircase_group,
                  gene_fdr_list = gene_fdr_list,
                  baseline = config$phenotype$baseline,
                  contrast = config$phenotype$contrast,
                  coverage_threshold = selected_threshold / n_datasets_total)

# ============================================================
# Save Summary
# ============================================================

summary_file <- file.path(output_dir, "summary.txt")
summary_conn <- file(summary_file, "w")
writeLines(c(
  "Phase 2B: Direct Merge Improvements Summary",
  "==========================================",
  "",
  paste("Date:", Sys.time()),
  paste("Config:", config_path),
  "",
  paste("Comparison:", config$phenotype$contrast, "vs", config$phenotype$baseline),
  paste("Datasets:", paste(datasets, collapse = ", ")),
  if (!is.null(config$coverage$min_datasets)) {
    paste("Coverage mode: fixed min_datasets =", selected_threshold, "of", n_datasets_total)
  } else {
    paste("Coverage mode: max_imputation_allowed =", max_imp_allowed,
          "-> selected min_datasets =", selected_threshold, "of", n_datasets_total)
  },
  "",
  "Gene Recovery: skipped (compare_gene_recovery disabled)",
  "",
  "Method Comparison:",
  capture.output(print(comparison_df)),
  "",
  "Imputed-gene coverage per comparison group:",
  imputed_coverage$lines,
  "",
  if (!is.null(validation_results)) {
    c("Imputation Validation:",
      capture.output(print(aggregate(correlation ~ method, validation_results, mean))))
  } else {
    "Imputation validation: not run"
  }
), summary_conn)
close(summary_conn)

cat("\nSaved summary to:", summary_file, "\n")

# Close logging
close_logging(log_file)

cat("\n=== Imputed-gene coverage per comparison group (recap) ===\n")
for (line in imputed_coverage$lines) cat(line, "\n", sep = "")

cat("\n============================================================\n")
cat("  Phase 2B Complete!\n")
cat("============================================================\n\n")

cat("Output files in", output_dir, ":\n")
for (f in list.files(output_dir, pattern = "\\.(tsv|csv|txt)$")) {
  cat("  ", f, "\n")
}
cat("\n")

# Merging datasets that measure different gene sets, and the gene-coverage
# safety filter.

#' Create combined expression matrix with NAs for missing genes
#'
#' @param exprs_list Named list of expression matrices (genes x samples)
#' @param min_datasets Minimum number of datasets a gene must appear in (integer, 1..N)
#' @return List with: matrix (combined with NAs), gene_info, sample_info
#' @export
create_incomplete_matrix <- function(exprs_list, min_datasets = 1L) {

  n_datasets <- length(exprs_list)
  cat("Creating incomplete matrix from", n_datasets, "datasets\n")

  # Get all genes across datasets
  all_genes <- unique(unlist(lapply(exprs_list, rownames)))
  cat("  Total unique genes:", length(all_genes), "\n")

  # Calculate gene presence (vectorized for speed)
  # Create a presence matrix: genes x datasets
  gene_sets <- lapply(exprs_list, rownames)
  presence_matrix <- sapply(gene_sets, function(genes) all_genes %in% genes)
  rownames(presence_matrix) <- all_genes

  gene_presence <- rowSums(presence_matrix)
  gene_coverage <- gene_presence / n_datasets

  # Filter genes by minimum dataset count
  genes_to_keep <- names(gene_presence)[gene_presence >= min_datasets]
  cat("  Genes in >=", min_datasets, "of", n_datasets, "datasets:",
      length(genes_to_keep), "\n")

  # Get all samples
  all_samples <- unlist(lapply(exprs_list, colnames))
  cat("  Total samples:", length(all_samples), "\n")

  # Create sample-to-dataset mapping
  sample_dataset <- rep(names(exprs_list), sapply(exprs_list, ncol))
  names(sample_dataset) <- all_samples

  # Create incomplete matrix
  combined <- matrix(
    NA,
    nrow = length(genes_to_keep),
    ncol = length(all_samples),
    dimnames = list(genes_to_keep, all_samples)
  )

  # Fill in values
  for (ds_name in names(exprs_list)) {
    exprs <- exprs_list[[ds_name]]
    common_genes <- intersect(genes_to_keep, rownames(exprs))
    combined[common_genes, colnames(exprs)] <- as.matrix(exprs[common_genes, ])
  }

  # Calculate missingness
  n_missing <- sum(is.na(combined))
  pct_missing <- round(100 * n_missing / length(combined), 1)
  cat("  Missing values:", n_missing, "(", pct_missing, "%)\n")

  # Gene info
  gene_info <- data.frame(
    gene = genes_to_keep,
    n_datasets = gene_presence[genes_to_keep],
    coverage = gene_coverage[genes_to_keep],
    n_missing = rowSums(is.na(combined)),
    stringsAsFactors = FALSE
  )

  # Sample info
  sample_info <- data.frame(
    sample = all_samples,
    dataset = sample_dataset[all_samples],
    n_missing = colSums(is.na(combined)),
    stringsAsFactors = FALSE
  )

  list(
    matrix = combined,
    gene_info = gene_info,
    sample_info = sample_info,
    n_datasets = n_datasets,
    min_datasets = min_datasets
  )
}


# Inner join: only the genes every dataset measures, no missing values.
merge_common_genes <- function(exprs_list) {
  common <- Reduce(intersect, lapply(exprs_list, rownames))
  # unname() keeps the dataset name out of the column names
  do.call(cbind, unname(lapply(exprs_list, function(e) e[common, ])))
}


# For every gene that is missing from at least one dataset (and so will be
# imputed there): in how many datasets it is measured, and how many measured
# samples it has in each comparison group. A gene with no measured sample in
# one group would have that group's expression entirely synthesised.
imputed_gene_group_coverage <- function(exprs_list, phenodata, config) {
  baseline <- config$phenotype$baseline
  contrast <- config$phenotype$contrast
  gene_sets <- lapply(exprs_list, rownames)
  genes <- setdiff(Reduce(union, gene_sets), Reduce(intersect, gene_sets))
  if (length(genes) == 0)
    return(data.frame(gene = character(0), in_both_groups = logical(0)))
  measured <- sapply(gene_sets, function(s) genes %in% s)          # genes x datasets
  if (!is.matrix(measured)) measured <- matrix(measured, nrow = length(genes))
  group_of <- setNames(phenodata[[config$phenotype$group_column]], phenodata[[SAMPLE_ID]])
  per_dataset <- function(group) sapply(exprs_list, function(e)
    sum(group_of[colnames(e)] == group, na.rm = TRUE))
  n_base <- as.integer(measured %*% per_dataset(baseline))
  n_cont <- as.integer(measured %*% per_dataset(contrast))
  coverage <- data.frame(gene = genes,
                         n_datasets_with_gene = as.integer(rowSums(measured)),
                         n_base, n_cont,
                         in_both_groups = n_base >= 1 & n_cont >= 1,
                         stringsAsFactors = FALSE)
  colnames(coverage)[3:4] <- paste0("non_imputed_", make.names(c(baseline, contrast)))
  coverage
}

# Safety filter: remove the genes that have no measured sample in one of the
# two comparison groups.
drop_genes <- function(exprs_list, genes) {
  lapply(exprs_list, function(e) e[!(rownames(e) %in% genes), , drop = FALSE])
}

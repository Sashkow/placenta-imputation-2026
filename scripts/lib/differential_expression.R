# Differential expression with limma.

library(limma)

# Design matrix: comparison group (baseline first), optionally the dataset of
# origin, then covariates.
de_design <- function(pdata, config, covariates = config$phenotype$covariates %||% c(),
                      batch_in_design = FALSE) {
  ph <- config$phenotype
  design_data <- data.frame(group = factor(pdata[[ph$group_column]],
                                           levels = c(ph$baseline, ph$contrast)),
                            stringsAsFactors = FALSE)
  if (batch_in_design) {
    batch <- droplevels(as.factor(pdata[[DATASET_ID]]))
    if (nlevels(batch) > 1) design_data$batch <- batch
  }
  for (cov in covariates) design_data[[cov]] <- pdata[[cov]]
  model.matrix(as.formula(paste("~", paste(colnames(design_data), collapse = " + "))),
               data = design_data)
}

# Per-cell weights: `weight` for imputed cells, 1 for measured ones.
# NULL (an unweighted fit) when the weight is 1 or nothing was imputed.
imputed_cell_weights <- function(imputed_mask, weight) {
  if (weight == 1.0 || !any(imputed_mask)) return(NULL)
  W <- matrix(1.0, nrow(imputed_mask), ncol(imputed_mask), dimnames = dimnames(imputed_mask))
  W[imputed_mask] <- weight
  W
}

# limma fit of the group effect (second design column); one row per gene,
# in the row order of `exprs`.
fit_limma <- function(exprs, design, weights = NULL) {
  fit <- if (is.null(weights)) lmFit(exprs, design) else lmFit(exprs, design, weights = weights)
  table <- topTable(eBayes(fit), coef = 2, number = Inf, sort.by = "none")
  table$gene <- rownames(table)
  table[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}

# Differentially expressed: FDR below the threshold and |logFC| above it.
is_significant <- function(de, fdr = 0.05, logfc = 1.0) {
  !is.na(de$adj.P.Val) & !is.na(de$logFC) & de$adj.P.Val < fdr & abs(de$logFC) > logfc
}

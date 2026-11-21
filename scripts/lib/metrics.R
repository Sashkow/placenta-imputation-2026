# Agreement metrics between two differential-expression results.

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

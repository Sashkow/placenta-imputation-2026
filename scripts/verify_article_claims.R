#!/usr/bin/env Rscript
# Compares pipeline output against quantitative claims in the article.
# Run from the companion repo root: Rscript scripts/verify_article_claims.R
# Works with git-tracked data only (no pipeline re-run required).

library(tools)

base_dir <- getwd()
pipeline_dir <- file.path(base_dir, "data", "pipeline")
main_dir <- file.path(pipeline_dir, "main")
sens_dir <- file.path(pipeline_dir, "sensitivity")
balanced_dir <- file.path(pipeline_dir, "balanced")
validation_dir <- file.path(pipeline_dir, "validation")
enrichment_dir <- file.path(main_dir, "enrichment")

results <- data.frame(
  section = character(),
  claim = character(),
  article_value = character(),
  reproduced_value = character(),
  status = character(),
  stringsAsFactors = FALSE
)

add_result <- function(section, claim, article_val, reproduced_val, tol = NULL) {
  if (is.null(tol)) {
    status <- if (article_val == reproduced_val) "MATCH" else "MISMATCH"
  } else {
    status <- if (abs(as.numeric(article_val) - as.numeric(reproduced_val)) <= tol) "MATCH" else "MISMATCH"
  }
  results <<- rbind(results, data.frame(
    section = section,
    claim = claim,
    article_value = as.character(article_val),
    reproduced_value = as.character(reproduced_val),
    status = status,
    stringsAsFactors = FALSE
  ))
}

add_not_verified <- function(section, claim, reason) {
  results <<- rbind(results, data.frame(
    section = section,
    claim = claim,
    article_value = "—",
    reproduced_value = "—",
    status = paste0("NOT VERIFIED (", reason, ")"),
    stringsAsFactors = FALSE
  ))
}

# --- 1. Gene recovery (Table 1) ---

# Use exprs_none_combat_ref for intersection count (same 8260 genes, git-tracked)
none_combat_ref_file <- file.path(main_dir, "exprs_none_combat_ref.tsv")
soft_expr_file <- file.path(main_dir, "exprs_imputed_softimpute.tsv")
none_imputed_file <- file.path(main_dir, "exprs_imputed_none.tsv")

if (file.exists(none_imputed_file)) {
  none_expr <- read.delim(none_imputed_file, row.names = 1, check.names = FALSE)
  int_genes <- rownames(none_expr)
} else if (file.exists(none_combat_ref_file)) {
  none_expr <- read.delim(none_combat_ref_file, row.names = 1, check.names = FALSE)
  int_genes <- rownames(none_expr)
} else {
  int_genes <- NULL
}

if (!is.null(int_genes)) {
  add_result("Table 1", "Intersection gene count (k=6)", "8260", length(int_genes))
} else {
  add_not_verified("Table 1", "Intersection gene count (k=6)", "no expression matrix")
}

if (file.exists(soft_expr_file)) {
  soft_expr <- read.delim(soft_expr_file, row.names = 1, check.names = FALSE)
  add_result("Table 1", "Union gene count (k=1)", "17531", nrow(soft_expr))
} else {
  add_not_verified("Table 1", "Union gene count (k=1)", "exprs_imputed_softimpute.tsv missing")
}

summary_file <- file.path(main_dir, "summary.txt")
if (file.exists(summary_file)) {
  summary_txt <- readLines(summary_file)
  drop_line <- grep("Dropping.*imputed genes", summary_txt, value = TRUE)
  n_dropped <- as.integer(gsub(".*Dropping (\\d+).*", "\\1", drop_line))
  add_result("Table 1", "Safety-drop genes removed", "298", n_dropped)
} else {
  add_not_verified("Table 1", "Safety-drop genes removed", "summary.txt missing")
}

# --- 2. DEG counts (Table 4) ---

degs_soft_ref_file <- file.path(main_dir, "difexp_significant_softimpute_combat_ref.tsv")
degs_soft_ref <- read.delim(degs_soft_ref_file)
add_result("Table 4", "SoftImpute+ComBat-ref total DEGs", "538", nrow(degs_soft_ref))
add_result("Table 4", "SoftImpute+ComBat-ref up", "446",
           sum(degs_soft_ref$logFC > 0))
add_result("Table 4", "SoftImpute+ComBat-ref down", "92",
           sum(degs_soft_ref$logFC < 0))

degs_none_ref <- read.delim(file.path(main_dir,
                            "difexp_significant_none_combat_ref.tsv"))
add_result("Table 4", "Intersection-only (none+ref) total DEGs", "277",
           nrow(degs_none_ref))
add_result("Table 4", "Intersection-only up", "224",
           sum(degs_none_ref$logFC > 0))
add_result("Table 4", "Intersection-only down", "53",
           sum(degs_none_ref$logFC < 0))

degs_soft_file <- file.path(main_dir, "difexp_significant_softimpute_combat.tsv")
if (file.exists(degs_soft_file)) {
  degs_soft <- read.delim(degs_soft_file)
  add_result("Table 4", "SoftImpute+ComBat (no ref) total DEGs", "504",
             nrow(degs_soft))
  add_result("Table 4", "SoftImpute+ComBat up", "421",
             sum(degs_soft$logFC > 0))
  add_result("Table 4", "SoftImpute+ComBat down", "83",
             sum(degs_soft$logFC < 0))
} else {
  add_not_verified("Table 4", "SoftImpute+ComBat (no ref) DEGs",
                   "difexp_significant_softimpute_combat.tsv missing")
}

# --- 3. Imputation validation (Table 3) ---

iv_file <- file.path(main_dir, "imputation_validation.csv")
if (file.exists(iv_file)) {
  iv <- read.csv(iv_file)
  mean_r <- mean(iv$correlation)
  add_result("Table 3", "Block-mask softImpute r", "0.816",
             round(mean_r, 3), tol = 0.001)
  add_result("Table 3", "Block-mask RMSE", "1.52",
             round(mean(iv$rmse), 2), tol = 0.015)
  add_result("Table 3", "Block-mask MAE", "1.17",
             round(mean(iv$mae), 2), tol = 0.015)
} else {
  add_not_verified("Table 3", "Block-mask imputation metrics",
                   "imputation_validation.csv missing")
}

# --- 4. Subsampling validation (Section 3.4) ---

test1 <- read.delim(file.path(validation_dir,
                    "test1_first_trim_subsample.tsv"))
med_j_10 <- median(test1$jaccard_vs_full[test1$N_1st == 10])
add_result("Section 3.4", "Jaccard at N=10", "0.54",
           round(med_j_10, 2), tol = 0.01)

# --- 5. Split-half (Section 3.4) ---

test3 <- read.delim(file.path(validation_dir, "test3_split_half.tsv"))
med_j_ab <- median(test3$jaccard_a_vs_b)

add_result("Section 3.4", "Split-half median Jaccard between halves", "0.45",
           round(med_j_ab, 2), tol = 0.05)

# --- 6. Balanced reference (Section 3.4) ---

bal_ref_file <- file.path(balanced_dir,
                          "difexp_significant_softimpute_combat_ref.tsv")
if (file.exists(bal_ref_file)) {
  degs_balanced <- read.delim(bal_ref_file)
} else {
  bal_file <- file.path(balanced_dir, "difexp_softimpute_combat_ref.tsv")
  full_bal <- read.delim(bal_file)
  degs_balanced <- full_bal[abs(full_bal$logFC) >= 1 & full_bal$adj.P.Val < 0.05, ]
}

add_result("Section 3.4", "Balanced reference DEG count", "484",
           nrow(degs_balanced))

full_genes <- degs_soft_ref$gene
bal_genes <- degs_balanced$gene
jac_bal <- length(intersect(full_genes, bal_genes)) /
           length(union(full_genes, bal_genes))
ret_bal <- length(intersect(full_genes, bal_genes)) / length(full_genes)

add_result("Section 3.4", "Balanced vs full Jaccard", "0.841",
           round(jac_bal, 3), tol = 0.001)
add_result("Section 3.4", "Balanced retention of full DEGs", "86.8%",
           paste0(round(ret_bal * 100, 1), "%"))

# --- 7. ComBat sensitivity (Section 3.2) ---

if (file.exists(file.path(sens_dir, "exprs_softimpute_combat_ref.tsv")) &&
    !is.null(int_genes)) {
  main_full <- read.delim(file.path(main_dir,
                          "exprs_softimpute_combat_ref.tsv"),
                          row.names = 1, check.names = FALSE)
  sens_full <- read.delim(file.path(sens_dir,
                          "exprs_softimpute_combat_ref.tsv"),
                          row.names = 1, check.names = FALSE)
  shared_cols <- intersect(colnames(main_full), colnames(sens_full))

  m_vec <- as.vector(as.matrix(main_full[int_genes, shared_cols]))
  s_vec <- as.vector(as.matrix(sens_full[int_genes, shared_cols]))
  sens_r <- cor(m_vec, s_vec)
  sens_mad <- mean(abs(m_vec - s_vec))
  gene_mae <- rowMeans(abs(main_full[int_genes, shared_cols] -
                           sens_full[int_genes, shared_cols]))
  med_gene_mae <- median(gene_mae)

  add_result("Section 3.2", "Sensitivity Pearson r (8260 genes)", "0.9994",
             round(sens_r, 4), tol = 0.0005)
  add_result("Section 3.2", "Sensitivity mean abs diff", "0.031",
             round(sens_mad, 3), tol = 0.005)
  add_result("Section 3.2", "Sensitivity median per-gene MAE", "0.027",
             round(med_gene_mae, 3), tol = 0.005)
} else {
  add_not_verified("Section 3.2", "Sensitivity analysis", "data missing")
}

# --- 8. Enrichment (Section 3.5) ---

if (dir.exists(enrichment_dir)) {
  count_sig <- function(f) {
    d <- read.csv(f)
    sum(d$qvalue < 0.05, na.rm = TRUE)
  }

  f538_go <- file.path(enrichment_dir, "enrichment_full_538_GO_BP.csv")
  f538_kegg <- file.path(enrichment_dir, "enrichment_full_538_KEGG.csv")
  f277_go <- file.path(enrichment_dir, "enrichment_intersection_277_GO_BP.csv")
  f277_kegg <- file.path(enrichment_dir, "enrichment_intersection_277_KEGG.csv")
  f262_go <- file.path(enrichment_dir, "enrichment_gained_262_GO_BP.csv")
  f262_kegg <- file.path(enrichment_dir, "enrichment_gained_262_KEGG.csv")

  if (file.exists(f538_go)) {
    add_result("Enrichment", "538-DEG GO BP terms (q<0.05)", "678",
               count_sig(f538_go))
    add_result("Enrichment", "538-DEG KEGG pathways (q<0.05)", "43",
               count_sig(f538_kegg))
    add_result("Enrichment", "277-DEG GO BP terms", "292",
               count_sig(f277_go))
    add_result("Enrichment", "277-DEG KEGG pathways", "25",
               count_sig(f277_kegg))
    add_result("Enrichment", "262-DEG GO BP terms", "240",
               count_sig(f262_go))
    add_result("Enrichment", "262-DEG KEGG pathways", "24",
               count_sig(f262_kegg))
  } else {
    add_not_verified("Enrichment", "All enrichment counts",
                     "enrichment CSVs missing")
  }
} else {
  add_not_verified("Enrichment", "All enrichment counts",
                   "enrichment dir missing")
}

# --- Build report ---

n_match <- sum(results$status == "MATCH")
n_mismatch <- sum(results$status == "MISMATCH")
n_not_verified <- sum(grepl("NOT VERIFIED", results$status))
n_total <- nrow(results)

report <- c(
  "# Replication Report: Companion Repository vs Article Claims",
  "",
  paste0("**Generated:** ", Sys.Date()),
  paste0("**R version:** ", R.version.string),
  paste0("**Platform:** ", R.version$platform),
  "",
  "## Environment",
  ""
)

if (file.exists(file.path(base_dir, "renv.lock"))) {
  lock <- jsonlite::fromJSON(file.path(base_dir, "renv.lock"))
  report <- c(report, paste0("**renv R version:** ", lock$R$Version))
} else {
  report <- c(report, "**renv:** not present")
}

report <- c(report,
  "",
  "## Summary",
  "",
  paste0("- **Total claims verified:** ", n_total),
  paste0("- **Matches:** ", n_match),
  paste0("- **Mismatches:** ", n_mismatch),
  paste0("- **Not verified:** ", n_not_verified),
  "",
  "## Comparison Table",
  "",
  "| Section | Claim | Article | Reproduced | Status |",
  "|---------|-------|---------|------------|--------|"
)

for (i in seq_len(nrow(results))) {
  report <- c(report, sprintf("| %s | %s | %s | %s | %s |",
    results$section[i], results$claim[i], results$article_value[i],
    results$reproduced_value[i], results$status[i]))
}

report <- c(report, "", "## Reproduction Outcome", "")
if (n_mismatch == 0 && n_not_verified == 0) {
  report <- c(report,
    "All verified quantitative claims in the article are reproduced by",
    "the companion repository pipeline.")
} else if (n_mismatch == 0) {
  report <- c(report,
    paste0("All ", n_match, " verified claims match the article. ",
           n_not_verified, " claim(s) could not be verified due to missing data."))
} else {
  report <- c(report,
    paste0(n_mismatch, " claim(s) show discrepancies:"),
    "")
  mismatches <- results[results$status == "MISMATCH", ]
  for (i in seq_len(nrow(mismatches))) {
    report <- c(report, sprintf("- **%s — %s**: article says %s, reproduced %s",
      mismatches$section[i], mismatches$claim[i],
      mismatches$article_value[i], mismatches$reproduced_value[i]))
  }
}

report_path <- file.path(base_dir, "replication_report.md")
writeLines(report, report_path)
cat("Report written to:", report_path, "\n")
cat("\nQuick summary:", n_match, "match,", n_mismatch, "mismatch,",
    n_not_verified, "not verified out of", n_total, "claims\n")

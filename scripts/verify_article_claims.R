#!/usr/bin/env Rscript
# Compares pipeline output against quantitative claims in the article.
# Run from the companion repo root: Rscript scripts/verify_article_claims.R
# Works with git-tracked data only (no pipeline re-run required).

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
    status <- if (as.character(article_val) == as.character(reproduced_val)) "MATCH" else "MISMATCH"
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
  n_union_genes <- length(readLines(soft_expr_file)) - 1L
  add_result("Table 1", "Union gene count (k=1)", "17531", n_union_genes)
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

check_deg_counts <- function(section, label, file, expected_total, expected_up,
                             expected_down, optional = FALSE) {
  if (optional && !file.exists(file)) {
    add_not_verified(section, paste0(label, " DEGs"), paste0(basename(file), " missing"))
    return(NULL)
  }
  degs <- read.delim(file)
  add_result(section, paste0(label, " total DEGs"), expected_total, nrow(degs))
  add_result(section, paste0(label, " up"), expected_up, sum(degs$logFC > 0))
  add_result(section, paste0(label, " down"), expected_down, sum(degs$logFC < 0))
  degs
}

degs_soft_ref <- check_deg_counts(
  "Table 4", "SoftImpute+ComBat-ref",
  file.path(main_dir, "difexp_significant_softimpute_combat_ref.tsv"),
  "538", "446", "92")

check_deg_counts(
  "Table 4", "Intersection-only (none+ref)",
  file.path(main_dir, "difexp_significant_none_combat_ref.tsv"),
  "277", "224", "53")

check_deg_counts(
  "Table 4", "SoftImpute+ComBat (no ref)",
  file.path(main_dir, "difexp_significant_softimpute_combat.tsv"),
  "504", "421", "83", optional = TRUE)

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
shared_genes <- intersect(full_genes, bal_genes)
jac_bal <- length(shared_genes) / length(union(full_genes, bal_genes))
ret_bal <- length(shared_genes) / length(full_genes)

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

  m_sub <- main_full[int_genes, shared_cols]
  s_sub <- sens_full[int_genes, shared_cols]
  m_vec <- as.vector(as.matrix(m_sub))
  s_vec <- as.vector(as.matrix(s_sub))
  sens_r <- cor(m_vec, s_vec)
  gene_mae <- rowMeans(abs(m_sub - s_sub))
  sens_mad <- mean(gene_mae)
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

  enrichment_checks <- data.frame(
    label = c("538-DEG GO BP terms (q<0.05)", "538-DEG KEGG pathways (q<0.05)",
              "277-DEG GO BP terms", "277-DEG KEGG pathways",
              "262-DEG GO BP terms", "262-DEG KEGG pathways"),
    file = c("enrichment_full_538_GO_BP.csv", "enrichment_full_538_KEGG.csv",
             "enrichment_intersection_277_GO_BP.csv", "enrichment_intersection_277_KEGG.csv",
             "enrichment_gained_262_GO_BP.csv", "enrichment_gained_262_KEGG.csv"),
    expected = c("678", "43", "292", "25", "240", "24"),
    stringsAsFactors = FALSE
  )

  first_file <- file.path(enrichment_dir, enrichment_checks$file[1])
  if (file.exists(first_file)) {
    for (i in seq_len(nrow(enrichment_checks))) {
      add_result("Enrichment", enrichment_checks$label[i],
                 enrichment_checks$expected[i],
                 count_sig(file.path(enrichment_dir, enrichment_checks$file[i])))
    }
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

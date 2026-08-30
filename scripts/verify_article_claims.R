#!/usr/bin/env Rscript
# Compares pipeline output against quantitative claims in the article.
# Run from the companion repo root: Rscript scripts/verify_article_claims.R
# Works with git-tracked data only (no pipeline re-run required).

base_dir <- getwd()
pipeline_dir <- file.path(base_dir, "data", "pipeline")
main_dir <- file.path(pipeline_dir, "main")
balanced_dir <- file.path(pipeline_dir, "balanced")
validation_dir <- file.path(pipeline_dir, "validation")
enrichment_dir <- file.path(main_dir, "enrichment")

results <- data.frame(
  section = character(),
  claim = character(),
  article_value = character(),
  reproduced_value = character(),
  tolerance = character(),
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
    tolerance = if (is.null(tol)) "exact" else format(tol, scientific = FALSE),
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
    tolerance = "—",
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
  "530", "436", "94")

check_deg_counts(
  "Table 4", "Intersection-only (none+ref)",
  file.path(main_dir, "difexp_significant_none_combat_ref.tsv"),
  "277", "224", "53")

check_deg_counts(
  "Table 4", "SoftImpute+ComBat (no ref)",
  file.path(main_dir, "difexp_significant_softimpute_combat.tsv"),
  "525", "433", "92", optional = TRUE)

# --- 3. Imputation accuracy (Supp Table S1, unified holdout) ---

hs1 <- "data/pipeline/holdout/main/table_s1.csv"
if (file.exists(hs1)) {
  t1 <- read.csv(hs1, stringsAsFactors = FALSE)
  si <- function(sc, col) t1[[col]][t1$scheme == sc & t1$method == "softimpute"]
  add_result("Supp S1", "softImpute random-cell r", "0.993", round(si("random_cells", "r_mean"), 3), tol = 0.0005)
  add_result("Supp S1", "softImpute gene-dataset block r", "0.817", round(si("gene_dataset_block", "r_mean"), 3), tol = 0.0005)
  add_result("Supp S1", "softImpute block RMSE", "1.52", round(si("gene_dataset_block", "rmse_mean"), 2), tol = 0.005)
  add_result("Supp S1", "softImpute block MAE", "1.16", round(si("gene_dataset_block", "mae_mean"), 2), tol = 0.005)
  add_result("Supp S1", "softImpute worst-case block r", "0.824", round(si("progressive_tax_block", "r_mean"), 3), tol = 0.0005)
  add_result("Supp S1", "KNN block r", "0.683", round(t1$r_mean[t1$scheme == "gene_dataset_block" & t1$method == "knn"], 3), tol = 0.0005)
  add_result("Supp S1", "imputePCA worst-case r", "0.436", round(t1$r_mean[t1$scheme == "progressive_tax_block" & t1$method == "missmda"], 3), tol = 0.0005)
} else {
  add_not_verified("Supp S1", "Imputation accuracy (unified holdout)", "data/pipeline/holdout/main/table_s1.csv missing")
}

# --- 4. Subsampling validation (Section 3.4) ---

test1 <- read.delim(file.path(validation_dir,
                    "test1_first_trim_subsample.tsv"))
med_j_10 <- median(test1$jaccard_vs_full[test1$N_1st == 10])
add_result("Section 3.4", "Jaccard at N=10", "0.52",
           round(med_j_10, 2), tol = 0.01)

# --- 5. Split-half (Section 3.4) ---

test3 <- read.delim(file.path(validation_dir, "test3_split_half.tsv"))
med_j_ab <- median(test3$jaccard_a_vs_b)

# Tolerance tightened from 0.05 to 0.01: the loose value masked the fact that
# the article prose quoted an older validation run (0.45) while the shipped
# data gives 0.402.
add_result("Section 3.4", "Split-half median Jaccard between halves", "0.445",
           round(med_j_ab, 3), tol = 0.01)
add_result("Supp S3", "Split-half DEG retention", "86.1",
           round(100 * median(c(test3$overlap_a_vs_full, test3$overlap_b_vs_full)), 1),
           tol = 0.1)
add_result("Supp S3", "Split-half logFC CCC between halves", "0.763",
           round(median(test3$logfc_ccc_a_vs_b), 3), tol = 0.01)

# --- 4b. Convergence numbers synced to the canonical run ---
add_result("Supp S3", "Subsampling DEG retention at N=10", "90.6",
           round(100 * median(test1$overlap_vs_full[test1$N_1st == 10]), 1), tol = 0.1)
add_result("Supp S3", "Subsampling DEG retention at N=30", "92.7",
           round(100 * median(test1$overlap_vs_full[test1$N_1st == 30]), 1), tol = 0.1)
add_result("Supp S3", "Subsampling logFC CCC at N=10", "0.89",
           round(median(test1$logfc_ccc_vs_full[test1$N_1st == 10]), 2), tol = 0.01)
add_result("Supp S3", "Subsampling logFC CCC at N=100", "1.00",
           round(median(test1$logfc_ccc_vs_full[test1$N_1st == 100]), 2), tol = 0.01)
add_result("Section 3.4", "Jaccard at N=100", "0.98",
           round(median(test1$jaccard_vs_full[test1$N_1st == 100]), 2), tol = 0.01)

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

add_result("Section 3.4", "Balanced reference DEG count", "474",
           nrow(degs_balanced))

full_genes <- degs_soft_ref$gene
bal_genes <- degs_balanced$gene
shared_genes <- intersect(full_genes, bal_genes)
jac_bal <- length(shared_genes) / length(union(full_genes, bal_genes))
ret_bal <- length(shared_genes) / length(full_genes)

add_result("Section 3.4", "Balanced vs full Jaccard", "0.839",
           round(jac_bal, 3), tol = 0.001)
add_result("Section 3.4", "Balanced retention of full DEGs", "86.4%",
           paste0(round(ret_bal * 100, 1), "%"))

# --- 7. Batch-correction distortion check (Results 2.3) ---
# Article claim: post-ComBat-ref values for the 8,260 shared genes are
# practically identical between the intersection-only and softImpute runs.
# Both matrices come from the main pipeline run.

inter_file <- file.path(main_dir, "exprs_none_combat_ref.tsv")
soft_file <- file.path(main_dir, "exprs_softimpute_combat_ref.tsv")

if (file.exists(inter_file) && file.exists(soft_file)) {
  inter_post <- read.delim(inter_file, row.names = 1, check.names = FALSE)
  soft_post <- read.delim(soft_file, row.names = 1, check.names = FALSE)
  shared_genes <- intersect(rownames(inter_post), rownames(soft_post))
  i_sub <- as.matrix(inter_post[shared_genes, ])
  s_sub <- as.matrix(soft_post[shared_genes, colnames(inter_post)])

  add_result("Results 2.3", "Shared genes compared", "8260",
             length(shared_genes))
  add_result("Results 2.3", "Post-ComBat Pearson r (all cells)", "0.9999",
             round(cor(as.vector(i_sub), as.vector(s_sub)), 4), tol = 0.0005)
  gene_mae <- rowMeans(abs(i_sub - s_sub))
  add_result("Results 2.3", "Mean abs diff", "0.01",
             round(mean(abs(i_sub - s_sub)), 3), tol = 0.005)
  add_result("Results 2.3", "Median per-gene MAE", "0.009",
             round(median(gene_mae), 3), tol = 0.005)
} else {
  add_not_verified("Results 2.3", "Batch-correction distortion check",
                   "post-ComBat matrices missing")
}

# --- 8. Enrichment (Section 3.5) ---

if (dir.exists(enrichment_dir)) {
  count_sig <- function(f) {
    d <- read.csv(f)
    sum(d$qvalue < 0.05, na.rm = TRUE)
  }

  enrichment_checks <- data.frame(
    label = c("530-DEG GO BP terms (q<0.05)", "530-DEG KEGG pathways (q<0.05)",
              "277-DEG GO BP terms", "277-DEG KEGG pathways",
              "253-DEG GO BP terms", "253-DEG KEGG pathways"),
    file = c("enrichment_full_530_GO_BP.csv", "enrichment_full_530_KEGG.csv",
             "enrichment_intersection_277_GO_BP.csv", "enrichment_intersection_277_KEGG.csv",
             "enrichment_gained_253_GO_BP.csv", "enrichment_gained_253_KEGG.csv"),
    expected = c("659", "39", "292", "25", "250", "17"),
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

# --- 9. External concordance with Lykhenko 2021 (Discussion 3.2) ---
# "Gained" DEGs = significant in softImpute+ComBat-ref but not in the
# intersection-only run. Their direction and FDR in the prior study come
# from the FULL protein-coding limma table (most are below the 2021
# significance threshold, so the filtered DEG list cannot back this claim).

lykhenko_full_file <- file.path(base_dir, "data", "references",
                                "lykhenko_2021_full_protein_coding.csv")
none_sig_file <- file.path(main_dir, "difexp_significant_none_combat_ref.tsv")

if (file.exists(lykhenko_full_file) && !is.null(degs_soft_ref) &&
    file.exists(none_sig_file)) {
  lykhenko_full <- read.csv(lykhenko_full_file)
  degs_int <- read.delim(none_sig_file)
  gained <- degs_soft_ref[!(degs_soft_ref$gene %in% degs_int$gene), ]
  add_result("Discussion 3.2", "Gained DEGs (530 set minus 277 set)",
             "253", nrow(gained))
  matched <- merge(gained, lykhenko_full, by.x = "gene", by.y = "ENTREZID")
  add_result("Discussion 3.2", "Gained genes present in Lykhenko 2021 limma table",
             "233", nrow(matched))
  same_dir <- sign(matched$logFC.x) == sign(matched$logFC.y)
  add_result("Discussion 3.2", "Same direction of change (%)",
             "89.3", round(100 * mean(same_dir), 1), tol = 0.05)
  add_result("Discussion 3.2", "logFC Pearson r vs Lykhenko 2021",
             "0.654", round(cor(matched$logFC.x, matched$logFC.y), 3), tol = 0.0051)
  add_result("Discussion 3.2", "Already FDR-significant in Lykhenko 2021 (%)",
             "62.7", round(100 * mean(matched$adj.P.Val.y < 0.05), 1), tol = 0.05)
  ## Supp S4: composition of the 27 direction-discordant genes
  opp <- matched[!same_dir, ]
  add_result("Supp S4", "Direction-discordant gained genes",
             "25", nrow(opp))
  add_result("Supp S4", "Discordant genes FDR-significant in 2021",
             "0", sum(opp$adj.P.Val.y < 0.05))
  add_result("Supp S4", "Median 2021 |logFC| of discordant genes",
             "0.13", round(median(abs(opp$logFC.y)), 2), tol = 0.005)
  add_result("Supp S4", "Median 2021 |logFC| of concordant genes",
             "0.60", round(median(abs(matched$logFC.y[same_dir])), 2), tol = 0.005)
  sig21 <- matched[matched$adj.P.Val.y < 0.05, ]
  add_result("Supp S4", "2021-significant matched genes agreeing in direction",
             "146/146", paste0(sum(sign(sig21$logFC.x) == sign(sig21$logFC.y)),
                               "/", nrow(sig21)))
} else {
  add_not_verified("Discussion 3.2", "External concordance with Lykhenko 2021",
                   "lykhenko_2021_full_protein_coding.csv or DEG tables missing")
}

# --- 9. Enrichment q-values and term counts (Section 3.3) ---
# Previously unchecked; the article's q-values had drifted from the shipped CSVs.

enr <- function(f) {
  p <- file.path(enrichment_dir, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}
q_of <- function(d, desc) {
  if (is.null(d)) return(NA_real_)
  i <- which(d$Description == desc)
  if (!length(i)) NA_real_ else d$qvalue[i[1]]
}
fmt_q <- function(x) if (is.na(x)) "NA" else signif(x, 2)

go_full  <- enr("enrichment_full_530_GO_BP.csv")
kegg_full <- enr("enrichment_full_530_KEGG.csv")
go_int   <- enr("enrichment_intersection_277_GO_BP.csv")
kegg_int <- enr("enrichment_intersection_277_KEGG.csv")
go_gain  <- enr("enrichment_gained_253_GO_BP.csv")
kegg_gain <- enr("enrichment_gained_253_KEGG.csv")

if (!is.null(go_full) && !is.null(kegg_full)) {
  add_result("Section 3.3", "GO BP terms (full 530)", "659", nrow(go_full))
  add_result("Section 3.3", "KEGG pathways (full 530)", "39", nrow(kegg_full))
  add_result("Section 3.3", "GO BP terms (intersection 277)", "292", nrow(go_int))
  add_result("Section 3.3", "KEGG pathways (intersection 277)", "25", nrow(kegg_int))
  add_result("Section 3.3", "GO BP terms (gained 253)", "250", nrow(go_gain))
  add_result("Section 3.3", "KEGG pathways (gained 253)", "17", nrow(kegg_gain))

  add_result("Section 3.3", "q: positive regulation of cytokine production (full)",
             "7e-14", fmt_q(q_of(go_full, "positive regulation of cytokine production")),
             tol = 1e-15)
  add_result("Section 3.3", "q: leukocyte migration (full)",
             "6.7e-12", fmt_q(q_of(go_full, "leukocyte migration")), tol = 1e-12)
  add_result("Section 3.3", "q: chemotaxis (full)",
             "2.8e-11", fmt_q(q_of(go_full, "chemotaxis")), tol = 1e-12)
  add_result("Section 3.3", "q: humoral immune response (full)",
             "5.3e-11", fmt_q(q_of(go_full, "humoral immune response")), tol = 1e-11)
  add_result("Section 3.3", "q: hsa05150 S. aureus infection (full)",
             "2.2e-10", fmt_q(q_of(kegg_full, "Staphylococcus aureus infection")), tol = 1e-11)
  add_result("Section 3.3", "q: hsa04610 complement and coagulation (full)",
             "1e-09", fmt_q(q_of(kegg_full, "Complement and coagulation cascades")), tol = 1e-10)
  add_result("Section 3.3", "q: hsa04514 cell adhesion molecules (full)",
             "7.7e-06", fmt_q(q_of(kegg_full, "Cell adhesion molecule (CAM) interaction")), tol = 1e-7)
  add_result("Section 3.3", "q: chemotaxis (gained)",
             "2.4e-07", fmt_q(q_of(go_gain, "chemotaxis")), tol = 1e-7)
  add_result("Section 3.3", "q: leukocyte migration (gained)",
             "2.4e-07", fmt_q(q_of(go_gain, "leukocyte migration")), tol = 1e-7)
  add_result("Section 3.3", "q: hsa04145 phagocytosis (gained)",
             "0.0023", fmt_q(q_of(kegg_gain, "Phagocytosis")), tol = 1e-4)

  shared_go <- length(intersect(go_gain$Description, go_int$Description))
  shared_kegg <- length(intersect(kegg_gain$Description, kegg_int$Description))
  add_result("Section 3.4", "Gained GO BP terms shared with intersection", "83", shared_go)
  add_result("Section 3.4", "Gained KEGG pathways shared with intersection", "7", shared_kegg)
} else {
  add_not_verified("Section 3.3", "Enrichment q-values",
                   "data/pipeline/main/enrichment/ CSVs missing")
}

# --- 10. Gained-DEG composition (Section 3.3) ---
if (!is.null(degs_soft_ref) && file.exists(none_sig_file) &&
    file.exists(file.path(main_dir, "exprs_imputed_none.tsv"))) {
  int_only <- rownames(read.delim(file.path(main_dir, "exprs_imputed_none.tsv"),
                                  row.names = 1, check.names = FALSE))
  degs_int2 <- read.delim(none_sig_file)
  gained2 <- setdiff(degs_soft_ref$gene, degs_int2$gene)
  add_result("Section 3.3", "Gained DEGs that are intersection genes", "10",
             sum(gained2 %in% int_only))
  add_result("Section 3.3", "Gained DEGs testable only via imputation", "243",
             sum(!(gained2 %in% int_only)))
}

# --- 11. qPCR benchmark recovery (Supp S6) ---
qpcr_file <- file.path(base_dir, "data", "references", "qpcr_benchmark_genes.csv")
full_de_file <- file.path(main_dir, "difexp_softimpute_combat_ref.tsv")
if (file.exists(qpcr_file) && file.exists(full_de_file)) {
  bench <- read.csv(qpcr_file, stringsAsFactors = FALSE)
  full_de <- read.delim(full_de_file, stringsAsFactors = FALSE)
  idx <- match(as.character(bench$entrez_id), as.character(full_de$gene))
  hit <- !is.na(idx) & full_de$adj.P.Val[idx] < 0.05 &
         abs(full_de$logFC[idx]) > 1
  add_result("Supp S6", "qPCR benchmark genes recovered as DEGs", "11", sum(hit, na.rm = TRUE))
} else {
  add_not_verified("Supp S6", "qPCR benchmark recovery",
                   "qpcr_benchmark_genes.csv or full DE table missing")
}

# --- 12. GA-matched Prater re-analysis (Supp S5) ---
ga_file <- file.path(pipeline_dir, "ga_matched_prater", "prater_concordance.csv")
if (file.exists(ga_file)) {
  ga <- read.csv(ga_file)
  add_result("Supp S5", "GA-matched testable genes", "11761", ga$genes_tested[1])
  add_result("Supp S5", "GA-matched DEGs", "602", ga$degs[1])
  add_result("Supp S5", "GA-matched shared genes with Prater", "1992", ga$shared_genes[1])
  add_result("Supp S5", "GA-matched Prater Pearson r", "0.808",
             round(ga$prater_r[1], 3), tol = 0.002)
  add_result("Supp S5", "GA-matched Prater CCC", "0.755",
             round(ga$prater_ccc[1], 3), tol = 0.002)
} else {
  add_not_verified("Supp S5", "GA-matched Prater re-analysis",
                   "run scripts/ga_matched_prater.R first")
}

# --- 13. ComBat covariate sweep (Supp Tables S2/S3) ---
sweep_file <- file.path(pipeline_dir, "covariate_sweep", "table_s2_covariate_comparison.csv")
if (file.exists(sweep_file)) {
  sw <- read.csv(sweep_file, stringsAsFactors = FALSE)
  want <- c(categorical = 436, linear = 530, poly2 = 486, ns3 = 513)
  for (v in names(want)) {
    r <- sw[sw$variant == v, ]
    if (nrow(r)) add_result("Supp S6", paste0("Covariate sweep DEGs (", v, ")"),
                            want[[v]], r$degs[1])
  }
  wq <- c(categorical = 9, linear = 11, poly2 = 8, ns3 = 10)
  for (v in names(wq)) {
    r <- sw[sw$variant == v, ]
    if (nrow(r)) add_result("Supp S6", paste0("Covariate sweep qPCR recovered (", v, ")"),
                            wq[[v]], r$qpcr_recovered[1])
  }
} else {
  add_not_verified("Supp S6", "ComBat covariate sweep (Tables S2/S3)",
                   "run scripts/covariate_sweep.R first")
}

# --- 13b. Subsampling stability of the four covariate specs (Supp Table S4) ---
jac_file <- file.path(pipeline_dir, "covariate_sweep", "table_s3_covariate_stability.csv")
if (file.exists(jac_file)) {
  jc <- read.csv(jac_file, stringsAsFactors = FALSE)
  published <- list(
    "10"  = c(categorical = 0.453, linear = 0.481, poly2 = 0.429, ns3 = 0.443),
    "30"  = c(categorical = 0.625, linear = 0.647, poly2 = 0.565, ns3 = 0.573),
    "100" = c(categorical = 0.967, linear = 0.972, poly2 = 0.961, ns3 = 0.957))
  for (n in names(published)) {
    row <- jc[jc$N_1st == as.integer(n), ]
    if (!nrow(row)) next
    for (v in names(published[[n]])) {
      if (!v %in% colnames(row)) next
      add_result("Supp S6", sprintf("Within-size Jaccard N=%s (%s)", n, v),
                 published[[n]][[v]], round(row[[v]][1], 3), tol = 0.001)
    }
  }
  # the property the covariate choice actually rests on
  all_best <- all(jc$linear > jc$categorical & jc$linear > jc$poly2 & jc$linear > jc$ns3)
  add_result("Supp S6", "Linear GA most reproducible at every subsample size",
             "TRUE", as.character(all_best))
} else {
  add_not_verified("Supp S6", "Covariate subsampling stability (Table S4)",
                   "run the four config_validation_covariate_*.yaml test1 runs first")
}

# --- 14. Coverage patterns / reference-batch finding (test 4.1) ---
cov_file <- file.path(pipeline_dir, "test4", "coverage_patterns_gained.csv")
if (file.exists(cov_file)) {
  cp <- read.csv(cov_file, stringsAsFactors = FALSE)
  add_result("Section 3.3", "Gained DEGs carrying the ComBat reference batch",
             "245", sum(cp$n_genes[cp$has_ref_batch]))
  add_result("Section 3.3", "Gained DEGs lacking the ComBat reference batch",
             "8", sum(cp$n_genes[!cp$has_ref_batch]))
} else {
  add_not_verified("Section 3.3", "Coverage patterns of gained DEGs",
                   "run scripts/test4_coverage_patterns.R first")
}

# --- 15. Unified holdout validation (main-text holdout table; Supp S3-S6) ---

hd <- "data/pipeline/holdout"
cf <- file.path(hd, "main", "confusion_full_vs_masked.csv"); bdf <- file.path(hd, "main", "by_dataset.csv")
if (file.exists(cf) && file.exists(bdf)) {
  cfd <- read.csv(cf, stringsAsFactors = FALSE); bd <- read.csv(bdf, stringsAsFactors = FALSE)
  cv <- function(sc, st, col) cfd[[col]][cfd$scheme == sc & cfd$stratum == st]
  bv <- function(v, h, col) bd[[col]][bd$variant == v & bd$hidden == h]
  add_result("Holdout", "Block masking: masked genes per run", "5157", round(cv("gene_dataset_block", "overall", "genes")), tol = 0.5)
  add_result("Holdout", "Block masking: true DEGs lost (%)", "24.8", round(cv("gene_dataset_block", "overall", "lost_pct"), 1), tol = 0.05)
  add_result("Holdout", "Block masking: non-DEGs gained (%)", "0.21", round(cv("gene_dataset_block", "overall", "gained_pct"), 2), tol = 0.005)
  add_result("Holdout", "Block masking: genes gaining DEG status per run", "10", round(cv("gene_dataset_block", "overall", "not_to_deg")), tol = 0.5)
  add_result("Holdout", "Block masking, GSE100051 visible: lost (%)", "6.4", round(cv("gene_dataset_block", "ref_visible", "lost_pct"), 1), tol = 0.05)
  add_result("Holdout", "Block masking, GSE100051 hidden: lost (%)", "87.1", round(cv("gene_dataset_block", "ref_hidden", "lost_pct"), 1), tol = 0.05)
  add_result("Holdout", "Block masking, GSE100051 hidden: true DEGs lost per run", "36", round(cv("gene_dataset_block", "ref_hidden", "deg_to_not")), tol = 0.5)
  add_result("Holdout", "Block masking, GSE100051 hidden: |logFC| ratio", "0.48", round(cv("gene_dataset_block", "ref_hidden", "ratio"), 2), tol = 0.005)
  add_result("Holdout", "Random-cell masking: lost (%)", "4.2", round(cv("random_cells", "overall", "lost_pct"), 1), tol = 0.05)
  add_result("Holdout", "Random-cell masking: gained (%)", "0.16", round(cv("random_cells", "overall", "gained_pct"), 2), tol = 0.005)
  add_result("Holdout", "Worst-case masking: lost (%)", "63.9", round(cv("progressive_tax_block", "overall", "lost_pct"), 1), tol = 0.05)
  add_result("Holdout", "Worst-case masking, GSE100051 hidden: lost (%)", "96.4", round(cv("progressive_tax_block", "ref_hidden", "lost_pct"), 1), tol = 0.05)
  add_result("Holdout", "Hidden GSE100051 (ComBat-ref): lost (%)", "81", round(bv("ComBat-ref", "GSE100051", "lost_pct")), tol = 0.5)
  add_result("Holdout", "Hidden GSE100051 (ComBat-ref): |logFC| ratio", "0.52", round(bv("ComBat-ref", "GSE100051", "ratio"), 2), tol = 0.005)
  add_result("Holdout", "Hidden GSE100051 (ComBat-ref): oracle sensitivity (%)", "95", round(bv("ComBat-ref", "GSE100051", "sens_O")), tol = 0.5)
  add_result("Holdout", "Hidden GSE37901 (ComBat-ref): lost (%)", "23", round(bv("ComBat-ref", "GSE37901", "lost_pct")), tol = 0.5)
  add_result("Holdout", "Hidden GSE37901 (ComBat-ref): |logFC| ratio", "0.88", round(bv("ComBat-ref", "GSE37901", "ratio"), 2), tol = 0.005)
  add_result("Holdout", "Hidden GSE100051 (plain ComBat): lost (%)", "78", round(bv("plain ComBat", "GSE100051", "lost_pct")), tol = 0.5)
  add_result("Holdout", "Hidden GSE100051 (plain ComBat): |logFC| ratio", "0.52", round(bv("plain ComBat", "GSE100051", "ratio"), 2), tol = 0.005)
} else {
  add_not_verified("Holdout", "Unified holdout tables", "data/pipeline/holdout/main/{confusion_full_vs_masked,by_dataset}.csv missing")
}
w1f <- file.path(hd, "w1", "w1_summary.csv")
if (file.exists(w1f)) {
  w1 <- read.csv(w1f, stringsAsFactors = FALSE); wv <- function(n) w1$value[w1$name == n]
  add_result("Holdout", "Imputed cells at weight 1: full-run DEGs", "527", wv("degs_w1"))
  add_result("Holdout", "Imputed cells at weight 1: DEGs shared with weight 0", "505", wv("shared"))
  add_result("Holdout", "Imputed cells at weight 1: FDR-only genes", "9109", wv("fdr_only_w1"))
  add_result("Holdout", "Imputed cells at weight 0: FDR-only genes", "6794", wv("fdr_only_w0"))
} else {
  add_not_verified("Holdout", "Imputed-cell weight sensitivity", "data/pipeline/holdout/w1/w1_summary.csv missing")
}
pm <- file.path(hd, "plain", "reference", "softimpute_meta.csv")
if (file.exists(pm)) add_result("Holdout", "Plain ComBat full-run DEGs", "525", read.csv(pm)$n_deg[1]) else
  add_not_verified("Holdout", "Plain ComBat full run", "data/pipeline/holdout/plain/reference/softimpute_meta.csv missing")
ca <- file.path(hd, "main", "claims_article.csv")
if (file.exists(ca)) {
  st <- system2("Rscript", c("scripts/verify_holdout_claims.R"), stdout = TRUE, stderr = FALSE)
  add_result("Holdout", "Article number audit (verify_holdout_claims.R): claims needing attention", "0",
             as.integer(sub(".* ok, (\\d+) need attention", "\\1", tail(st, 1))))
} else {
  add_not_verified("Holdout", "Article number audit", "claims_article.csv missing")
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
  "| Section | Claim | Article | Reproduced | Tolerance | Status |",
  "|---------|-------|---------|------------|-----------|--------|"
)

for (i in seq_len(nrow(results))) {
  report <- c(report, sprintf("| %s | %s | %s | %s | %s | %s |",
    results$section[i], results$claim[i], results$article_value[i],
    results$reproduced_value[i], results$tolerance[i], results$status[i]))
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

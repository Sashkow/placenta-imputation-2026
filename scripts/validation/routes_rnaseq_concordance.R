#!/usr/bin/env Rscript
## Agreement with the RNA-seq study (Prater 2021, table of its DEGs) for three
## routes: softImpute + ComBat-ref (published tables), dissection + ComBat
## (HarmonizR idea with the article's model), and batch in limma.
## Routes without imputation, 4 of 4. Run from the repository root:
##   Rscript scripts/validation/routes_rnaseq_concordance.R --config=config/main.yaml --label=main
##   Rscript scripts/validation/routes_rnaseq_concordance.R --config=config/ga_matched_prater.yaml --label=ga_matched
suppressPackageStartupMessages(library(openxlsx))
HOLDOUT_DEFAULT_OUT <- "data/pipeline/routes"
source("scripts/validation/_holdout_common.R")
LABEL <- get_arg("label", "main")
run_dir <- config$paths$output
config$paths$reference_de <- file.path(run_dir, "difexp_softimpute_combat_ref.tsv")
D <- prepare_data(config); tools <- make_de_tools(D, config)
ph <- config$phenotype; covariates <- ph$covariates %||% character(0)

dissect <- function(X) {
  pres <- sapply(D$datasets, function(ds) rowSums(!is.na(X[, D$ds_cols[[ds]], drop = FALSE])) >= 2)
  key <- apply(pres, 1, function(r) paste(D$datasets[r], collapse = "|"))
  out <- X
  for (k in unique(key)) {
    dsets <- strsplit(k, "|", fixed = TRUE)[[1]]; if (length(dsets) < 2) next
    g <- rownames(X)[key == k]; cols <- unlist(D$ds_cols[dsets]); sub <- X[g, cols, drop = FALSE]
    ok <- complete.cases(sub) & apply(sub, 1, var) > 0; sub <- sub[ok, , drop = FALSE]
    if (nrow(sub) < 2) next
    pd_sub <- D$pd[colnames(sub), , drop = FALSE]
    mod <- combat_model_matrix(pd_sub, config, informative_covariates(pd_sub, covariates))
    ref <- if (REF_BATCH %in% dsets) REF_BATCH else NULL
    res <- tryCatch(suppressMessages(ComBat(sub, factor(pd_sub[[DATASET_ID]]), mod = mod,
                    par.prior = TRUE, prior.plots = FALSE, ref.batch = ref)), error = function(e) NULL)
    if (is.null(res)) { cat("pattern", k, ": ComBat with model failed, left uncorrected\n"); next }
    out[rownames(sub), cols] <- res
  }
  out
}
limma_na <- function(X, batch_in_design) {
  E <- X[, tools$keep_samples, drop = FALSE]; pd_de <- D$pd[colnames(E), , drop = FALSE]
  if (batch_in_design) {
    grp <- pd_de[[ph$group_column]]; ds <- D$ds_of_sample[colnames(E)]
    both <- sapply(D$datasets, function(d)
      rowSums(!is.na(E[, ds == d & grp == ph$baseline, drop = FALSE])) > 0 &
      rowSums(!is.na(E[, ds == d & grp == ph$contrast, drop = FALSE])) > 0)
    E <- E[rowSums(both) > 0, , drop = FALSE]
  }
  suppressWarnings(fit_limma(E, de_design(pd_de, config, informative_covariates(pd_de, covariates),
                                          batch_in_design = batch_in_design)))
}
des <- list(softimpute_combat_ref = read.delim(config$paths$reference_de, stringsAsFactors = FALSE),
            dissection_combat = limma_na(dissect(D$X), FALSE),
            batch_in_limma = limma_na(D$X, TRUE))

prater <- read.xlsx(loadWorkbook("data/references/prater_2021_supp_tables.xlsx"),
                    sheet = "T1 DEGs_results_table_l2fc1")
prater$entrez <- as.character(prater$entrezgene_id)
prater <- prater[!is.na(prater$entrez) & prater$entrez != "NA" & !duplicated(prater$entrez), ]
pl <- setNames(as.numeric(prater$log2FoldChange), prater$entrez)
pp <- setNames(as.numeric(prater$padj), prater$entrez)
psig <- names(pl)[!is.na(pp) & pp < FDR & abs(pl) > LFC]

ref_vis <- setNames(D$presence[, REF_BATCH], D$genes); kds <- setNames(D$n_ds_per_gene, D$genes)
## genes of the RNA-seq table that all three routes tested: one common set
common <- Reduce(intersect, c(lapply(des, function(d) as.character(d$gene[is.finite(d$logFC)])), list(names(pl))))
strata <- list(all_tested_by_route = NULL, common_set = common,
               common_ref_visible = common[ref_vis[common]], common_ref_hidden = common[!ref_vis[common]],
               common_all_datasets = common[kds[common] == length(D$datasets)])
rows <- list()
for (a in names(des)) {
  d <- des[[a]]; l <- setNames(d$logFC, as.character(d$gene)); sig <- as.character(d$gene)[sig_call(d$logFC, d$adj.P.Val)]
  for (s in names(strata)) {
    g <- if (is.null(strata[[s]])) intersect(names(l)[is.finite(l)], names(pl)) else strata[[s]]
    if (length(g) < 3) next
    ss <- intersect(intersect(sig, psig), g)
    rows[[length(rows) + 1]] <- data.frame(cohort = LABEL, route = a, stratum = s,
      genes_tested = sum(is.finite(l)), degs = length(sig), rnaseq_genes = length(g),
      r = cor(l[g], pl[g]), ccc = lin_ccc(l[g], pl[g]), slope = unname(coef(lm(l[g] ~ pl[g]))[2]),
      sign_agree_pct = 100 * mean(sign(l[g]) == sign(pl[g])),
      rnaseq_degs_in_set = length(intersect(psig, g)), confirmed_degs = length(ss),
      confirmed_same_direction = sum(sign(l[ss]) == sign(pl[ss])),
      degs_of_route_in_set_scope = sum(sig %in% g), stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, paste0("rnaseq_concordance_", LABEL, ".csv")), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)

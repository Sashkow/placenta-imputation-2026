#!/usr/bin/env Rscript
##
## Summaries of the unified holdout run:
##
##   output/table_s1.csv / .tex       imputation accuracy, scheme x imputer
##                                    (+ per-gene MAE by reference-batch stratum)
##   output/table_s2.csv / .tex       DE-call accuracy, scheme x imputer x stratum
##   output/table_s2_by_k.csv         the same by visible-dataset count
##   output/oracle_attribution.csv    ComBat-pathway vs sample-loss split (block schemes)
##   output/derived_sensitivity.csv   expected real-world sensitivity = census-weighted strata
##   output/reference_runs.tex        DEGs per imputer side table
##   output/s1_check_vs_companion.csv softImpute rows vs the companion's cv_comparators rerun
##
## Usage: Rscript scripts/holdout_tables.R [--out=output]

source(file.path({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "."
}, "holdout_common.R"))

read_all <- function(dir, pattern, exclude = NULL) {
  fs <- list.files(file.path(OUT_DIR, dir), pattern = pattern, full.names = TRUE)
  if (!is.null(exclude)) fs <- fs[!grepl(exclude, basename(fs))]
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
}
msd <- function(x, d = 3) {
  x <- x[is.finite(x)]
  if (!length(x)) return("--")
  if (length(x) == 1) sprintf(paste0("%.", d, "f"), x) else
    sprintf(paste0("$%.", d, "f \\pm %.", d, "f$"), mean(x), sd(x))
}
pct <- function(x, d = 1) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), 100 * x), "--")
fmt <- function(x, d = 3) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "--")
order_by <- function(df, col, levels) df[order(match(df[[col]], levels)), ]

schemes_all <- unlist(HCFG$mask_types)
methods_all <- unlist(HCFG$methods)

## ======================================================== values (Table S1)
vals <- read_all("values", "\\.csv$", exclude = "_per_gene")
pg   <- read_all("values", "_per_gene\\.csv$")
stopifnot(!is.null(vals))

s1 <- do.call(rbind, lapply(split(vals, list(vals$scheme, vals$method), drop = TRUE), function(g) {
  ok <- g$converged %in% TRUE
  data.frame(scheme = g$scheme[1], method = g$method[1],
             scheme_label = SCHEME_LABELS[g$scheme[1]], method_label = METHOD_LABELS[g$method[1]],
             n_repeats = nrow(g), converged = sum(ok), cells_masked = round(mean(g$n_masked)),
             r_mean = mean(g$pearson_r[ok]), r_sd = sd(g$pearson_r[ok]),
             rmse_mean = mean(g$rmse[ok]), rmse_sd = sd(g$rmse[ok]),
             mae_mean = mean(g$mae[ok]), mae_sd = sd(g$mae[ok]),
             imputation_seconds = round(mean(g$imputation_seconds)),
             combat_ok = sum(g$combat_status %in% "ok"),
             stringsAsFactors = FALSE, row.names = NULL)
}))

## per-gene MAE by reference-batch stratum (median per repeat, then mean +/- SD)
if (!is.null(pg)) {
  pgs <- do.call(rbind, lapply(split(pg, list(pg$scheme, pg$method, pg$rep), drop = TRUE), function(g)
    data.frame(scheme = g$scheme[1], method = g$method[1], rep = g$rep[1],
               n_genes = nrow(g),
               n_ref_visible = sum(g$ref_batch_visible), n_ref_hidden = sum(!g$ref_batch_visible),
               med_mae_all = median(g$mae, na.rm = TRUE),
               med_mae_ref_visible = if (any(g$ref_batch_visible)) median(g$mae[g$ref_batch_visible], na.rm = TRUE) else NA,
               med_mae_ref_hidden  = if (any(!g$ref_batch_visible)) median(g$mae[!g$ref_batch_visible], na.rm = TRUE) else NA,
               med_r_all = median(g$r, na.rm = TRUE),
               stringsAsFactors = FALSE, row.names = NULL)))
  write.csv(pgs, file.path(OUT_DIR, "values_per_gene_by_repeat.csv"), row.names = FALSE)
  agg <- do.call(rbind, lapply(split(pgs, list(pgs$scheme, pgs$method), drop = TRUE), function(g)
    data.frame(scheme = g$scheme[1], method = g$method[1],
               pg_n_genes = round(mean(g$n_genes)),
               pg_n_ref_visible = round(mean(g$n_ref_visible)), pg_n_ref_hidden = round(mean(g$n_ref_hidden)),
               pg_med_mae_all = mean(g$med_mae_all), pg_med_mae_all_sd = sd(g$med_mae_all),
               pg_med_mae_ref_visible = mean(g$med_mae_ref_visible), pg_med_mae_ref_visible_sd = sd(g$med_mae_ref_visible),
               pg_med_mae_ref_hidden = mean(g$med_mae_ref_hidden), pg_med_mae_ref_hidden_sd = sd(g$med_mae_ref_hidden),
               stringsAsFactors = FALSE, row.names = NULL)))
  s1 <- merge(s1, agg, by = c("scheme", "method"), all.x = TRUE)
}
s1 <- s1[order(match(s1$scheme, schemes_all), match(s1$method, methods_all)), ]
write.csv(s1, file.path(OUT_DIR, "table_s1.csv"), row.names = FALSE)

## LaTeX Table S1
tex <- c("\\begin{tabular}{@{}llrccc@{}}", "\\toprule",
         "\\textbf{Masking} & \\textbf{Method} & \\textbf{Cells masked} & \\textbf{Pearson $r$} & \\textbf{RMSE} & \\textbf{MAE} \\\\",
         "\\midrule")
for (sc in unique(s1$scheme)) {
  d <- s1[s1$scheme == sc, ]
  for (i in seq_len(nrow(d))) {
    lab <- if (i == 1) SCHEME_LABELS[sc] else ""
    g <- vals[vals$scheme == sc & vals$method == d$method[i] & vals$converged %in% TRUE, ]
    tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s \\\\", lab, d$method_label[i],
                          format(d$cells_masked[i], big.mark = "{,}"),
                          msd(g$pearson_r), msd(g$rmse), msd(g$mae)))
  }
  if (sc != tail(unique(s1$scheme), 1)) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(OUT_DIR, "table_s1.tex"))

## per-gene stratified table (block schemes)
if (!is.null(pg)) {
  d <- s1[s1$scheme != "random_cells" & !is.na(s1$pg_n_genes), ]
  tex <- c("\\begin{tabular}{@{}llrrcc@{}}", "\\toprule",
           "\\textbf{Masking} & \\textbf{Method} & \\multicolumn{2}{c}{\\textbf{Genes}} & \\multicolumn{2}{c}{\\textbf{Median per-gene MAE}} \\\\",
           "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
           " & & \\textbf{GSE100051 vis.} & \\textbf{GSE100051 hid.} & \\textbf{GSE100051 vis.} & \\textbf{GSE100051 hid.} \\\\",
           "\\midrule")
  for (sc in unique(d$scheme)) {
    dd <- d[d$scheme == sc, ]
    for (i in seq_len(nrow(dd))) {
      g <- pgs[pgs$scheme == sc & pgs$method == dd$method[i], ]
      tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s \\\\",
                            if (i == 1) SCHEME_LABELS[sc] else "", dd$method_label[i],
                            format(dd$pg_n_ref_visible[i], big.mark = "{,}"),
                            format(dd$pg_n_ref_hidden[i], big.mark = "{,}"),
                            msd(g$med_mae_ref_visible, 2), msd(g$med_mae_ref_hidden, 2)))
    }
    if (sc != tail(unique(d$scheme), 1)) tex <- c(tex, "\\addlinespace")
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}")
  writeLines(tex, file.path(OUT_DIR, "table_s1_per_gene.tex"))
}

## softImpute rows vs the companion's mask-matched rerun
cvp <- file.path(COMPANION, config$paths$cv_comparators)
if (file.exists(cvp) && any(s1$method == "softimpute")) {
  cv <- read.csv(cvp, stringsAsFactors = FALSE)
  cv$r <- as.numeric(sub(" .*", "", cv$pearson_r))
  chk <- s1[s1$method == "softimpute", c("scheme", "r_mean", "rmse_mean")]
  chk$companion_scheme <- ifelse(chk$scheme == "progressive_tax_block", "gene_dataset_block_hard", chk$scheme)
  chk$companion_r <- cv$r[match(paste(chk$companion_scheme, "softimpute"), paste(cv$mask_type, cv$method))]
  chk$diff <- chk$r_mean - chk$companion_r
  chk$note <- ifelse(chk$scheme == "progressive_tax_block",
                     "companion ran gene_dataset_block_hard under the worst-case label; not the same scheme", "")
  write.csv(chk, file.path(OUT_DIR, "s1_check_vs_companion.csv"), row.names = FALSE)
  cat("\n=== softImpute r vs companion cv_comparators rerun ===\n"); print(chk, row.names = FALSE, digits = 4)
}

## ======================================================== calls (Table S2)
de <- read_all("de", "\\.csv$")
if (is.null(de)) { cat("No DE outputs yet; stopping after Table S1.\n"); quit(save = "no") }
de$tested <- !de$dropped_masked & !de$nonfinite_masked & is.finite(de$logFC_full)

summarise_stratum <- function(d, label) {
  tp <- d$sig_full %in% TRUE
  tested <- d$tested
  okO <- tested & is.finite(d$logFC_oracle)
  n_true <- sum(tp)
  rec_tested <- sum(tp & tested & d$sig_masked %in% TRUE)
  ## per-repeat sensitivity for the SD
  sr <- sapply(split(d, d$rep), function(x) { p <- x$sig_full %in% TRUE
    if (sum(p)) mean(x$sig_masked[p] %in% TRUE) else NA })
  data.frame(
    stratum = label, n_scored = nrow(d), n_tested = sum(tested),
    n_dropped = sum(d$dropped_masked), n_nonfinite = sum(d$nonfinite_masked),
    n_true_deg = n_true, n_true_deg_tested = sum(tp & tested), n_recovered = rec_tested,
    sensitivity_tested = if (sum(tp & tested)) rec_tested / sum(tp & tested) else NA,
    sensitivity_all = if (n_true) rec_tested / n_true else NA,
    sensitivity_sd_over_repeats = if (sum(!is.na(sr)) > 1) sd(sr, na.rm = TRUE) else NA,
    n_true_nondeg_tested = sum(!tp & tested),
    n_false_positive = sum(!tp & tested & d$sig_masked %in% TRUE),
    specificity = if (sum(!tp & tested)) mean(!(d$sig_masked[!tp & tested] %in% TRUE)) else NA,
    direction_flip = if (sum(tp & tested)) mean(sign(d$logFC_masked[tp & tested]) != sign(d$logFC_full[tp & tested])) else NA,
    logfc_mae_MF = if (any(tested)) mean(abs(d$logFC_masked[tested] - d$logFC_full[tested])) else NA,
    logfc_medae_MF = if (any(tested)) median(abs(d$logFC_masked[tested] - d$logFC_full[tested])) else NA,
    logfc_r_MF = if (sum(tested) > 2) cor(d$logFC_masked[tested], d$logFC_full[tested]) else NA,
    logfc_ccc_MF = if (sum(tested) > 2) lin_ccc(d$logFC_masked[tested], d$logFC_full[tested]) else NA,
    deg_logfc_shrink = { p <- tp & tested
      if (sum(p)) median(abs(d$logFC_masked[p]) / abs(d$logFC_full[p])) else NA },
    n_oracle = sum(okO),
    sensitivity_oracle = if (sum(tp & okO)) mean(d$sig_oracle[tp & okO] %in% TRUE) else NA,
    logfc_mae_OF = if (any(okO)) mean(abs(d$logFC_oracle[okO] - d$logFC_full[okO])) else NA,
    logfc_mae_MO = if (any(okO)) mean(abs(d$logFC_masked[okO] - d$logFC_oracle[okO])) else NA,
    stringsAsFactors = FALSE, row.names = NULL)
}

strata_of <- function(d) {
  out <- list(summarise_stratum(d, "all"))
  for (rb in c(TRUE, FALSE)) {
    x <- d[d$ref_batch_visible == rb, ]
    out[[length(out) + 1]] <- if (nrow(x)) summarise_stratum(x, if (rb) "ref_visible" else "ref_hidden") else
      { z <- summarise_stratum(d[0, ], if (rb) "ref_visible" else "ref_hidden"); z }
  }
  do.call(rbind, out)
}

s2 <- do.call(rbind, lapply(split(de, list(de$scheme, de$method), drop = TRUE), function(g)
  cbind(scheme = g$scheme[1], method = g$method[1],
        scheme_label = SCHEME_LABELS[g$scheme[1]], method_label = METHOD_LABELS[g$method[1]],
        strata_of(g), stringsAsFactors = FALSE)))
s2 <- s2[order(match(s2$scheme, schemes_all), match(s2$method, methods_all),
               match(s2$stratum, c("all", "ref_visible", "ref_hidden"))), ]
rownames(s2) <- NULL
write.csv(s2, file.path(OUT_DIR, "table_s2.csv"), row.names = FALSE)

s2k <- do.call(rbind, lapply(split(de, list(de$scheme, de$method, de$k_visible, de$ref_batch_visible), drop = TRUE),
  function(g) cbind(scheme = g$scheme[1], method = g$method[1], k_visible = g$k_visible[1],
                    ref_batch_visible = g$ref_batch_visible[1],
                    summarise_stratum(g, paste0("k=", g$k_visible[1])), stringsAsFactors = FALSE)))
s2k <- s2k[order(match(s2k$scheme, schemes_all), match(s2k$method, methods_all), s2k$k_visible, !s2k$ref_batch_visible), ]
rownames(s2k) <- NULL
write.csv(s2k, file.path(OUT_DIR, "table_s2_by_k.csv"), row.names = FALSE)

## oracle attribution (block schemes): share of the M-vs-F |logFC| discrepancy
## carried by the ComBat pathway (M vs O) versus sample loss (O vs F)
oa <- s2[s2$stratum == "all" & is.finite(s2$logfc_mae_MO), ]
oa <- data.frame(scheme = oa$scheme, method = oa$method, n_genes = oa$n_oracle,
                 mae_MF = oa$logfc_mae_MF, mae_OF = oa$logfc_mae_OF, mae_MO = oa$logfc_mae_MO,
                 ratio_combat_to_sampleloss = oa$logfc_mae_MO / oa$logfc_mae_OF,
                 share_combat_pathway = oa$logfc_mae_MO / (oa$logfc_mae_MO + oa$logfc_mae_OF),
                 sensitivity_masked = oa$sensitivity_tested, sensitivity_oracle = oa$sensitivity_oracle,
                 stringsAsFactors = FALSE)
write.csv(oa, file.path(OUT_DIR, "oracle_attribution.csv"), row.names = FALSE)

## Three-arm comparison for softImpute: Masked vs Full, Oracle vs Full, Masked vs
## Oracle, per scheme and reference-batch stratum (table_s3_arms)
arms <- s2[s2$method == "softimpute", ]
arms_out <- data.frame(scheme = arms$scheme, stratum = arms$stratum, n_scored = arms$n_scored, n_dropped = arms$n_dropped,
                       n_true_deg = arms$n_true_deg, r_MF = arms$logfc_r_MF, mae_MF = arms$logfc_mae_MF,
                       mae_OF = arms$logfc_mae_OF, mae_MO = arms$logfc_mae_MO,
                       sens_M = arms$sensitivity_all, sens_O = arms$sensitivity_oracle, spec_M = arms$specificity,
                       stringsAsFactors = FALSE)
write.csv(arms_out, file.path(OUT_DIR, "table_s3_arms.csv"), row.names = FALSE)
STRATUM_LABELS <- c(all = "all masked genes", ref_visible = "GSE100051 visible", ref_hidden = "GSE100051 hidden")
tex <- c("\\begingroup\\setlength{\\tabcolsep}{3.5pt}", "\\begin{tabular}{@{}lrrrrrrrrr@{}}", "\\toprule",
         "\\textbf{Stratum} & \\textbf{Genes} & \\textbf{True} & \\multicolumn{4}{c}{\\textbf{log fold-change agreement between runs}} & \\multicolumn{2}{c}{\\textbf{Sens.\\ \\%}} & \\textbf{Spec.\\ \\%} \\\\",
         "\\cmidrule(lr){4-7}\\cmidrule(lr){8-9}",
         " & & \\textbf{DEGs} & \\textbf{$r$ M,F} & \\textbf{MAE M,F} & \\textbf{MAE O,F} & \\textbf{MAE M,O} & \\textbf{M} & \\textbf{O} & \\textbf{M} \\\\",
         "\\midrule")
for (sc in unique(arms_out$scheme)) {
  d <- arms_out[arms_out$scheme == sc, ]
  a0 <- d[d$stratum == "all", ]
  tex <- c(tex, sprintf("\\multicolumn{10}{@{}l}{\\textit{%s} (%s dropped before testing)} \\\\[2pt]", SCHEME_LABELS[sc], format(a0$n_dropped, big.mark = "{,}")))
  for (st in c("all", "ref_visible", "ref_hidden")) {
    x <- d[d$stratum == st, ]
    tex <- c(tex, sprintf("%s & %s & %d & %s & %s & %s & %s & %s & %s & %s \\\\", STRATUM_LABELS[st],
                          format(x$n_scored, big.mark = "{,}"), x$n_true_deg, fmt(x$r_MF), fmt(x$mae_MF), fmt(x$mae_OF), fmt(x$mae_MO),
                          pct(x$sens_M), pct(x$sens_O), pct(x$spec_M)))
  }
  if (sc != tail(unique(arms_out$scheme), 1)) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\endgroup")
writeLines(tex, file.path(OUT_DIR, "table_s3_arms.tex"))

## LaTeX Table S2: one header line per scheme (gene-level results, dropped),
## then one row per method with the two reference-batch strata
tex <- c("\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
         "\\textbf{Method} & \\multicolumn{3}{c}{\\textbf{GSE100051 visible}} & \\multicolumn{3}{c}{\\textbf{GSE100051 hidden}} \\\\",
         "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
         " & \\textbf{true DEGs} & \\textbf{sens.\\ \\%} & \\textbf{spec.\\ \\%} & \\textbf{true DEGs} & \\textbf{sens.\\ \\%} & \\textbf{spec.\\ \\%} \\\\",
         "\\midrule")
for (sc in unique(s2$scheme)) {
  d <- s2[s2$scheme == sc, ]
  a0 <- d[d$stratum == "all", ][1, ]
  tex <- c(tex, sprintf("\\multicolumn{7}{@{}l}{\\textit{%s: %s gene-level results, %s dropped before testing}} \\\\[2pt]",
                        SCHEME_LABELS[sc], format(a0$n_scored, big.mark = "{,}"), format(a0$n_dropped, big.mark = "{,}")))
  for (m in unique(d$method)) {
    v <- d[d$method == m & d$stratum == "ref_visible", ]
    h <- d[d$method == m & d$stratum == "ref_hidden", ]
    tex <- c(tex, sprintf("%s & %d & %s & %s & %d & %s & %s \\\\", METHOD_LABELS[m],
                          v$n_true_deg, pct(v$sensitivity_all), pct(v$specificity),
                          h$n_true_deg, pct(h$sensitivity_all), pct(h$specificity)))
  }
  if (sc != tail(unique(s2$scheme), 1)) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(OUT_DIR, "table_s2.tex"))

## ======================================================== derived sensitivity
census <- read.csv(file.path(COMPANION, config$paths$census), stringsAsFactors = FALSE)
ni <- census[census$set == "non_intersection", ]
p_ref <- sum(ni$n_genes[ni$has_ref_batch]) / sum(ni$n_genes)
rows <- list()
for (m in unique(s2$method)) for (sc in setdiff(unique(s2$scheme), "random_cells")) {
  v <- s2[s2$scheme == sc & s2$method == m & s2$stratum == "ref_visible", ]
  h <- s2[s2$scheme == sc & s2$method == m & s2$stratum == "ref_hidden", ]
  for (kind in c("sensitivity_all", "sensitivity_tested")) {
    sv <- v[[kind]]; sh <- h[[kind]]
    rows[[length(rows) + 1]] <- data.frame(
      method = m, scheme = sc, sensitivity_kind = kind,
      real_gene_frac_ref_visible = p_ref, real_gene_frac_ref_hidden = 1 - p_ref,
      census_source = config$paths$census,
      sens_ref_visible = sv, sens_ref_hidden = sh,
      n_true_deg_ref_visible = v$n_true_deg, n_true_deg_ref_hidden = h$n_true_deg,
      expected_sensitivity = p_ref * sv + (1 - p_ref) * ifelse(is.finite(sh), sh, 0),
      stringsAsFactors = FALSE)
  }
}
ds <- do.call(rbind, rows)
write.csv(ds, file.path(OUT_DIR, "derived_sensitivity.csv"), row.names = FALSE)

## ======================================================== reference side table
rp <- file.path(OUT_DIR, "reference_runs_summary.csv")
if (file.exists(rp)) {
  rr <- read.csv(rp, stringsAsFactors = FALSE)
  rr <- rr[order(match(rr$method, methods_all)), ]
  tex <- c("\\begin{tabular}{@{}lrrrr@{}}", "\\toprule",
           "\\textbf{Imputer} & \\textbf{DEGs} & \\textbf{shared with softImpute} & \\textbf{Jaccard} & \\textbf{logFC CCC} \\\\",
           "\\midrule",
           sprintf("%s & %d & %d & %s & %s \\\\", METHOD_LABELS[rr$method], rr$n_deg, rr$shared_with_softimpute,
                   fmt(rr$jaccard_vs_softimpute), fmt(rr$logfc_ccc_vs_softimpute)),
           "\\bottomrule", "\\end{tabular}")
  writeLines(tex, file.path(OUT_DIR, "reference_runs.tex"))
}

## ======================================================== console summary
cat("\n=== Table S1 (values) ===\n")
print(s1[, c("scheme_label", "method_label", "cells_masked", "r_mean", "rmse_mean", "mae_mean", "converged", "combat_ok")],
      row.names = FALSE, digits = 3)
cat("\n=== Table S2 (calls), reference-batch strata ===\n")
print(s2[, c("scheme_label", "method_label", "stratum", "n_scored", "n_dropped", "n_true_deg", "n_recovered",
             "sensitivity_all", "specificity", "logfc_mae_MF", "deg_logfc_shrink")], row.names = FALSE, digits = 3)
cat("\n=== Oracle attribution ===\n"); print(oa, row.names = FALSE, digits = 3)
cat("\n=== Derived real-world sensitivity (softImpute) ===\n")
print(ds[ds$method == "softimpute", c("scheme", "sensitivity_kind", "real_gene_frac_ref_visible",
                                      "sens_ref_visible", "sens_ref_hidden", "expected_sensitivity")],
      row.names = FALSE, digits = 3)
cat("\nWritten to:", OUT_DIR, "\n")

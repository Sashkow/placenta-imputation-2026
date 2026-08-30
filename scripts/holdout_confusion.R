#!/usr/bin/env Rscript
##
## Full vs Masked 2x2 DEG-call matrices (stayed DEG / DEG -> not / not -> DEG /
## stayed not) for one imputer, per masking scheme, overall and by GSE100051
## presence. Counts are computed per repeat and reported as mean +/- SD over
## the three repeats, i.e. on the scale of one run. Dropped genes (no real
## value left in one trimester) count as "not DEG" in the masked arm.
##
## Outputs (in the run's output directory):
##   confusion_full_vs_masked.csv   all schemes x strata (+ oracle sensitivity, |logFC| ratio, CCC M vs F)
##   confusion.tex                  supplement table: three schemes x {all, GSE100051 visible, hidden}
##   confusion_main.tex             main-text table: block scheme strata + by-dataset rows
##                                  (requires by_dataset.csv from holdout_by_dataset.R)
##   w1_summary.csv (--w1_ref=...)  Full(w=1) vs Full(w=0) DEG lists and FDR-only counts
##
## Usage: Rscript scripts/holdout_confusion.R [--method=softimpute] [--out=output]
##        Rscript scripts/holdout_confusion.R --config=config/config_holdout_w1.yaml --w1_ref=output/reference/softimpute_difexp.csv

source(file.path({ fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "." }, "holdout_common.R"))
METHOD <- get_arg("method", "softimpute")
STRATUM_LABELS <- c(overall = "all masked genes", ref_visible = "GSE100051 visible", ref_hidden = "GSE100051 hidden")

D <- prepare_data(config)
rows <- list()
for (sc in MASK_TYPES) {
  per <- list()
  for (rep in 1:3) {
    de <- read.csv(file.path(OUT_DIR, "de", sprintf("%s_rep%d_%s.csv", sc, rep, METHOD)), stringsAsFactors = FALSE)
    de$rep <- rep; per[[rep]] <- de
  }
  d <- do.call(rbind, per); d <- d[is.finite(d$logFC_full), ]
  d$full <- d$sig_full %in% TRUE; d$masked <- d$sig_masked %in% TRUE; d$oracle <- d$sig_oracle %in% TRUE
  d$tested <- !d$dropped_masked & is.finite(d$logFC_masked)
  strata <- list(overall = rep(TRUE, nrow(d)), ref_visible = d$ref_batch_visible, ref_hidden = !d$ref_batch_visible)
  for (st in names(strata)) {
    pr <- do.call(rbind, lapply(1:3, function(r) { g <- d[strata[[st]] & d$rep == r, ]; tp <- g$full
      data.frame(genes = nrow(g), dropped = sum(g$dropped_masked), true_deg = sum(tp),
                 stayed_deg = sum(tp & g$masked), deg_to_not = sum(tp & !g$masked),
                 not_to_deg = sum(!tp & g$masked), stayed_not = sum(!tp & !g$masked),
                 ratio = { p <- tp & g$tested; if (sum(p)) median(abs(g$logFC_masked[p]) / abs(g$logFC_full[p])) else NA },
                 ccc_MF = if (sum(g$tested) > 2) lin_ccc(g$logFC_masked[g$tested], g$logFC_full[g$tested]) else NA,
                 sens_O = if (sum(tp)) 100 * mean(g$oracle[tp]) else NA) }))
    m <- colMeans(pr, na.rm = TRUE); sdv <- apply(pr, 2, sd, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(scheme = sc, stratum = st, n_repeats = nrow(pr),
      as.list(setNames(c(rbind(m, sdv)), paste0(rep(names(m), each = 2), c("", "_sd")))),
      lost_pct = 100 * m["deg_to_not"] / (m["stayed_deg"] + m["deg_to_not"]),
      gained_pct = 100 * m["not_to_deg"] / (m["not_to_deg"] + m["stayed_not"]), row.names = NULL) }
}
res <- do.call(rbind, rows); write.csv(res, file.path(OUT_DIR, "confusion_full_vs_masked.csv"), row.names = FALSE)

pm <- function(m, s, d = 0) ifelse(is.finite(m), sprintf(paste0("%.", d, "f $\\pm$ %.", d, "f"), m, s), "--")
pc <- function(x, d = 1) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "--")
big <- function(x) format(round(x), big.mark = "{,}")

## ---- supplement table: three schemes
tex <- c("\\begingroup\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{@{}lrrrrrrr@{}}", "\\toprule",
         "\\textbf{Stratum} & \\textbf{Genes} & \\textbf{Stayed DEG} & \\textbf{DEG $\\to$ not} & \\textbf{Not $\\to$ DEG} & \\textbf{Stayed not} & \\textbf{Lost \\%} & \\textbf{Gained \\%} \\\\",
         "\\midrule")
for (sc in MASK_TYPES) {
  d <- res[res$scheme == sc, ]; a0 <- d[d$stratum == "overall", ]
  tex <- c(tex, sprintf("\\multicolumn{8}{@{}l}{\\textit{%s} (%s genes per run; %s dropped before testing)} \\\\[2pt]",
                        SCHEME_LABELS[sc], big(a0$genes), big(a0$dropped)))
  for (st in c("overall", "ref_visible", "ref_hidden")) { x <- d[d$stratum == st, ]
    tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\", STRATUM_LABELS[st], big(x$genes),
                          pm(x$stayed_deg, x$stayed_deg_sd), pm(x$deg_to_not, x$deg_to_not_sd), pm(x$not_to_deg, x$not_to_deg_sd),
                          pm(x$stayed_not, x$stayed_not_sd), pc(x$lost_pct), pc(x$gained_pct, 2))) }
  if (sc != tail(MASK_TYPES, 1)) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\endgroup"); writeLines(tex, file.path(OUT_DIR, "confusion.tex"))

## ---- main-text table: block scheme strata + by-dataset rows
bd_file <- file.path(OUT_DIR, "by_dataset.csv")
if (file.exists(bd_file) && "gene_dataset_block" %in% MASK_TYPES) {
  bd <- read.csv(bd_file, stringsAsFactors = FALSE); bd <- bd[bd$variant == bd$variant[1], ]
  d <- res[res$scheme == "gene_dataset_block", ]
  tex <- c("\\begingroup\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}lrrrrrrrr@{}}", "\\toprule",
           "\\textbf{Genes} & \\textbf{$n$} & \\textbf{True DEGs} & \\textbf{Stayed} & \\textbf{Lost} & \\textbf{Gained} & \\textbf{Lost \\%} & \\textbf{CCC} & \\textbf{Oracle \\%} \\\\",
           "\\midrule")
  for (st in c("overall", "ref_visible", "ref_hidden")) { x <- d[d$stratum == st, ]
    tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", STRATUM_LABELS[st], big(x$genes), pm(x$true_deg, x$true_deg_sd),
                          pm(x$stayed_deg, x$stayed_deg_sd), pm(x$deg_to_not, x$deg_to_not_sd), pm(x$not_to_deg, x$not_to_deg_sd),
                          pc(x$lost_pct), pc(x$ccc_MF, 2), pc(x$sens_O))) }
  n_one <- sum(bd$genes); n_all <- d$genes[d$stratum == "overall"]
  tex <- c(tex, "\\addlinespace", sprintf("\\multicolumn{9}{@{}l}{\\textit{genes that lost a single block (%s of %s per run), by the dataset of that block (1st/2nd-trim.\\ samples)}} \\\\[2pt]", big(n_one), big(n_all)))
  for (j in seq_len(nrow(bd))) tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", bd$label[j], big(bd$genes[j]),
    pm(bd$true_deg[j], bd$true_deg_sd[j]), pm(bd$stayed[j], bd$stayed_sd[j]), pm(bd$lost[j], bd$lost_sd[j]), pm(bd$gained[j], bd$gained_sd[j]),
    pc(bd$lost_pct[j]), pc(bd$ccc_MF[j], 2), pc(bd$sens_O[j])))
  tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\endgroup"); writeLines(tex, file.path(OUT_DIR, "confusion_main.tex"))
}

## ---- w = 1 vs w = 0 summary (only when a w=0 reference DE table is given)
w1_ref <- get_arg("w1_ref")
if (!is.null(w1_ref)) {
  f0 <- read.csv(if (grepl("^/", w1_ref)) w1_ref else file.path(ALT_DIR, w1_ref), stringsAsFactors = FALSE)
  f1 <- read.csv(file.path(OUT_DIR, "reference", paste0(METHOD, "_difexp.csv")), stringsAsFactors = FALSE)
  S0 <- f0$gene[f0$sig]; S1 <- f1$gene[f1$sig]; sh <- intersect(f0$gene, f1$gene)
  only0 <- setdiff(S0, S1); only1 <- setdiff(S1, S0)
  w1 <- data.frame(name = c("degs_w0", "degs_w1", "shared", "only_w0", "only_w1", "only_w0_fdr_sig_in_w1", "only_w1_fdr_sig_in_w0",
                            "fdr_only_w0", "fdr_only_w1", "logfc_r", "median_adjp_ratio_w1_over_w0_lfc_gt1"),
                   value = c(length(S0), length(S1), length(intersect(S0, S1)), length(only0), length(only1),
                             sum(f1$adj.P.Val[match(only0, f1$gene)] < FDR), sum(f0$adj.P.Val[match(only1, f0$gene)] < FDR),
                             sum(f0$adj.P.Val < FDR), sum(f1$adj.P.Val < FDR),
                             cor(f0$logFC[match(sh, f0$gene)], f1$logFC[match(sh, f1$gene)]),
                             { i0 <- match(f1$gene, f0$gene); ok <- abs(f1$logFC) > LFC & abs(f0$logFC[i0]) > LFC; median((f1$adj.P.Val / f0$adj.P.Val[i0])[ok]) }))
  write.csv(w1, file.path(OUT_DIR, "w1_summary.csv"), row.names = FALSE); cat("\n=== w=1 vs w=0 full runs ===\n"); print(w1, row.names = FALSE)
}

cat(sprintf("\n=== %s: Full vs Masked (mean ± SD over 3 repeats, per-run scale) ===\n", METHOD_LABELS[METHOD]))
print(res[, c("scheme", "stratum", "genes", "dropped", "true_deg", "stayed_deg", "deg_to_not", "not_to_deg", "lost_pct", "gained_pct", "ratio", "sens_O")], row.names = FALSE, digits = 3)
cat("\nWritten to:", OUT_DIR, "\n")

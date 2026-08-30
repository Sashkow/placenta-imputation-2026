#!/usr/bin/env Rscript
##
## Split the gene-dataset block holdout by WHICH dataset was hidden.
##
## For genes with exactly one dataset hidden (k_visible == 5) the Full/Masked/
## Oracle calls are grouped by the hidden dataset, per repeat, and reported as
## mean +/- SD over repeats: true DEGs, stayed, lost, gained, lost %, median
## masked/full |logFC| ratio among true DEGs, oracle sensitivity, Masked-vs-
## Oracle logFC MAE. This separates "second-trimester samples were lost"
## (GSE37901, a second-trimester-only batch) from "the dominant batch was
## imputed" (GSE100051: 49 of 117 samples, 7 of 15 second-trimester samples).
##
## Several output directories can be compared side by side (e.g. ComBat-ref in
## output/ and plain ComBat in output_plain/); each has its own Full arm.
##
## Usage:
##   Rscript scripts/holdout_by_dataset.R --dirs=output,output_plain --labels=ComBat-ref,plain\ ComBat
## Writes by_dataset.csv and by_dataset.tex into the FIRST directory.

source(file.path({ fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "." }, "holdout_common.R"))

dirs <- split_arg("dirs", "output"); labels <- split_arg("labels", dirs)
dirs <- sapply(dirs, function(d) if (grepl("^/", d)) d else file.path(ALT_DIR, d))
METHOD <- get_arg("method", "softimpute"); SCHEME <- get_arg("scheme", "gene_dataset_block")

D <- prepare_data(config)
tr <- table(D$ds_of_sample, D$grp_of_sample)
ds_label <- setNames(sprintf("%s (%d/%d)", rownames(tr), tr[, config$phenotype$baseline], tr[, config$phenotype$contrast]), rownames(tr))

per_dir <- list()
for (i in seq_along(dirs)) {
  rows <- list()
  for (rep in 1:3) {
    mk <- readRDS(file.path(dirs[i], "masks", sprintf("%s_rep%d.rds", SCHEME, rep)))
    M <- matrix(FALSE, nrow(D$X), ncol(D$X), dimnames = dimnames(D$X)); M[mk$mask_idx] <- TRUE
    hid <- sapply(D$datasets, function(ds) rowSums(M[, D$ds_cols[[ds]], drop = FALSE]) > 0)
    de <- read.csv(file.path(dirs[i], "de", sprintf("%s_rep%d_%s.csv", SCHEME, rep, METHOD)), stringsAsFactors = FALSE)
    de$gene <- as.character(de$gene); one <- de[de$k_visible == 5 & is.finite(de$logFC_full), ]
    one$hidden_ds <- apply(hid[match(one$gene, rownames(hid)), , drop = FALSE], 1, function(r) D$datasets[r][1])
    one$F <- one$sig_full %in% TRUE; one$Mc <- one$sig_masked %in% TRUE; one$O <- one$sig_oracle %in% TRUE
    one$tested <- !one$dropped_masked & is.finite(one$logFC_masked)
    for (ds in D$datasets) { g <- one[one$hidden_ds == ds, ]; tp <- g$F
      rows[[length(rows) + 1]] <- data.frame(rep = rep, hidden = ds, genes = nrow(g), true_deg = sum(tp),
        stayed = sum(tp & g$Mc), lost = sum(tp & !g$Mc), gained = sum(!tp & g$Mc),
        lost_pct = if (sum(tp)) 100 * mean(!g$Mc[tp]) else NA,
        ratio = { p <- tp & g$tested; if (sum(p)) median(abs(g$logFC_masked[p]) / abs(g$logFC_full[p])) else NA },
        sens_O = if (sum(tp)) 100 * mean(g$O[tp]) else NA,
        ccc_MF = if (sum(g$tested) > 2) lin_ccc(g$logFC_masked[g$tested], g$logFC_full[g$tested]) else NA,
        mae_MO = mean(abs(g$logFC_masked[g$tested] - g$logFC_oracle[g$tested]), na.rm = TRUE)) }
  }
  pr <- do.call(rbind, rows)
  agg <- do.call(rbind, lapply(split(pr, pr$hidden), function(x) {
    m <- colMeans(x[, -(1:2)], na.rm = TRUE); s <- apply(x[, -(1:2)], 2, sd, na.rm = TRUE)
    data.frame(variant = labels[i], hidden = x$hidden[1], label = ds_label[x$hidden[1]],
               as.list(setNames(c(rbind(m, s)), paste0(rep(names(m), each = 2), c("", "_sd")))), row.names = NULL) }))
  per_dir[[i]] <- agg
}
res <- do.call(rbind, per_dir)
## order: by first variant's lost_pct descending
ord <- res[res$variant == labels[1], ]; ds_order <- ord$hidden[order(-ord$lost_pct)]
res <- res[order(match(res$variant, labels), match(res$hidden, ds_order)), ]
write.csv(res, file.path(dirs[1], "by_dataset.csv"), row.names = FALSE)

pm <- function(m, s, d = 0) ifelse(is.finite(m), sprintf(paste0("%.", d, "f $\\pm$ %.", d, "f"), m, s), "--")
pc <- function(x, d = 1) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "--")
## LaTeX: one block per variant; columns true DEGs | lost | gained | lost % | ratio | oracle sens %
tex <- c("\\begingroup\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
         "\\textbf{Hidden dataset (1st/2nd)} & \\textbf{True DEGs} & \\textbf{Lost} & \\textbf{Gained} & \\textbf{Lost \\%} & \\textbf{$|$logFC$|$ ratio} & \\textbf{Oracle sens.\\ \\%} \\\\",
         "\\midrule")
for (i in seq_along(labels)) {
  d <- res[res$variant == labels[i], ]
  if (length(labels) > 1) tex <- c(tex, sprintf("\\multicolumn{7}{@{}l}{\\textit{%s}} \\\\[2pt]", labels[i]))
  for (j in seq_len(nrow(d))) tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", d$label[j],
    pm(d$true_deg[j], d$true_deg_sd[j]), pm(d$lost[j], d$lost_sd[j]), pm(d$gained[j], d$gained_sd[j]),
    pc(d$lost_pct[j]), pc(d$ratio[j], 2), pc(d$sens_O[j])))
  if (i < length(labels)) tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\endgroup")
writeLines(tex, file.path(dirs[1], "by_dataset.tex"))

cat(sprintf("\n=== %s, %s: by single hidden dataset (mean ± SD over 3 repeats) ===\n", SCHEME, METHOD))
print(res[, c("variant", "label", "genes", "true_deg", "lost", "gained", "lost_pct", "ratio", "sens_O", "mae_MO")], row.names = FALSE, digits = 3)
cat("\nWritten:", file.path(dirs[1], "by_dataset.csv"), "\n")

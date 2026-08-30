#!/usr/bin/env Rscript
## Every DE-derived number written into the alternative main.tex / supplement.tex
## during the post-fix sweep, computed from the regenerated companion outputs and
## written as name/value rows to output/de_sweep_numbers.csv for the claims audit.
## Usage: Rscript scripts/de_sweep_numbers.R [--companion=PATH] [--out=output]
source(file.path({ fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "." }, "holdout_common.R"))
suppressMessages(library(openxlsx))
q <- function(f) read.delim(f, stringsAsFactors = FALSE)
ccc <- function(x, y) 2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
sig <- function(x) x$gene[x$adj.P.Val < FDR & abs(x$logFC) > LFC]
rows <- list(); add <- function(name, value) rows[[length(rows) + 1]] <<- data.frame(name = name, value = value)
P <- function(...) file.path(COMPANION, ...)
d <- q(P("data/pipeline/main/difexp_softimpute_combat_ref.tsv")); n <- q(P("data/pipeline/main/difexp_none_combat_ref.tsv"))
nc <- q(P("data/pipeline/main/difexp_softimpute_combat.tsv")); b <- q(P("data/pipeline/balanced/difexp_softimpute_combat_ref.tsv"))
o <- q(P("data/pipeline/batch_in_limma/difexp_softimpute_batch_in_limma.tsv"))
d$gene <- as.character(d$gene); n$gene <- as.character(n$gene); b$gene <- as.character(b$gene); o$gene <- as.character(o$gene); nc$gene <- as.character(nc$gene)
S <- sig(d); N <- sig(n); B <- sig(b); O <- sig(o); NC <- sig(nc)
add("degs_full", length(S)); add("degs_full_up", sum(d$logFC[d$gene %in% S] > 0)); add("degs_full_down", sum(d$logFC[d$gene %in% S] < 0))
add("degs_combat_noref", length(NC)); add("degs_combat_noref_up", sum(nc$logFC[nc$gene %in% NC] > 0)); add("degs_combat_noref_down", sum(nc$logFC[nc$gene %in% NC] < 0))
add("degs_intersection", length(N)); gained <- setdiff(S, N); add("gained", length(gained))
add("gained_crossing", length(intersect(gained, n$gene))); add("gained_outside", length(setdiff(gained, n$gene))); add("lost_vs_intersection", length(setdiff(N, S)))
cr <- intersect(gained, n$gene); add("crossing_all_fdr_sig_in_intersection", sum(n$adj.P.Val[match(cr, n$gene)] < FDR))
add("crossing_min_abs_logfc_intersection", min(abs(n$logFC[match(cr, n$gene)]))); add("crossing_max_abs_logfc_intersection", max(abs(n$logFC[match(cr, n$gene)])))
sh <- intersect(d$gene, n$gene); x <- d$logFC[match(sh, d$gene)]; y <- n$logFC[match(sh, n$gene)]
add("shared_genes", length(sh)); add("shared_logfc_r", cor(x, y)); add("shared_logfc_ccc", ccc(x, y)); add("shared_same_direction_pct", 100 * mean(sign(x) == sign(y)))
## balanced reference
add("balanced_degs", length(B)); add("full_in_balanced", length(intersect(S, B))); add("full_in_balanced_pct", 100 * length(intersect(S, B)) / length(S))
add("balanced_jaccard", length(intersect(S, B)) / length(union(S, B))); add("balanced_in_full_pct", 100 * length(intersect(S, B)) / length(B))
u <- union(S, B); add("balanced_union_ccc", ccc(d$logFC[match(u, d$gene)], b$logFC[match(u, b$gene)]))
miss <- setdiff(S, B); bm <- b[match(miss, b$gene), ]; dm <- d[match(miss, d$gene), ]
add("full_only", length(miss)); add("full_only_same_dir", sum(sign(bm$logFC) == sign(dm$logFC))); add("full_only_fdr_sig", sum(bm$adj.P.Val < FDR))
add("full_only_median_abs_logfc", median(abs(bm$logFC))); add("full_only_exceed_0.8", sum(abs(bm$logFC) > 0.8)); add("full_only_fail_both", sum(bm$adj.P.Val >= FDR & abs(bm$logFC) <= LFC))
bo <- setdiff(B, S); add("balanced_only", length(bo)); add("balanced_only_same_dir", sum(sign(b$logFC[match(bo, b$gene)]) == sign(d$logFC[match(bo, d$gene)])))
add("full_fdr_only", sum(d$adj.P.Val < FDR)); add("balanced_fdr_only", sum(b$adj.P.Val < FDR))
## one-step
sho <- intersect(d$gene, o$gene); add("onestep_degs", length(O)); add("onestep_shared", length(intersect(S, O))); add("onestep_logfc_r", cor(d$logFC[match(sho, d$gene)], o$logFC[match(sho, o$gene)]))
t2 <- setdiff(S, O); ob <- o[match(t2, o$gene), ]; add("twostep_only", length(t2)); add("twostep_only_below_lfc", sum(abs(ob$logFC) <= LFC)); add("twostep_only_fail_fdr", sum(ob$adj.P.Val >= FDR))
## subsampling
t1 <- q(P("data/pipeline/validation/test1_first_trim_subsample.tsv")); t1 <- t1[t1$status == "ok", ]
for (N1 in c(10, 20, 30, 40, 100)) { z <- t1[t1$N_1st == N1, ]
  add(paste0("t1_ccc_N", N1), median(z$logfc_ccc_vs_full)); add(paste0("t1_jaccard_N", N1), median(z$jaccard_vs_full))
  add(paste0("t1_retention_N", N1), 100 * median(z$overlap_vs_full)); add(paste0("t1_retention_fdr_only_N", N1), 100 * median(z$overlap_fdr_only)); add(paste0("t1_ndeg_N", N1), median(z$n_deg)) }
t3 <- q(P("data/pipeline/validation/test3_split_half.tsv")); t3 <- t3[t3$status_a == "ok" & t3$status_b == "ok", ]
add("t3_jaccard", median(t3$jaccard_a_vs_b)); add("t3_ccc", median(t3$logfc_ccc_a_vs_b)); add("t3_retention", 100 * median(c(t3$overlap_a_vs_full, t3$overlap_b_vs_full)))
t1b <- q(P("data/pipeline/validation/test1b_vs_balanced.tsv")); t1b <- t1b[t1b$status == "ok", ]
add("t1b_balanced_retention_N10", 100 * median(t1b$overlap_vs_balanced[t1b$N_1st == 10]))
m1 <- aggregate(overlap_vs_full ~ N_1st, t1, median); m2 <- aggregate(overlap_vs_balanced ~ N_1st, t1b, median); add("retention_curves_max_abs_diff", max(abs(m1$overlap_vs_full - m2$overlap_vs_balanced)))
v <- q(P("data/pipeline/validation_batch_in_limma/test1_first_trim_subsample.tsv")); v <- v[v$status == "ok", ]
for (N1 in c(10, 40)) add(paste0("onestep_ccc_N", N1), median(v$logfc_ccc_vs_full[v$N_1st == N1]))
## ComBat sensitivity
a <- as.matrix(read.delim(P("data/pipeline/main/exprs_none_combat_ref.tsv"), row.names = 1, check.names = FALSE)); bb <- as.matrix(read.delim(P("data/pipeline/main/exprs_softimpute_combat_ref.tsv"), row.names = 1, check.names = FALSE))
g <- intersect(rownames(a), rownames(bb)); a <- a[g, colnames(bb)]; bb <- bb[g, ]; dd <- abs(a - bb); pg <- rowMeans(dd)
add("combat_cells", length(a)); add("combat_r", cor(as.vector(a), as.vector(bb))); add("combat_mean_abs_diff", mean(dd)); add("combat_median_gene_mae", median(pg))
ids <- c("3046", "6696", "3048", "885", "2162"); add("combat_named_genes_min_mae", min(pg[ids])); add("combat_named_genes_max_mae", max(pg[ids]))
## Lykhenko
ly <- read.csv(P("data/references/lykhenko_2021_full_protein_coding.csv")); m <- merge(d[d$gene %in% gained, ], ly, by.x = "gene", by.y = "ENTREZID"); same <- sign(m$logFC.x) == sign(m$logFC.y); s21 <- m$adj.P.Val.y < 0.05
add("ly_matched", nrow(m)); add("ly_matched_pct", 100 * nrow(m) / length(gained)); add("ly_absent", length(gained) - nrow(m)); add("ly_same_dir", sum(same)); add("ly_same_dir_pct", 100 * mean(same))
add("ly_r", cor(m$logFC.x, m$logFC.y)); add("ly_ccc", ccc(m$logFC.x, m$logFC.y)); add("ly_mean_abs_logfc_ours", mean(abs(m$logFC.x))); add("ly_mean_abs_logfc_2021", mean(abs(m$logFC.y))); add("ly_fdr_sig", sum(s21)); add("ly_fdr_sig_pct", 100 * mean(s21)); add("ly_discordant", sum(!same)); add("ly_discordant_pct", 100 * mean(!same))
add("ly_discordant_fdr_sig", sum(s21 & !same)); add("ly_discordant_median_abs_logfc", median(abs(m$logFC.y[!same]))); add("ly_concordant_median_abs_logfc", median(abs(m$logFC.y[same]))); add("ly_sig21_agreeing", sum(s21 & same)); add("ly_sig21", sum(s21))
## Prater
pr <- read.xlsx(loadWorkbook(P("data/references/prater_2021_supp_tables.xlsx")), sheet = 1); pr$entrez <- as.character(pr$entrezgene_id); pr <- pr[!is.na(pr$entrez) & pr$entrez != "NA", ]
plf <- setNames(as.numeric(pr$log2FoldChange), pr$entrez); ppa <- setNames(as.numeric(pr$padj), pr$entrez); psig <- pr$entrez[!is.na(ppa) & ppa < 0.05 & abs(plf) > 1]
for (lab in c("twostep", "onestep")) { dd2 <- if (lab == "twostep") d else o; SS <- if (lab == "twostep") S else O
  shp <- intersect(dd2$gene, names(plf)); x <- dd2$logFC[match(shp, dd2$gene)]; y <- plf[shp]; a1 <- shp %in% psig; b1 <- shp %in% SS
  add(paste0("prater_", lab, "_shared"), length(shp)); add(paste0("prater_", lab, "_r"), cor(x, y)); add(paste0("prater_", lab, "_ccc"), ccc(x, y))
  add(paste0("prater_", lab, "_lfc_subset_n"), sum(a1)); add(paste0("prater_", lab, "_lfc_subset_r"), cor(x[a1], y[a1])); add(paste0("prater_", lab, "_oursig_subset_n"), sum(b1)); add(paste0("prater_", lab, "_oursig_subset_r"), cor(x[b1], y[b1]))
  add(paste0("prater_", lab, "_replicated"), sum(SS %in% psig)); add(paste0("prater_", lab, "_replicated_pct"), 100 * mean(SS %in% psig)); add(paste0("prater_", lab, "_below_thr"), sum(SS %in% shp & !(SS %in% psig))); add(paste0("prater_", lab, "_absent"), sum(!(SS %in% shp)))
  add(paste0("prater_", lab, "_same_dir_pct_shared"), 100 * mean(sign(x) == sign(y))) }
gm <- read.csv(P("data/pipeline/ga_matched_prater/prater_concordance.csv")); for (nm in names(gm)) add(paste0("ga_matched_", nm), gm[[nm]][1])
## covariate sweep
cs <- read.csv(P("data/pipeline/covariate_sweep/table_s2_covariate_comparison.csv")); for (i in seq_len(nrow(cs))) for (nm in setdiff(names(cs), c("variant", "label"))) add(paste0("cov_", cs$variant[i], "_", nm), cs[[nm]][i])
js <- read.csv(P("data/pipeline/covariate_sweep/table_s3_covariate_stability.csv")); for (i in seq_len(nrow(js))) for (nm in setdiff(names(js), "N_1st")) add(paste0("covjac_N", js$N_1st[i], "_", nm), js[[nm]][i])
## enrichment
rd <- function(f) read.csv(P("data/pipeline/main/enrichment", f), stringsAsFactors = FALSE)
for (set in c("full_530", "intersection_277", "gained_253")) { gg <- rd(paste0("enrichment_", set, "_GO_BP.csv")); kk <- rd(paste0("enrichment_", set, "_KEGG.csv")); add(paste0("enr_", set, "_go"), nrow(gg)); add(paste0("enr_", set, "_kegg"), nrow(kk))
  for (id in c("GO:0001819", "GO:0050900", "GO:0006935", "GO:0006959", "GO:0002443")) if (id %in% gg$ID) add(paste0("enr_", set, "_", id), gg$qvalue[gg$ID == id])
  for (id in c("hsa05150", "hsa04610", "hsa04514", "hsa04145")) if (id %in% kk$ID) add(paste0("enr_", set, "_", id), kk$qvalue[kk$ID == id]) }
gg <- rd("enrichment_gained_253_GO_BP.csv"); ig <- rd("enrichment_intersection_277_GO_BP.csv"); gk <- rd("enrichment_gained_253_KEGG.csv"); ik <- rd("enrichment_intersection_277_KEGG.csv")
add("enr_gained_go_shared_with_intersection", length(intersect(gg$ID, ig$ID))); add("enr_gained_kegg_shared_with_intersection", length(intersect(gk$ID, ik$ID)))
out <- do.call(rbind, rows); write.csv(out, file.path(OUT_DIR, "de_sweep_numbers.csv"), row.names = FALSE); cat(nrow(out), "numbers written\n")

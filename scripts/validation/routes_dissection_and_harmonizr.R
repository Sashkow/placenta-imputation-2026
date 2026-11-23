#!/usr/bin/env Rscript
## Routes without imputation, 1 of 4: matrix dissection + ComBat (the idea of
## HarmonizR) against the article's softImpute + reference-batch ComBat, on the
## article's six datasets. The dissection here is OUR OWN code, written so that
## ComBat can be given the article's model and reference batch; the published
## HarmonizR package (arm hz_pkg) accepts neither.
##
## Run from the repository root; the other routes_*.R scripts add their arms to
## the tables this one writes:
##   Rscript scripts/validation/routes_dissection_and_harmonizr.R
##   Rscript scripts/validation/routes_batch_in_limma.R
##   Rscript scripts/validation/routes_bert.R
##   Rscript scripts/validation/routes_rnaseq_concordance.R --config=config/main.yaml --label=main
## Output: data/pipeline/routes/ (or --out=<dir>).
##
## Part A  full data: DE of each arm against the published softImpute run.
## Part B  holdout:   the stored block masks (data/pipeline/holdout/main/masks); DE calls of the masked
##                    genes against the full-data run, for each arm.
##
## Arms
##   hz_mod_ref  dissection by dataset membership; ComBat per sub-matrix with the
##               article's model matrix (gestational age + fetal sex) and
##               GSE100051 as reference when it is in the sub-matrix. If the
##               model is confounded with the datasets, ComBat is run without it.
##   hz_mod      the same without a reference batch (grand mean)
##   hz_pkg      the published HarmonizR package (ComBat mode 1, no covariates)
## All arms: limma on the measured cells only (NA cells are left NA).

HOLDOUT_DEFAULT_OUT <- "data/pipeline/routes"
source("scripts/validation/_holdout_common.R")
HOLDOUT_DIR <- "data/pipeline/holdout/main"

D <- prepare_data(config)
tools <- make_de_tools(D, config)
ph <- config$phenotype
covariates <- ph$covariates %||% character(0)

## ------------------------------------------------------------ dissection
## status per gene: "mod" (ComBat with the biological model), "nomod" (model
## confounded with batch, ComBat without it), "single" (one dataset, nothing to
## correct), "failed" (ComBat error, values left as they were)
hz_dissect <- function(X, use_ref = TRUE, use_mod = TRUE) {
  pres <- sapply(D$datasets, function(ds) rowSums(!is.na(X[, D$ds_cols[[ds]], drop = FALSE])) >= 2)
  key <- apply(pres, 1, function(r) paste(D$datasets[r], collapse = "|"))
  out <- X
  status <- setNames(rep("single", nrow(X)), rownames(X))
  groups <- split(rownames(X), key)
  info <- list()
  for (k in names(groups)) {
    g <- groups[[k]]
    dsets <- strsplit(k, "|", fixed = TRUE)[[1]]
    if (length(dsets) < 2) { status[g] <- if (length(dsets)) "single" else "empty"; next }
    cols <- unlist(D$ds_cols[dsets])
    sub <- X[g, cols, drop = FALSE]
    ok <- complete.cases(sub) & apply(sub, 1, var) > 0
    status[g[!ok]] <- "failed"
    g <- g[ok]; sub <- sub[ok, , drop = FALSE]
    if (nrow(sub) < 2) { status[g] <- "failed"; next }
    pd_sub <- D$pd[colnames(sub), , drop = FALSE]
    batch <- factor(pd_sub[[DATASET_ID]])
    ref <- if (use_ref && REF_BATCH %in% dsets) REF_BATCH else NULL
    res <- NULL; st <- "failed"
    if (use_mod) {
      mod <- combat_model_matrix(pd_sub, config, informative_covariates(pd_sub, covariates))
      res <- tryCatch(suppressMessages(ComBat(sub, batch, mod = mod, par.prior = TRUE,
                                              prior.plots = FALSE, ref.batch = ref)),
                      error = function(e) NULL)
      if (!is.null(res)) st <- "mod"
    }
    if (is.null(res)) {
      res <- tryCatch(suppressMessages(ComBat(sub, batch, mod = NULL, par.prior = TRUE,
                                              prior.plots = FALSE, ref.batch = ref)),
                      error = function(e) NULL)
      if (!is.null(res)) st <- "nomod"
    }
    if (!is.null(res) && all(is.finite(res))) out[g, cols] <- res else st <- "failed"
    status[g] <- st
    tr <- table(pd_sub[[DATASET_ID]], pd_sub[[ph$group_column]])
    info[[k]] <- data.frame(pattern = k, n_datasets = length(dsets), n_genes = length(g),
                            status = st, ref_in_pattern = REF_BATCH %in% dsets,
                            datasets_with_both_groups = sum(rowSums(tr > 0) >= 2),
                            stringsAsFactors = FALSE)
  }
  list(matrix = out, status = status, groups = do.call(rbind, info))
}

hz_package <- function(X) {
  desc <- data.frame(ID = colnames(X), sample = seq_len(ncol(X)),
                     batch = as.integer(factor(D$ds_of_sample[colnames(X)], levels = D$datasets)),
                     stringsAsFactors = FALSE)
  res <- as.matrix(HarmonizR::harmonizR(data_as_input = as.data.frame(X), description_as_input = desc,
                                        algorithm = "ComBat", ComBat_mode = 1, output_file = FALSE,
                                        verbosity = 0, cores = 1))
  out <- matrix(NA_real_, nrow(X), ncol(X), dimnames = dimnames(X))
  g <- intersect(rownames(res), rownames(X))
  out[g, colnames(res)] <- res[g, ]
  list(matrix = out, status = setNames(ifelse(rownames(X) %in% g, "pkg", "pkg_dropped"), rownames(X)),
       groups = NULL)
}

## limma on measured cells only (NA stays NA), the article's design
de_na <- function(mat) {
  E <- mat[, tools$keep_samples, drop = FALSE]
  pd_de <- D$pd[colnames(E), , drop = FALSE]
  E <- E[rowSums(!is.na(E)) >= 3, , drop = FALSE]
  de <- suppressWarnings(fit_limma(E, de_design(pd_de, config, informative_covariates(pd_de, covariates))))
  data.frame(gene = de$gene, logFC = de$logFC, adj.P.Val = de$adj.P.Val, stringsAsFactors = FALSE)
}

ARMS <- list(hz_mod_ref = function(X) hz_dissect(X, TRUE, TRUE),
             hz_mod     = function(X) hz_dissect(X, FALSE, TRUE),
             hz_pkg     = hz_package)

## ---------------------------------------------------------------- Part A
ref_de <- read.delim("data/pipeline/main/difexp_softimpute_combat_ref.tsv", stringsAsFactors = FALSE)
ref_de$sig <- sig_call(ref_de$logFC, ref_de$adj.P.Val)
none_de <- read.delim("data/pipeline/main/difexp_none_combat_ref.tsv", stringsAsFactors = FALSE)
log_msg("Reference softImpute + ComBat-ref: ", nrow(ref_de), " genes, ", sum(ref_de$sig), " DEGs")

full <- list(); summ_A <- list(); pergene <- data.frame(gene = D$genes, k = D$n_ds_per_gene,
                                                        ref_visible = D$presence[, REF_BATCH],
                                                        stringsAsFactors = FALSE)
pergene$logFC_softimpute <- ref_de$logFC[match(pergene$gene, ref_de$gene)]
pergene$padj_softimpute  <- ref_de$adj.P.Val[match(pergene$gene, ref_de$gene)]
pergene$sig_softimpute   <- ref_de$sig[match(pergene$gene, ref_de$gene)]
for (a in names(ARMS)) {
  t0 <- Sys.time()
  r <- ARMS[[a]](D$X)
  de <- de_na(r$matrix)
  i <- match(pergene$gene, de$gene)
  pergene[[paste0("logFC_", a)]]  <- de$logFC[i]
  pergene[[paste0("padj_", a)]]   <- de$adj.P.Val[i]
  pergene[[paste0("sig_", a)]]    <- sig_call(de$logFC[i], de$adj.P.Val[i])
  pergene[[paste0("status_", a)]] <- r$status[pergene$gene]
  if (!is.null(r$groups)) write.csv(r$groups, file.path(OUT_DIR, paste0("groups_", a, ".csv")), row.names = FALSE)
  saveRDS(r$matrix, file.path(OUT_DIR, paste0("corrected_", a, ".rds")))
  log_msg(a, sprintf(": %d genes tested, %d DEGs (%.0fs)", sum(!is.na(de$logFC)),
                     sum(pergene[[paste0("sig_", a)]], na.rm = TRUE),
                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
write.csv(pergene, file.path(OUT_DIR, "partA_per_gene.csv"), row.names = FALSE)

strata <- list(all = rep(TRUE, nrow(pergene)), shared_all6 = pergene$k == 6,
               partial_ref_visible = pergene$k < 6 & pergene$ref_visible,
               partial_ref_hidden = pergene$k < 6 & !pergene$ref_visible)
rows <- list()
for (a in names(ARMS)) for (s in names(strata)) {
  p <- pergene[strata[[s]], ]
  sa <- p[[paste0("sig_", a)]] %in% TRUE; sr <- p$sig_softimpute %in% TRUE
  ok <- is.finite(p[[paste0("logFC_", a)]]) & is.finite(p$logFC_softimpute)
  rows[[length(rows) + 1]] <- data.frame(
    arm = a, stratum = s, genes = nrow(p), tested = sum(!is.na(p[[paste0("logFC_", a)]])),
    deg_arm = sum(sa), deg_softimpute = sum(sr), both = sum(sa & sr),
    only_softimpute = sum(sr & !sa), only_arm = sum(sa & !sr),
    jaccard = sum(sa & sr) / max(1, sum(sa | sr)),
    logFC_r = cor(p[[paste0("logFC_", a)]][ok], p$logFC_softimpute[ok]),
    logFC_slope = unname(coef(lm(p[[paste0("logFC_", a)]][ok] ~ p$logFC_softimpute[ok]))[2]),
    opposite_sign_among_both = sum(sa & sr & sign(p[[paste0("logFC_", a)]]) != sign(p$logFC_softimpute)),
    stringsAsFactors = FALSE)
}
partA <- do.call(rbind, rows)
write.csv(partA, file.path(OUT_DIR, "partA_summary.csv"), row.names = FALSE)
print(partA, digits = 3)

st <- as.data.frame(table(status = pergene$status_hz_mod_ref, k = pergene$k,
                          ref_visible = pergene$ref_visible), stringsAsFactors = FALSE)
st <- st[st$Freq > 0, ]
write.csv(st, file.path(OUT_DIR, "partA_status_census.csv"), row.names = FALSE)
print(st)

## residual dataset effect: share of variance explained by dataset, per gene,
## among first-trimester samples of the genes all six datasets measure
ft <- D$grp_of_sample == ph$baseline
r2_dataset <- function(mat) {
  g <- pergene$gene[pergene$k == 6]
  E <- mat[g, ft, drop = FALSE]; b <- factor(D$ds_of_sample[ft])
  apply(E, 1, function(y) { o <- !is.na(y); if (sum(o) < 10) return(NA)
    summary(lm(y[o] ~ b[o]))$r.squared })
}
soft_corr <- as.matrix(read.delim("data/pipeline/main/exprs_softimpute_combat_ref.tsv",
                                  row.names = 1, check.names = FALSE))
resid <- data.frame(arm = c("uncorrected", "softimpute_combat_ref", names(ARMS)),
                    median_r2_dataset = c(median(r2_dataset(D$X), na.rm = TRUE),
                                          median(r2_dataset(soft_corr[, colnames(D$X)]), na.rm = TRUE),
                                          sapply(names(ARMS), function(a) median(r2_dataset(
                                            readRDS(file.path(OUT_DIR, paste0("corrected_", a, ".rds")))), na.rm = TRUE))))
write.csv(resid, file.path(OUT_DIR, "partA_residual_dataset_effect.csv"), row.names = FALSE)
print(resid, digits = 3)

## ---------------------------------------------------------------- Part B
ref_F <- read.csv(file.path(HOLDOUT_DIR, "reference", "softimpute_difexp.csv"), stringsAsFactors = FALSE)
rowsB <- list()
for (scheme in c("gene_dataset_block", "progressive_tax_block")) for (rep in 1:3) {
  tag <- sprintf("%s_rep%d", scheme, rep)
  mk <- readRDS(file.path(HOLDOUT_DIR, "masks", paste0(tag, ".rds")))
  stopifnot(mk$n_genes == nrow(D$X), mk$n_samples == ncol(D$X))
  vis <- mk$gene_visibility
  Xm <- D$X; Xm[mk$mask_idx] <- NA_real_
  keep_g <- !vis$lost_group; scored <- vis$n_masked_cells > 0
  soft <- read.csv(file.path(HOLDOUT_DIR, "de", paste0(tag, "_softimpute.csv")), stringsAsFactors = FALSE)
  calls <- data.frame(gene = D$genes[scored], ref_visible = vis$ref_batch_visible[scored],
                      k_visible = vis$k_visible[scored], lost_group = vis$lost_group[scored],
                      stringsAsFactors = FALSE)
  fi <- match(calls$gene, ref_F$gene)
  calls$logFC_full <- ref_F$logFC[fi]; calls$sig_full <- sig_call(ref_F$logFC[fi], ref_F$adj.P.Val[fi])
  si <- match(calls$gene, soft$gene)
  calls$logFC_softimpute <- soft$logFC_masked[si]; calls$sig_softimpute <- soft$sig_masked[si] %in% TRUE
  for (a in names(ARMS)) {
    r <- ARMS[[a]](Xm[keep_g, , drop = FALSE])
    de <- de_na(r$matrix)
    i <- match(calls$gene, de$gene)
    calls[[paste0("logFC_", a)]] <- de$logFC[i]
    calls[[paste0("sig_", a)]]   <- sig_call(de$logFC[i], de$adj.P.Val[i])
    calls[[paste0("status_", a)]] <- r$status[calls$gene]
  }
  write.csv(cbind(scheme = scheme, rep = rep, calls), file.path(OUT_DIR, paste0("partB_calls_", tag, ".csv")),
            row.names = FALSE)
  for (a in c("softimpute", names(ARMS))) for (s in c("overall", "ref_visible", "ref_hidden")) {
    p <- calls[switch(s, overall = TRUE, ref_visible = calls$ref_visible, ref_hidden = !calls$ref_visible), ]
    sa <- p[[paste0("sig_", a)]] %in% TRUE; tr <- p$sig_full %in% TRUE
    ok <- is.finite(p[[paste0("logFC_", a)]]) & is.finite(p$logFC_full)
    rowsB[[length(rowsB) + 1]] <- data.frame(
      scheme = scheme, rep = rep, arm = a, stratum = s, genes = nrow(p), true_deg = sum(tr),
      recovered = sum(sa & tr), lost = sum(tr & !sa), false_new = sum(sa & !tr),
      logFC_ccc = if (sum(ok) > 2) lin_ccc(p[[paste0("logFC_", a)]][ok], p$logFC_full[ok]) else NA,
      logFC_slope = if (sum(ok) > 2) unname(coef(lm(p[[paste0("logFC_", a)]][ok] ~ p$logFC_full[ok]))[2]) else NA,
      stringsAsFactors = FALSE)
  }
  log_msg(tag, " done")
}
partB_rep <- do.call(rbind, rowsB)
write.csv(partB_rep, file.path(OUT_DIR, "partB_by_repeat.csv"), row.names = FALSE)
num <- c("genes", "true_deg", "recovered", "lost", "false_new", "logFC_ccc", "logFC_slope")
partB <- aggregate(partB_rep[, num], partB_rep[, c("arm", "stratum", "scheme")], mean)
partB$lost_pct <- 100 * partB$lost / pmax(partB$true_deg, 1e-9)
partB <- partB[order(partB$scheme, partB$stratum, partB$arm), ]
write.csv(partB, file.path(OUT_DIR, "partB_summary.csv"), row.names = FALSE)
print(partB, digits = 3, row.names = FALSE)
cat("\nDone. Outputs in", OUT_DIR, "\n")

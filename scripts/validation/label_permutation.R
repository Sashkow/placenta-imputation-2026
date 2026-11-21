#!/usr/bin/env Rscript
# Label-permutation null: shuffle the trimester labels (together with each
# sample's gestational age, which enters ComBat as the biological covariate)
# across all 117 samples, keeping 102 first / 15 second trimester, then re-run
# batch correction + DE. With no real trimester signal, the pipeline should call
# (almost) no DEGs; the count of the real run is compared with this null.
#
# Fetal sex stays attached to its sample. Datasets, the gene set and the
# reused imputed matrix are unchanged.
#
# Run from repo root:
#   Rscript scripts/validation/label_permutation.R [--n=200] [--cores=10]
#   (article run: 200 permutations)

args <- commandArgs(trailingOnly = TRUE)
arg_int <- function(name, default) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) as.integer(sub(".*=", "", a)) else default
}
N_PERM <- arg_int("n", 200)
CORES  <- arg_int("cores", 10)

source("scripts/validation/_lodo_permutation_common.R")

check_full_run()

pd0 <- DATA$phenodata
ref_sig <- DATA$ref_sig_genes
one_perm <- function(b) {
  set.seed(20261113 + b)
  idx <- sample.int(nrow(pd0))
  pd <- pd0
  pd[[grp_col]]  <- pd0[[grp_col]][idx]
  pd$ga_weeks    <- pd0$ga_weeks[idx]
  r <- run_on(DATA$exprs_list, pd)
  if (r$status != "ok")
    return(data.frame(perm = b, status = r$status, n_deg = NA, n_fdr_only = NA,
                      n_p001 = NA, min_fdr = NA, overlap_ref = NA,
                      second_from_GSE100051 = NA, second_from_GSE9984 = NA))
  de <- r$de_results
  sec <- pd[[ds_col]][pd[[grp_col]] == config$phenotype$contrast]
  data.frame(perm = b, status = "ok",
             n_deg = r$n_sig, n_fdr_only = r$n_sig_fdr_only,
             n_p001 = sum(de$P.Value < 0.001, na.rm = TRUE),
             min_fdr = min(de$adj.P.Val, na.rm = TRUE),
             overlap_ref = length(intersect(r$sig_genes, ref_sig)),
             second_from_GSE100051 = sum(sec == "GSE100051"),
             second_from_GSE9984 = sum(sec == "GSE9984"))
}

t0 <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(N_PERM), one_perm, mc.cores = CORES))
cat(sprintf("%d permutations in %.1f min\n", N_PERM,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.table(res, file.path(OUT_DIR, "permutation_null.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

ok <- res[res$status == "ok", ]
n_obs <- length(ref_sig)
q <- function(x) paste(sprintf("%g", quantile(x, c(0, .5, .95, 1))), collapse = " / ")
cat(sprintf("ok runs: %d of %d\n", nrow(ok), N_PERM))
cat("DEGs (FDR<0.05 & |logFC|>1)  min / median / 95%% / max:", q(ok$n_deg), "\n")
cat("FDR<0.05 only                min / median / 95%% / max:", q(ok$n_fdr_only), "\n")
cat("p<0.001 genes (null exp. ~", round(0.001 * nrow(DATA$ref_de)), ") min / median / 95%% / max:", q(ok$n_p001), "\n")
cat(sprintf("runs with >=1 DEG: %d (%.1f%%); runs with >= %d DEGs: %d; empirical p = %.4f\n",
            sum(ok$n_deg > 0), 100 * mean(ok$n_deg > 0), n_obs,
            sum(ok$n_deg >= n_obs), (1 + sum(ok$n_deg >= n_obs)) / (nrow(ok) + 1)))
cat("Saved:", file.path(OUT_DIR, "permutation_null.tsv"), "\n")

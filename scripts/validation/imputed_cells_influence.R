#!/usr/bin/env Rscript
##
## How much do the imputed cells change the result on partially measured genes?
##
## Imputed cells have weight 0 in limma, so they can act only through ComBat's
## parameters. This script measures that: for every imputer (stored full-data
## reference runs of the holdout, data/pipeline/holdout/main/reference) and for softImpute fitted with
## other RNG seeds, it compares with the softImpute run of the article
##   - the ComBat-corrected values of the MEASURED cells, and
##   - the DEG calls and log fold changes,
## separately for genes measured in all datasets, partially measured genes the
## reference dataset measures, and genes it does not measure.
##
## Run from the repository root:
##   Rscript scripts/validation/imputed_cells_influence.R [--seeds=1,2,3,4,5]
## Output: data/pipeline/imputed_influence/ (or --out=<dir>).

HOLDOUT_DEFAULT_OUT <- "data/pipeline/imputed_influence"
source("scripts/validation/_holdout_common.R")

REF_DIR <- file.path(config$paths$output, "reference")   # stored full-data runs
SEEDS <- as.integer(split_arg("seeds", c("1", "2", "3", "4", "5")))
dir.create(file.path(OUT_DIR, "seeds"), recursive = TRUE, showWarnings = FALSE)

D <- prepare_data(config)
tools <- make_de_tools(D, config)

n_missing_ds <- length(D$datasets) - D$n_ds_per_gene
stratum <- ifelse(n_missing_ds == 0, "all_datasets",
           ifelse(D$presence[, REF_BATCH], "partial_ref_visible", "partial_ref_hidden"))
names(stratum) <- D$genes
cat("Genes per stratum:\n"); print(table(stratum))

base_de  <- read.csv(file.path(REF_DIR, "softimpute_difexp.csv"), stringsAsFactors = FALSE)
base_cor <- readRDS(file.path(REF_DIR, "softimpute_corrected.rds"))[D$genes, D$samples]
measured <- !D$orig_mask

## --------------------------------------------- softImpute with other seeds
seed_runs <- list()
for (s in SEEDS) {
  f_de <- file.path(OUT_DIR, "seeds", sprintf("softimpute_seed%d_difexp.csv", s))
  f_co <- file.path(OUT_DIR, "seeds", sprintf("softimpute_seed%d_corrected.rds", s))
  if (!file.exists(f_de) || !file.exists(f_co) || FORCE) {
    log_msg("softImpute seed ", s)
    cfg <- method_cfg("softimpute"); cfg$seed <- s
    imp <- run_imputer("softimpute", D$incomplete, cfg)
    corrected <- tools$combat_correct(imp$matrix)
    de <- tools$de_fit(corrected, D$orig_mask)
    de$sig <- sig_call(de$logFC, de$adj.P.Val)
    saveRDS(corrected, f_co); write.csv(de, f_de, row.names = FALSE)
  }
  seed_runs[[sprintf("softimpute_seed%d", s)]] <- list(de = f_de, corrected = f_co)
}

## ------------------------------------------------------------- comparison
stored <- setdiff(names(METHOD_LABELS), "softimpute")
stored <- stored[file.exists(file.path(REF_DIR, paste0(stored, "_difexp.csv")))]
runs <- c(setNames(lapply(stored, function(m) list(
            de = file.path(REF_DIR, paste0(m, "_difexp.csv")),
            corrected = file.path(REF_DIR, paste0(m, "_corrected.rds")))), stored),
          seed_runs)

compare <- function(name, r) {
  de  <- read.csv(r$de, stringsAsFactors = FALSE)
  co  <- readRDS(r$corrected)[D$genes, D$samples]
  de  <- de[match(base_de$gene, de$gene), ]
  do.call(rbind, lapply(c("all_datasets", "partial_ref_visible", "partial_ref_hidden"), function(st) {
    g  <- D$genes[stratum == st]
    m  <- measured[g, , drop = FALSE]
    x  <- base_cor[g, , drop = FALSE][m]; y <- co[g, , drop = FALSE][m]
    d  <- abs(x - y)
    ig <- base_de$gene %in% g
    b_sig <- base_de$gene[ig & base_de$sig]; m_sig <- de$gene[ig & de$sig]
    ok <- ig & is.finite(base_de$logFC) & is.finite(de$logFC)
    data.frame(run = name, stratum = st, n_genes = length(g), n_measured_cells = length(x),
               value_r = cor(x, y), value_mean_abs_diff = mean(d),
               value_p99_abs_diff = unname(quantile(d, 0.99)), value_max_abs_diff = max(d),
               logfc_r = cor(base_de$logFC[ok], de$logFC[ok]),
               logfc_max_abs_diff = max(abs(base_de$logFC[ok] - de$logFC[ok])),
               deg_softimpute = length(b_sig), deg_run = length(m_sig),
               deg_shared = length(intersect(b_sig, m_sig)),
               deg_jaccard = if (length(union(b_sig, m_sig)))
                 length(intersect(b_sig, m_sig)) / length(union(b_sig, m_sig)) else NA,
               stringsAsFactors = FALSE)
  }))
}
res <- do.call(rbind, Map(compare, names(runs), runs))
res$label <- ifelse(res$run %in% names(METHOD_LABELS), METHOD_LABELS[res$run], res$run)
write.csv(res, file.path(OUT_DIR, "imputed_influence.csv"), row.names = FALSE)
options(width = 200)
print(res[, c("label", "stratum", "n_genes", "value_r", "value_mean_abs_diff", "value_max_abs_diff",
              "logfc_r", "deg_softimpute", "deg_run", "deg_shared", "deg_jaccard")],
      row.names = FALSE, digits = 4)
cat("\nWritten:", file.path(OUT_DIR, "imputed_influence.csv"), "\n")

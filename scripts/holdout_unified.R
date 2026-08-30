#!/usr/bin/env Rscript
##
## Unified holdout validation: one mask, two readouts.
##
## For each masking scheme x repeat, ONE mask of `leave_out_fraction` of the
## observed cells is drawn with the companion's own mask functions and seeds
## (set.seed(rep * 123), as validate_imputation() does) and saved. Then, for
## every imputer:
##
##   values readout  the hidden cells are imputed and compared with the truth
##                   (pooled r / RMSE / MAE, and per-gene MAE so that accuracy
##                   can be stratified by what the gene was left with);
##
##   calls readout   the imputed matrix continues through ComBat-ref and limma
##                   with weight 0 on masked + originally-missing cells (arm M,
##                   the production treatment) and its DE calls are compared
##                   with the same imputer's full-data reference run (arm F).
##                   For the block schemes an oracle arm O -- limma on F's
##                   corrected matrix with M's weights -- separates the ComBat
##                   pathway from plain loss of samples.
##
## Every (scheme, repeat, imputer) iteration is checkpointed; rerunning skips
## finished iterations unless --force is given.
##
## Usage:
##   Rscript scripts/holdout_unified.R
##   Rscript scripts/holdout_unified.R --schemes=random_cells --methods=softimpute,knn \
##           --n_repeats=1 --frac=0.02 --out=output_smoke

source(file.path({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "."
}, "holdout_common.R"))

cat("=== Unified holdout validation ===\n")
cat("Companion:", COMPANION, "\nOutput:", OUT_DIR, "\n")
cat("Schemes:", paste(MASK_TYPES, collapse = ", "), "\nMethods:", paste(METHODS, collapse = ", "),
    "\nRepeats:", N_REPEATS, " Leave-out fraction:", LEAVE_OUT, "\n\n")

MIN_OBS   <- as.integer(HCFG$min_obs_per_gene %||% 4)
SEED_MULT <- as.integer(HCFG$seed_multiplier %||% 123)
MIN_CELLS <- as.integer(HCFG$min_cells_per_gene_metric %||% 3)
ORACLE_SCHEMES <- unlist(HCFG$oracle_schemes) %||% c("gene_dataset_block", "progressive_tax_block")

D <- prepare_data(config)
tools <- make_de_tools(D, config)
X <- D$X
## worst-case floor: the supplement's published definition uses n_min = 2; the
## matrix now also holds k = 1 genes, so the floor is a config setting
WORST_N_MIN <- as.integer(HCFG$worst_case_n_min %||% min(D$n_ds_per_gene))
observed_idx <- which(!is.na(X))
n_mask_target <- round(length(observed_idx) * LEAVE_OUT)
n_datasets_total <- length(D$datasets)
row_idx_map <- setNames(seq_len(nrow(X)), D$genes)
cat("Observed cells:", length(observed_idx), " target masked per repeat:", n_mask_target, "\n\n")

## reference runs are needed for the calls readout
refs <- lapply(setNames(METHODS, METHODS), function(m)
  tryCatch(load_reference(m), error = function(e) { log_msg(conditionMessage(e)); NULL }))
missing_ref <- names(refs)[sapply(refs, is.null)]
if (length(missing_ref))
  stop("Missing reference runs for: ", paste(missing_ref, collapse = ", "))

## --------------------------------------------------------------- masks
draw_mask <- function(scheme, rep) {
  set.seed(rep * SEED_MULT)
  if (scheme == "random_cells") {
    m <- mask_random_cells(observed_idx, n_mask_target)
  } else if (scheme == "gene_dataset_block") {
    ## (gene, dataset) pairs of genes seen in all datasets, as validate_imputation()
    eligible <- D$n_ds_per_gene == n_datasets_total
    obs_pairs <- do.call(rbind, lapply(D$datasets, function(ds) {
      g <- D$genes[D$presence[, ds] & eligible]
      if (length(g)) data.frame(gene = g, dataset = ds, stringsAsFactors = FALSE) else NULL
    }))
    rownames(obs_pairs) <- NULL
    m <- mask_gene_dataset_block(obs_pairs, X, D$ds_cols, row_idx_map, n_mask_target, MIN_OBS)
  } else if (scheme == "progressive_tax_block") {
    ## the supplement defines the worst case at the coverage floor n_min = 2;
    ## the current matrix also holds k = 1 genes, so the floor is configurable
    n_min <- WORST_N_MIN
    m <- mask_progressive_tax_block(X, D$presence, D$n_ds_per_gene, D$ds_cols, n_min, n_mask_target)
  } else stop("unknown scheme: ", scheme)
  m
}

get_mask <- function(scheme, rep) {
  f <- file.path(OUT_DIR, "masks", sprintf("%s_rep%d.rds", scheme, rep))
  if (file.exists(f) && !FORCE) {
    mk <- readRDS(f)
    if (identical(mk$leave_out_fraction, LEAVE_OUT) && identical(mk$n_genes, nrow(X)) &&
        identical(mk$worst_case_n_min %||% NA_integer_, WORST_N_MIN %||% NA_integer_))
      return(mk)
    log_msg("stored mask ", basename(f), " was drawn under different settings; redrawing")
  }
  m <- draw_mask(scheme, rep)
  new_mask <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
  new_mask[m$mask_idx] <- TRUE
  vis <- gene_visibility(D, new_mask, config)
  mk <- list(scheme = scheme, rep = rep, seed = rep * SEED_MULT,
             leave_out_fraction = LEAVE_OUT, n_genes = nrow(X), n_samples = ncol(X),
             worst_case_n_min = WORST_N_MIN,
             mask_idx = m$mask_idx, n_blocks = m$n_blocks,
             gene_visibility = vis)
  saveRDS(mk, f)
  mk
}

## --------------------------------------------------------- per-gene values
per_gene_values <- function(new_mask, filled, vis) {
  sel <- which(vis$n_masked_cells >= MIN_CELLS)
  if (!length(sel)) return(NULL)
  out <- lapply(sel, function(i) {
    idx <- new_mask[i, ]
    t <- X[i, idx]; f <- filled[i, idx]
    c(mae = mean(abs(t - f)), rmse = sqrt(mean((t - f)^2)),
      r = suppressWarnings(cor(t, f)), truth_sd = sd(t))
  })
  out <- do.call(rbind, out)
  data.frame(gene = vis$gene[sel], n_masked_cells = vis$n_masked_cells[sel],
             k_visible = vis$k_visible[sel], ref_batch_visible = vis$ref_batch_visible[sel],
             out, stringsAsFactors = FALSE, row.names = NULL)
}

## ------------------------------------------------------------------ loop
for (scheme in MASK_TYPES) {
  for (rep in seq_len(N_REPEATS)) {
    tag <- sprintf("%s_rep%d", scheme, rep)
    todo <- METHODS[!file.exists(file.path(OUT_DIR, "checkpoints", paste0(tag, "_", METHODS, ".done"))) | FORCE]
    if (!length(todo)) { log_msg(tag, ": all methods done, skipping"); next }

    cat("\n########## ", tag, " ##########\n", sep = "")
    mk <- get_mask(scheme, rep)
    mask_idx <- mk$mask_idx
    vis <- mk$gene_visibility
    new_mask <- matrix(FALSE, nrow(X), ncol(X), dimnames = dimnames(X))
    new_mask[mask_idx] <- TRUE
    mask_all <- D$orig_mask | new_mask
    truth <- X[mask_idx]
    log_msg(tag, ": ", length(mask_idx), " cells masked in ", sum(vis$n_masked_cells > 0),
            " genes (blocks: ", mk$n_blocks, "); total missing now ",
            sprintf("%.1f%%", 100 * mean(mask_all)),
            "; genes losing a whole group: ", sum(vis$lost_group))

    X_masked <- X; X_masked[mask_idx] <- NA_real_
    inc_masked <- D$incomplete; inc_masked$matrix <- X_masked

    ## genes the production pipeline would drop before ComBat/limma
    keep_g <- !vis$lost_group
    scored <- if (scheme == "random_cells") rep(TRUE, nrow(X)) else vis$n_masked_cells > 0

    for (m in todo) {
      it <- paste0(tag, "_", m)
      cat("\n--- ", it, " ---\n", sep = "")
      t0 <- Sys.time()
      imp <- tryCatch(run_imputer(m, inc_masked, method_cfg(m)),
                      error = function(e) { log_msg(it, ": imputation FAILED: ", conditionMessage(e)); NULL })
      t_imp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      ## ---- values readout
      if (is.null(imp)) {
        vrow <- data.frame(scheme = scheme, rep = rep, method = m, n_masked = length(mask_idx),
                           n_blocks = mk$n_blocks, pearson_r = NA, rmse = NA, mae = NA,
                           converged = FALSE, imputation_seconds = t_imp,
                           combat_status = "not_run", stringsAsFactors = FALSE)
        write.csv(vrow, file.path(OUT_DIR, "values", paste0(it, ".csv")), row.names = FALSE)
        file.create(file.path(OUT_DIR, "checkpoints", paste0(it, ".done")))
        next
      }
      filled <- imp$matrix[mask_idx]
      vrow <- data.frame(scheme = scheme, rep = rep, method = m, n_masked = length(mask_idx),
                         n_blocks = mk$n_blocks,
                         pearson_r = cor(truth, filled, use = "complete.obs"),
                         rmse = sqrt(mean((truth - filled)^2, na.rm = TRUE)),
                         mae = mean(abs(truth - filled), na.rm = TRUE),
                         converged = TRUE, imputation_seconds = t_imp,
                         combat_status = NA_character_, stringsAsFactors = FALSE)
      pg <- per_gene_values(new_mask, imp$matrix, vis)
      if (!is.null(pg)) {
        pg <- cbind(scheme = scheme, rep = rep, method = m, pg)
        write.csv(pg, file.path(OUT_DIR, "values", paste0(it, "_per_gene.csv")), row.names = FALSE)
      }
      log_msg(it, sprintf(": values r=%.3f RMSE=%.3f MAE=%.3f (%.0fs)",
                          vrow$pearson_r, vrow$rmse, vrow$mae, t_imp))

      ## ---- calls readout (arm M)
      corrected <- tryCatch(tools$combat_correct(imp$matrix[keep_g, , drop = FALSE]),
                            error = function(e) e)
      if (inherits(corrected, "error")) {
        vrow$combat_status <- paste0("error: ", conditionMessage(corrected))
        log_msg(it, ": ComBat FAILED: ", conditionMessage(corrected))
        write.csv(vrow, file.path(OUT_DIR, "values", paste0(it, ".csv")), row.names = FALSE)
        file.create(file.path(OUT_DIR, "checkpoints", paste0(it, ".done")))
        next
      }
      vrow$combat_status <- "ok"
      write.csv(vrow, file.path(OUT_DIR, "values", paste0(it, ".csv")), row.names = FALSE)

      de_M <- tools$de_fit(corrected, mask_all)
      de_F <- refs[[m]]$de
      de_O <- if (scheme %in% ORACLE_SCHEMES)
        tools$de_fit(refs[[m]]$corrected[keep_g, , drop = FALSE], mask_all) else NULL

      g <- D$genes[scored]
      fi <- match(g, de_F$gene); mi <- match(g, de_M$gene)
      rows <- data.frame(
        scheme = scheme, rep = rep, method = m, gene = g,
        n_masked_cells = vis$n_masked_cells[scored], k_visible = vis$k_visible[scored],
        ref_batch_visible = vis$ref_batch_visible[scored], lost_group = vis$lost_group[scored],
        logFC_full = de_F$logFC[fi], padj_full = de_F$adj.P.Val[fi],
        logFC_masked = de_M$logFC[mi], padj_masked = de_M$adj.P.Val[mi],
        stringsAsFactors = FALSE, row.names = NULL)
      if (!is.null(de_O)) {
        oi <- match(g, de_O$gene)
        rows$logFC_oracle <- de_O$logFC[oi]; rows$padj_oracle <- de_O$adj.P.Val[oi]
      } else {
        rows$logFC_oracle <- NA_real_; rows$padj_oracle <- NA_real_
      }
      rows$sig_full   <- sig_call(rows$logFC_full,   rows$padj_full)
      rows$sig_masked <- sig_call(rows$logFC_masked, rows$padj_masked)
      rows$sig_oracle <- sig_call(rows$logFC_oracle, rows$padj_oracle)
      rows$dropped_masked   <- rows$lost_group | is.na(rows$logFC_masked)
      rows$nonfinite_masked <- !rows$dropped_masked & !is.finite(rows$logFC_masked)
      rows$nonfinite_full   <- !is.finite(rows$logFC_full)
      write.csv(rows, file.path(OUT_DIR, "de", paste0(it, ".csv")), row.names = FALSE)

      tp <- rows$sig_full %in% TRUE
      log_msg(it, sprintf(": calls  scored=%d trueDEG=%d recovered=%d  spec=%.4f  dropped=%d nonfinite=%d",
                          nrow(rows), sum(tp), sum(rows$sig_masked[tp] %in% TRUE),
                          mean(!(rows$sig_masked[!tp] %in% TRUE)), sum(rows$dropped_masked),
                          sum(rows$nonfinite_masked)))
      file.create(file.path(OUT_DIR, "checkpoints", paste0(it, ".done")))
      rm(imp, corrected, de_M, de_O, rows); invisible(gc())
    }
  }
}

cat("\nAll requested iterations done. Outputs in", OUT_DIR, "\n")

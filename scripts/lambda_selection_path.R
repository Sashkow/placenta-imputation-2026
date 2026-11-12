#!/usr/bin/env Rscript
##
## Hold-out selection of the softImpute lambda (Mazumder, Hastie & Tibshirani
## 2010; one-SE rule of Hastie, Tibshirani & Friedman 2009, Sec. 7.10).
##
## The production warm-start path -- 8 geometric steps from lambda0 to
## 0.01 * lambda0, each fit started from the previous one (softimpute_path() in
## scripts/pipeline/imputation.R) -- is fitted on masked copies of the merged
## matrix, one mask per scheme and seed, and the hidden cells are scored at
## every step. The selected step is the largest-lambda step whose mean
## gene--dataset block RMSE over the seeds lies within one standard error of the
## lowest mean. That step is the `path_stop_step` of every pipeline config.
##
## The selection masks use their own seeds (config lambda_selection.seeds),
## distinct from the Table S1/S3 holdout masks; their overlap with those masks
## is reported in mask_overlap.csv.
##
## Outputs (data/pipeline/lambda_selection/):
##   masks/mask_<scheme>_<seed>.rds   hidden cells (linear indices into the merged matrix)
##   path_by_mask.csv                 per mask and step: lambda, rank, r, RMSE, MAE
##   path_full.csv                    the unmasked path: lambda, rank per step
##   path.csv                         per scheme and step: mean +/- SE over seeds
##   selected_step.txt                the chosen step and how it was chosen
##   mask_overlap.csv                 shared cells with the Table S1/S3 masks
##
## Usage:
##   Rscript scripts/lambda_selection_path.R [--config=config/config_lambda_selection.yaml] [--force]

HOLDOUT_DEFAULT_CONFIG <- "config/config_lambda_selection.yaml"
source(file.path({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "scripts"
}, "holdout_common.R"))
suppressPackageStartupMessages(library(parallel))

LCFG <- config$lambda_selection
SI   <- coerce_numeric(config$imputation$softimpute)
SEEDS      <- as.integer(unlist(LCFG$seeds))
SCHEMES    <- unlist(LCFG$mask_types)
SELECT_ON  <- LCFG$select_on %||% "gene_dataset_block"
LEAVE_OUT  <- as.numeric(LCFG$leave_out_fraction %||% 0.1)
MIN_OBS    <- as.integer(LCFG$min_obs_per_gene %||% 4)
N_STEPS    <- as.integer(LCFG$n_steps %||% 8)
MIN_RATIO  <- as.numeric(LCFG$min_ratio %||% 0.01)
CORES      <- as.integer(LCFG$cores %||% 1)

cat("=== softImpute lambda selection ===\nOutput:", OUT_DIR, "\nSeeds:", SEEDS,
    " schemes:", SCHEMES, " rank_max:", SI$rank_max, "\n\n")

D <- prepare_data(config)
X <- D$X
observed_idx <- which(!is.na(X))
n_mask_target <- round(length(observed_idx) * LEAVE_OUT)
row_idx_map <- setNames(seq_len(nrow(X)), D$genes)

## ------------------------------------------------------------------ masks
## same mask functions and eligibility as holdout_unified.R, own seeds
draw_mask <- function(scheme, seed) {
  set.seed(seed)
  if (scheme == "random_cells") {
    m <- mask_random_cells(observed_idx, n_mask_target)
  } else if (scheme == "gene_dataset_block") {
    eligible <- D$n_ds_per_gene == length(D$datasets)
    obs_pairs <- do.call(rbind, lapply(D$datasets, function(ds) {
      g <- D$genes[D$presence[, ds] & eligible]
      if (length(g)) data.frame(gene = g, dataset = ds, stringsAsFactors = FALSE) else NULL
    }))
    m <- mask_gene_dataset_block(obs_pairs, X, D$ds_cols, row_idx_map, n_mask_target, MIN_OBS)
  } else stop("unsupported scheme: ", scheme)
  m
}
get_mask <- function(scheme, seed) {
  f <- file.path(OUT_DIR, "masks", sprintf("mask_%s_%d.rds", scheme, seed))
  if (file.exists(f) && !FORCE) {
    mk <- readRDS(f)
    if (identical(mk$n_genes, nrow(X)) && identical(mk$leave_out_fraction, LEAVE_OUT)) return(mk)
  }
  m <- draw_mask(scheme, seed)
  mk <- list(scheme = scheme, seed = seed, leave_out_fraction = LEAVE_OUT,
             n_genes = nrow(X), n_samples = ncol(X),
             mask_idx = m$mask_idx, n_blocks = m$n_blocks)
  saveRDS(mk, f)
  mk
}
runs <- list(list(scheme = "full", seed = NA_integer_, mask_idx = integer(0)))
for (sc in SCHEMES) for (sd in SEEDS) runs[[length(runs) + 1]] <- get_mask(sc, sd)

## -------------------------------------------------------------- the paths
fit_run <- function(run) {
  Y <- X
  if (length(run$mask_idx)) Y[run$mask_idx] <- NA
  truth <- X[run$mask_idx]
  set.seed(SI$seed %||% 42)            # as impute_softimpute()
  Yc <- biScale(Y, maxit = 20, thresh = 1e-3)
  rows <- list()
  t0 <- Sys.time()
  path <- softimpute_path(Yc, rank_max = SI$rank_max, thresh = SI$thresh, maxit = SI$maxit,
                          type = SI$type %||% "als", n_steps = N_STEPS, min_ratio = MIN_RATIO,
    on_step = function(step, lam, fit) {
      sc <- c(r = NA_real_, rmse = NA_real_, mae = NA_real_)
      if (length(run$mask_idx)) {
        pred <- complete(Yc, fit)[run$mask_idx]
        sc <- c(r = cor(pred, truth), rmse = sqrt(mean((pred - truth)^2)),
                mae = mean(abs(pred - truth)))
      }
      rows[[length(rows) + 1]] <<- data.frame(
        scheme = run$scheme, seed = run$seed, step = step, lambda = lam,
        lambda_ratio = MIN_RATIO^(step / (N_STEPS - 1)), rank = softimpute_rank(fit),
        rank_cap_reached = softimpute_rank(fit) >= SI$rank_max,
        r = sc[["r"]], rmse = sc[["rmse"]], mae = sc[["mae"]],
        minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    })
  do.call(rbind, rows)
}
log_msg("fitting ", length(runs), " paths on ", CORES, " cores")
res <- mclapply(runs, function(run) tryCatch(fit_run(run), error = function(e) conditionMessage(e)),
                mc.cores = CORES, mc.preschedule = FALSE)
bad <- !sapply(res, is.data.frame)
if (any(bad)) stop("path failed: ", paste(unlist(res[bad]), collapse = "; "))
all_rows <- do.call(rbind, res)

full <- all_rows[all_rows$scheme == "full", c("step", "lambda", "lambda_ratio", "rank", "rank_cap_reached")]
by_mask <- all_rows[all_rows$scheme != "full", ]
write.csv(full, file.path(OUT_DIR, "path_full.csv"), row.names = FALSE)
write.csv(by_mask, file.path(OUT_DIR, "path_by_mask.csv"), row.names = FALSE)

## ---------------------------------------------------------- summary + rule
se <- function(x) sd(x) / sqrt(length(x))
summ <- do.call(rbind, lapply(split(by_mask, list(by_mask$scheme, by_mask$step), drop = TRUE), function(g)
  data.frame(scheme = g$scheme[1], step = g$step[1], lambda_ratio = g$lambda_ratio[1],
             lambda_full = full$lambda[full$step == g$step[1]],
             rank_full = full$rank[full$step == g$step[1]],
             rank_masked_mean = mean(g$rank), n_masks = nrow(g),
             rmse_mean = mean(g$rmse), rmse_se = se(g$rmse),
             r_mean = mean(g$r), r_se = se(g$r), mae_mean = mean(g$mae))))
summ <- summ[order(summ$scheme, summ$step), ]

sel <- summ[summ$scheme == SELECT_ON, ]
best <- sel[which.min(sel$rmse_mean), ]
limit <- best$rmse_mean + best$rmse_se
k_star <- min(sel$step[sel$rmse_mean <= limit])   # smallest step = largest lambda
summ$selected <- summ$scheme == SELECT_ON & summ$step == k_star
summ$min_error <- summ$scheme == SELECT_ON & summ$step == best$step
write.csv(summ, file.path(OUT_DIR, "path.csv"), row.names = FALSE)

msg <- c(sprintf("selected_step: %d", k_star),
         sprintf("rule: one-SE on %s RMSE over %d masks (seeds %s)", SELECT_ON, length(SEEDS),
                 paste(SEEDS, collapse = ", ")),
         sprintf("min_error_step: %d (mean RMSE %.4f, SE %.4f); one-SE limit %.4f",
                 best$step, best$rmse_mean, best$rmse_se, limit),
         sprintf("selected: lambda = %.4f = lambda0 * %.4f on the full matrix, rank %d (rank_max %d)",
                 full$lambda[full$step == k_star], MIN_RATIO^(k_star / (N_STEPS - 1)),
                 full$rank[full$step == k_star], as.integer(SI$rank_max)))
writeLines(msg, file.path(OUT_DIR, "selected_step.txt"))
cat("\n", paste(msg, collapse = "\n"), "\n\n", sep = "")
print(summ, row.names = FALSE, digits = 4)

## ------------------------------------------- overlap with Table S1/S3 masks
hm_dir <- file.path(ALT_DIR, config$paths$holdout_masks)
eligible_cells <- sum(!is.na(X[D$n_ds_per_gene == length(D$datasets), ]))
ov <- list()
for (run in runs[-1]) {
  for (f in list.files(hm_dir, pattern = paste0("^", run$scheme, "_rep[0-9]+\\.rds$"), full.names = TRUE)) {
    hm <- readRDS(f)
    if (!identical(hm$n_genes, nrow(X))) { log_msg("skip ", basename(f), ": different matrix"); next }
    pool <- if (run$scheme == "random_cells") length(observed_idx) else eligible_cells
    ov[[length(ov) + 1]] <- data.frame(
      scheme = run$scheme, selection_seed = run$seed, holdout_mask = basename(f),
      n_selection = length(run$mask_idx), n_holdout = length(hm$mask_idx),
      n_shared = length(intersect(run$mask_idx, hm$mask_idx)),
      shared_frac = length(intersect(run$mask_idx, hm$mask_idx)) / length(run$mask_idx),
      chance_frac = length(hm$mask_idx) / pool)
  }
}
if (length(ov)) {
  ov <- do.call(rbind, ov)
  write.csv(ov, file.path(OUT_DIR, "mask_overlap.csv"), row.names = FALSE)
  cat("\n=== overlap with the Table S1/S3 holdout masks ===\n"); print(ov, row.names = FALSE, digits = 3)
}

## holdout_common.R creates directories this script does not use
for (d in c("checkpoints", "reference", "values", "de")) {
  p <- file.path(OUT_DIR, d)
  if (dir.exists(p) && !length(list.files(p))) unlink(p, recursive = TRUE)
}
log_msg("done")

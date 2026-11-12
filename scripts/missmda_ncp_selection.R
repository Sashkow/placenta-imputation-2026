#!/usr/bin/env Rscript
##
## Choice of the imputePCA comparator's number of components (ncp).
##
## Two criteria, so the comparator is tuned as seriously as softImpute:
##   hold-out  imputePCA (missMDA, method "Regularized", package defaults
##             otherwise) on the gene--dataset block selection masks of
##             scripts/lambda_selection_path.R, over the ncp grid
##             lambda_selection.missmda_ncp, with the raw matrix as input (as
##             in the pipeline: genes centred and scaled) and with the
##             biScale()d matrix as input (gene AND sample offsets removed, as
##             softImpute gets it); chosen by the same one-SE rule as lambda.
##   GCV       missMDA's own criterion (estim_ncpPCA(method.cv = "gcv"),
##             Josse & Husson 2012) on the unmasked matrix for ncp 1..gcv_max,
##             computed with the package's formula, one fit per ncp in
##             parallel; the suggested ncp is the first local minimum, as
##             estim_ncpPCA returns it.
## imputePCA signals hitting its iteration cap only by a warning; it is
## recorded per fit (`capped`).
##
## Needs the selection masks: run scripts/lambda_selection_path.R first.
## Every fit is checkpointed under data/pipeline/lambda_selection/missmda/.
##
## Outputs (data/pipeline/lambda_selection/):
##   missmda_by_mask.csv  per mask, input and ncp: r, RMSE, MAE, minutes, capped
##   missmda_gcv.csv      GCV criterion per ncp
##   missmda_ncp.csv      per input and ncp: mean +/- SE block RMSE, selection flags
##   missmda_selected.txt chosen ncp and input, and the GCV suggestion
##
## Usage:
##   Rscript scripts/missmda_ncp_selection.R [--config=config/config_lambda_selection.yaml]

HOLDOUT_DEFAULT_CONFIG <- "config/config_lambda_selection.yaml"
source(file.path({
  fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else "scripts"
}, "holdout_common.R"))
suppressPackageStartupMessages({ library(parallel); library(missMDA) })

LCFG    <- config$lambda_selection
MM      <- config$imputation$missmda
SEEDS   <- as.integer(unlist(LCFG$seeds))
NCPS    <- as.integer(unlist(LCFG$missmda_ncp))
GCV_MAX <- as.integer(LCFG$missmda_gcv_max %||% 20)
CORES   <- as.integer(LCFG$missmda_cores %||% 1)
CK_DIR  <- file.path(OUT_DIR, "missmda")
dir.create(CK_DIR, showWarnings = FALSE)

D <- prepare_data(config)
X <- D$X

## imputePCA on genes x samples data (it wants samples in rows); warning-aware
fit_mm <- function(Y, ncp, scale) {
  capped <- FALSE
  t0 <- Sys.time()
  res <- withCallingHandlers(
    imputePCA(t(Y), ncp = ncp, scale = scale, method = MM$method %||% "Regularized"),
    warning = function(w) {
      if (grepl("Stopped after", conditionMessage(w))) capped <<- TRUE
      invokeRestart("muffleWarning")
    })
  list(completed = t(res$completeObs), fitted = t(res$fittedX), capped = capped,
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

## completion with raw or biScale()d input, on the original scale
complete_mm <- function(Y, ncp, input) {
  if (input == "raw") return(fit_mm(Y, ncp, isTRUE(MM$scale %||% TRUE)))
  Yc <- biScale(Y, maxit = 20, thresh = 1e-3)
  f <- fit_mm(unclass(Yc), ncp, FALSE)
  rw <- attr(Yc, "biScale:row"); cl <- attr(Yc, "biScale:column")
  f$completed <- outer(rw$center, cl$center, "+") + outer(rw$scale, cl$scale) * f$completed
  f
}

jobs <- list()
for (sd in SEEDS) for (input in c("raw", "biscale")) for (q in NCPS)
  jobs[[length(jobs) + 1]] <- list(kind = "holdout", seed = sd, input = input, ncp = q)
for (q in seq_len(GCV_MAX))
  jobs[[length(jobs) + 1]] <- list(kind = "gcv", seed = NA, input = "raw", ncp = q)

run_job <- function(j) {
  tag <- if (j$kind == "gcv") sprintf("gcv_ncp%02d", j$ncp)
         else sprintf("holdout_%d_%s_ncp%02d", j$seed, j$input, j$ncp)
  ck <- file.path(CK_DIR, paste0(tag, ".csv"))
  if (file.exists(ck) && !FORCE) return(read.csv(ck))
  if (j$kind == "gcv") {
    f <- fit_mm(X, j$ncp, isTRUE(MM$scale %||% TRUE))
    ## estim_ncpPCA(method.cv = "gcv"), samples x genes orientation
    Xs <- t(X); rec <- t(f$fitted); n <- nrow(Xs); p <- ncol(Xs); nna <- sum(is.na(Xs))
    crit <- mean(((n * p - nna) * (Xs - rec) / ((n - 1) * p - nna - j$ncp * (n + p - j$ncp - 1)))^2,
                 na.rm = TRUE)
    out <- data.frame(kind = "gcv", seed = NA, input = "raw", ncp = j$ncp, r = NA, rmse = NA,
                      mae = NA, gcv = crit, capped = f$capped, minutes = f$minutes)
  } else {
    mk <- readRDS(file.path(OUT_DIR, "masks", sprintf("mask_gene_dataset_block_%d.rds", j$seed)))
    if (!identical(mk$n_genes, nrow(X))) stop("selection mask drawn on a different matrix")
    Y <- X; Y[mk$mask_idx] <- NA
    f <- complete_mm(Y, j$ncp, j$input)
    pred <- f$completed[mk$mask_idx]; truth <- X[mk$mask_idx]
    out <- data.frame(kind = "holdout", seed = j$seed, input = j$input, ncp = j$ncp,
                      r = cor(pred, truth), rmse = sqrt(mean((pred - truth)^2)),
                      mae = mean(abs(pred - truth)), gcv = NA, capped = f$capped, minutes = f$minutes)
  }
  write.csv(out, ck, row.names = FALSE)
  log_msg(tag, sprintf(": %.1f min%s", out$minutes, if (out$capped) " (iteration cap)" else ""))
  out
}

log_msg(length(jobs), " imputePCA fits on ", CORES, " cores")
res <- mclapply(jobs, function(j) tryCatch(run_job(j), error = function(e) conditionMessage(e)),
                mc.cores = CORES, mc.preschedule = FALSE)
bad <- !sapply(res, is.data.frame)
if (any(bad)) stop("fits failed: ", paste(unique(unlist(res[bad])), collapse = "; "))
all_rows <- do.call(rbind, res)

hold <- all_rows[all_rows$kind == "holdout", ]
gcv  <- all_rows[all_rows$kind == "gcv", c("ncp", "gcv", "capped", "minutes")]
write.csv(hold, file.path(OUT_DIR, "missmda_by_mask.csv"), row.names = FALSE)
write.csv(gcv, file.path(OUT_DIR, "missmda_gcv.csv"), row.names = FALSE)

se <- function(x) sd(x) / sqrt(length(x))
summ <- do.call(rbind, lapply(split(hold, list(hold$input, hold$ncp), drop = TRUE), function(g)
  data.frame(input = g$input[1], ncp = g$ncp[1], n_masks = nrow(g),
             rmse_mean = mean(g$rmse), rmse_se = se(g$rmse), r_mean = mean(g$r),
             mae_mean = mean(g$mae), n_capped = sum(g$capped), minutes_mean = mean(g$minutes))))
summ <- summ[order(summ$input, summ$ncp), ]
best <- summ[which.min(summ$rmse_mean), ]                 # best input and ncp
sel_in <- summ[summ$input == best$input, ]
ncp_star <- min(sel_in$ncp[sel_in$rmse_mean <= best$rmse_mean + best$rmse_se])
summ$selected <- summ$input == best$input & summ$ncp == ncp_star
write.csv(summ, file.path(OUT_DIR, "missmda_ncp.csv"), row.names = FALSE)

g <- gcv[order(gcv$ncp), ]
gcv_ncp <- if (any(diff(g$gcv) > 0)) g$ncp[which(diff(g$gcv) > 0)[1]] else g$ncp[which.min(g$gcv)]
msg <- c(sprintf("selected_ncp: %d", ncp_star),
         sprintf("selected_input: %s", best$input),
         sprintf("rule: one-SE on gene_dataset_block RMSE over %d selection masks; min at %s input, ncp %d (%.4f +/- %.4f)",
                 length(SEEDS), best$input, best$ncp, best$rmse_mean, best$rmse_se),
         sprintf("gcv_ncp: %d (estim_ncpPCA rule: first local minimum over ncp 1..%d)", gcv_ncp, GCV_MAX),
         sprintf("fits at the iteration cap: %d of %d", sum(all_rows$capped), nrow(all_rows)))
writeLines(msg, file.path(OUT_DIR, "missmda_selected.txt"))
cat("\n", paste(msg, collapse = "\n"), "\n\n", sep = "")
print(summ, row.names = FALSE, digits = 4)
for (d in c("checkpoints", "reference", "values", "de")) {
  p <- file.path(OUT_DIR, d)
  if (dir.exists(p) && !length(list.files(p))) unlink(p, recursive = TRUE)
}
log_msg("done")

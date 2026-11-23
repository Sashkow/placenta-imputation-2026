#!/usr/bin/env Rscript
## Holdout calls of the dissection route on genes whose GSE100051 block is masked,
## split by whether the other dataset that holds both trimesters (GSE9984) still
## measures the gene. Reads the stored masks and the calls written by
## routes_dissection_and_harmonizr.R; computes nothing new.
## Routes without imputation, supplementary split. Run from the repository root:
##   Rscript scripts/validation/routes_holdout_by_two_group_dataset.R --config=config/main.yaml
## Output: data/pipeline/routes/partB_dissection_by_two_group_dataset.csv
HOLDOUT_DEFAULT_OUT <- "data/pipeline/routes"
source("scripts/validation/_holdout_common.R")
config$paths$reference_de <- file.path(config$paths$output, "difexp_softimpute_combat_ref.tsv")
D <- prepare_data(config)
SECOND <- "GSE9984"   # the only dataset besides the reference with both trimesters
rows <- list()
for (rep in 1:3) {
  tag <- sprintf("gene_dataset_block_rep%d", rep)
  mk <- readRDS(file.path("data/pipeline/holdout/main/masks", paste0(tag, ".rds")))
  new_mask <- matrix(FALSE, nrow(D$X), ncol(D$X), dimnames = dimnames(D$X)); new_mask[mk$mask_idx] <- TRUE
  second_masked <- setNames(rowSums(new_mask[, D$ds_of_sample[colnames(D$X)] == SECOND, drop = FALSE]) > 0, D$genes)
  calls <- read.csv(file.path(OUT_DIR, paste0("partB_calls_", tag, ".csv")), stringsAsFactors = FALSE)
  calls <- calls[!calls$ref_visible, ]
  calls$second_visible <- !second_masked[as.character(calls$gene)]
  ## the pattern that carries most real calls: only the three Affymetrix datasets remain
  left <- sapply(D$datasets, function(d) rowSums(!is.na(D$X[, D$ds_of_sample[colnames(D$X)] == d, drop = FALSE]) &
                                                   !new_mask[, D$ds_of_sample[colnames(D$X)] == d, drop = FALSE]) > 0)
  rownames(left) <- D$genes
  affy <- c("GSE122214", "GSE37901", "GSE9984")
  calls$affy_trio <- apply(left[as.character(calls$gene), , drop = FALSE], 1, function(r) setequal(colnames(left)[r], affy))
  rows[[rep]] <- data.frame(rep = rep, k_visible = calls$k_visible, second_visible = calls$second_visible,
                            true_deg = calls$sig_full %in% TRUE, call = calls$sig_hz_mod_ref %in% TRUE, affy_trio = calls$affy_trio)
}
x <- do.call(rbind, rows)
agg <- function(d, k) data.frame(k_visible = k, second_two_group_dataset_visible = d$second_visible[1], genes = nrow(d),
  true_degs = sum(d$true_deg), dissection_calls = sum(d$call), dissection_calls_true = sum(d$call & d$true_deg))
out <- do.call(rbind, c(
  lapply(split(x, list(x$k_visible, x$second_visible), drop = TRUE), function(d) agg(d, as.character(d$k_visible[1]))),
  lapply(split(x, x$second_visible), function(d) agg(d, "all"))))
out <- out[order(out$second_two_group_dataset_visible, out$k_visible), ]
out <- rbind(out, agg(x[x$affy_trio, ], "affymetrix_trio_only"))
write.csv(out, file.path(OUT_DIR, "partB_dissection_by_two_group_dataset.csv"), row.names = FALSE)
print(out, row.names = FALSE)

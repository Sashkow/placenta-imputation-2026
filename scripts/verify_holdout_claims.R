#!/usr/bin/env Rscript
##
## Audit every holdout/DE-derived number written into the article's main.tex /
## supplement.tex against the CSV it came from (companion port of the igea-r audit).
##
## output/claims_article.csv has one row per claim:
##   claim_id      short unique id
##   location      where the number appears (e.g. "main sec 2.2", "supp S5 Table S2")
##   claim_text    the sentence fragment as written
##   value         the number as written in the tex (plain number; percentages as 36.5)
##   source_file   CSV path relative to output/ (or absolute)
##   source_filter R expression selecting ONE row, evaluated inside the data frame,
##                 e.g. `scheme == "gene_dataset_block" & method == "softimpute" & stratum == "ref_visible"`
##   source_column column holding the value; may be an R expression over columns,
##                 e.g. `100 * sensitivity_all`
##   tolerance     absolute tolerance for the comparison
##   status        filled by this script: ok / MISMATCH / missing
##
## Usage: Rscript scripts/verify_holdout_claims.R [--out=data/pipeline/holdout/main]

`%||%` <- function(a, b) if (is.null(a)) b else a
.args <- commandArgs(TRUE)
get_arg <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), .args, value = TRUE)
  if (length(h)) sub(paste0("^--", n, "="), "", h[1]) else d }
script_dir <- { fa <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(normalizePath(sub("--file=", "", fa[1]))) else getwd() }
ALT_DIR <- normalizePath(file.path(script_dir, ".."))
OUT_DIR <- { p <- get_arg("out", "data/pipeline/holdout/main"); if (!grepl("^/", p)) file.path(ALT_DIR, p) else p }

claims <- read.csv(file.path(OUT_DIR, "claims_article.csv"), stringsAsFactors = FALSE,
                   colClasses = "character")
if (!nrow(claims)) { cat("No claims recorded yet.\n"); quit(save = "no") }

cache <- list()
load_src <- function(f) {
  p <- if (grepl("^/", f)) f else file.path(OUT_DIR, f)
  if (is.null(cache[[p]])) cache[[p]] <<- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  cache[[p]]
}

res <- lapply(seq_len(nrow(claims)), function(i) {
  cl <- claims[i, ]
  out <- cl; out$actual <- NA_character_; out$status <- "missing"
  src <- tryCatch(load_src(cl$source_file), error = function(e) NULL)
  if (is.null(src)) { out$status <- paste0("missing: ", cl$source_file); return(out) }
  d <- if (nzchar(cl$source_filter)) src[with(src, eval(parse(text = cl$source_filter))), , drop = FALSE] else src
  if (nrow(d) != 1) { out$status <- sprintf("filter selects %d rows", nrow(d)); return(out) }
  actual <- tryCatch(with(d, eval(parse(text = cl$source_column))), error = function(e) NA)
  actual <- suppressWarnings(as.numeric(actual))   # source columns may be character (mixed-type value columns)
  out$actual <- format(actual, digits = 6)
  tol <- as.numeric(cl$tolerance %||% 0)
  claimed <- as.numeric(gsub("[{},]", "", cl$value))
  out$status <- if (is.na(actual) || is.na(claimed)) "not numeric" else
    if (abs(actual - claimed) <= tol) "ok" else "MISMATCH"
  out
})
res <- do.call(rbind, res)
write.csv(res, file.path(OUT_DIR, "claims_article_checked.csv"), row.names = FALSE)

cat(sprintf("%-28s %-26s %12s %12s  %s\n", "claim_id", "location", "claimed", "actual", "status"))
for (i in seq_len(nrow(res)))
  cat(sprintf("%-28s %-26s %12s %12s  %s\n", res$claim_id[i], substr(res$location[i], 1, 26),
              res$value[i], res$actual[i], res$status[i]))
n_bad <- sum(res$status != "ok")
cat(sprintf("\n%d claims checked, %d ok, %d need attention\n", nrow(res), sum(res$status == "ok"), n_bad))
if (n_bad) quit(status = 1)

#!/usr/bin/env Rscript
#
# Regenerate the README "Key numbers" table from verified sources.
#
# Values come from two places, never from a literal in this file:
#   * replication_report.md    -- the "Reproduced" column, looked up by exact
#                                 claim label (written by verify_article_claims.R)
#   * de_sweep_numbers.csv     -- the full-window Prater comparison, which the
#                                 report does not cover (it carries only the
#                                 GA-matched variant); these rows are checked by
#                                 verify_holdout_claims.R
# Denominators that depend on an input file are read from that file.
#
# Usage: Rscript scripts/update_readme_key_numbers.R

report_path <- "replication_report.md"
readme_path <- "README.md"
sweep_path  <- "data/pipeline/holdout/main/de_sweep_numbers.csv"
bench_path  <- "data/references/qpcr_benchmark_genes.csv"

## --- report table -----------------------------------------------------------
## Claim labels may themselves contain pipes (the report writes |logFC|), so the
## row is split on " | " after the outer delimiters are stripped, and the field
## count is asserted rather than assumed.
lines <- readLines(report_path)
rows <- grep("^\\|.*\\|$", lines, value = TRUE)
rows <- rows[!grepl("^\\|[-| ]+\\|$", rows)]          # drop the header rule
fields <- lapply(rows, function(r) {
  f <- strsplit(sub("\\s*\\|\\s*$", "", sub("^\\s*\\|\\s*", "", r)), " \\| ")[[1]]
  trimws(f)
})
keep <- lengths(fields) == 6
if (any(!keep))
  warning("skipping ", sum(!keep), " malformed report row(s); expected 6 fields")
fields <- fields[keep]
fields <- Filter(function(f) f[1] != "Section", fields)   # drop the header row

claims  <- vapply(fields, function(f) f[2], character(1))
values  <- vapply(fields, function(f) f[4], character(1))
statuses <- vapply(fields, function(f) f[6], character(1))
lookup <- setNames(values, claims)
status_of <- setNames(statuses, claims)

v <- function(claim) {
  if (!claim %in% names(lookup))
    stop("claim label not present in ", report_path, ": ", claim)
  st <- unname(status_of[claim])
  if (!identical(st, "MATCH"))
    stop("claim '", claim, "' is ", st, " in ", report_path,
         " -- the verification failed, so the README must not quote it. ",
         "Investigate the number before regenerating the table.")
  unname(lookup[claim])
}

## --- data-sourced values ----------------------------------------------------
sweep <- read.csv(sweep_path, stringsAsFactors = FALSE)
s <- function(name) {
  hit <- sweep$value[sweep$name == name]
  if (length(hit) != 1) stop("expected exactly one '", name, "' row in ", sweep_path)
  as.numeric(hit)
}
n_benchmark <- nrow(read.csv(bench_path, stringsAsFactors = FALSE))

fmt <- function(x) prettyNum(x, big.mark = ",")

tbl <- c(
  "| Claim | Value | Verified by |",
  "|-------|-------|-------------|",
  sprintf("| Genes at intersection (k=6) | %s | replication_report.md |", fmt(v("Intersection gene count (k=6)"))),
  sprintf("| Genes after imputation | %s | replication_report.md |", fmt(v("Union gene count (k=1)"))),
  sprintf("| Safety-drop genes removed | %s | replication_report.md |", v("Safety-drop genes removed")),
  sprintf("| DEGs (softImpute + ComBat-ref) | %s (%s up / %s down) | replication_report.md |",
          v("SoftImpute+ComBat-ref total DEGs"), v("SoftImpute+ComBat-ref up"), v("SoftImpute+ComBat-ref down")),
  sprintf("| DEGs (intersection + ComBat-ref) | %s | replication_report.md |", v("Intersection-only (none+ref) total DEGs")),
  sprintf("| Gained DEGs | %s (%s testable only via imputation) | replication_report.md |",
          v("Gained DEGs (530 set minus 277 set)"), v("Gained DEGs testable only via imputation")),
  sprintf("| softImpute CV r (random / block / worst-case) | %s / %s / %s | replication_report.md |",
          v("softImpute random-cell r"), v("softImpute gene-dataset block r"), v("softImpute worst-case block r")),
  sprintf("| softImpute block RMSE / MAE | %s / %s | replication_report.md |",
          v("softImpute block RMSE"), v("softImpute block MAE")),
  sprintf("| Balanced reference DEGs / retention of full-run DEGs | %s / %s | replication_report.md |",
          v("Balanced reference DEG count"), v("Balanced retention of full DEGs")),
  sprintf("| Subsampling DEG retention N=10 / N=30 | %s%% / %s%% | replication_report.md |",
          v("Subsampling DEG retention at N=10"), v("Subsampling DEG retention at N=30")),
  sprintf("| Subsampling logFC CCC at N=10 | %s | replication_report.md |", v("Subsampling logFC CCC at N=10")),
  sprintf("| Split-half DEG retention / logFC CCC | %s%% / %s | replication_report.md |",
          v("Split-half DEG retention"), v("Split-half logFC CCC between halves")),
  sprintf("| GO BP terms / KEGG pathways (full 530) | %s / %s | replication_report.md |",
          v("GO BP terms (full 530)"), v("KEGG pathways (full 530)")),
  sprintf("| Gained genes in Lykhenko 2021 limma table | %s/%s | replication_report.md |",
          v("Gained genes present in Lykhenko 2021 limma table"), v("Gained DEGs (530 set minus 277 set)")),
  sprintf("| Same direction as Lykhenko 2021 / already FDR-significant there | %s%% / %s%% | replication_report.md |",
          v("Same direction of change (%)"), v("Already FDR-significant in Lykhenko 2021 (%)")),
  sprintf("| Holdout block masking: true DEGs lost / non-DEGs gained | %s%% / %s%% | replication_report.md |",
          v("Block masking: true DEGs lost (%)"), v("Block masking: non-DEGs gained (%)")),
  sprintf("| Holdout, GSE100051 hidden: true DEGs lost | %s%% | replication_report.md |",
          v("Block masking, GSE100051 hidden: lost (%)")),
  sprintf("| qPCR benchmark genes recovered | %s/%d | replication_report.md; denominator from %s |",
          v("qPCR benchmark genes recovered as DEGs"), n_benchmark, basename(bench_path)),
  sprintf("| GA-matched Prater r / CCC | %s / %s | replication_report.md |",
          v("GA-matched Prater Pearson r"), v("GA-matched Prater CCC")),
  sprintf("| Prater (full window) r / CCC / shared significant DEGs | %.3f / %.3f / %s | %s (checked by verify_holdout_claims.R) |",
          s("prater_twostep_r"), s("prater_twostep_ccc"),
          fmt(s("prater_twostep_replicated")), basename(sweep_path))
)

readme <- readLines(readme_path)
start <- grep("<!-- KEY_NUMBERS_START -->", readme)
end <- grep("<!-- KEY_NUMBERS_END -->", readme)
stopifnot(length(start) == 1, length(end) == 1, start < end)
out <- c(readme[1:start], "", tbl, "", readme[end:length(readme)])
writeLines(out, readme_path)
cat("Key numbers table regenerated:", length(tbl) - 2, "rows\n")

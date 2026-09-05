#!/usr/bin/env Rscript
#
# Regenerate the README "Key numbers" table from replication_report.md.
#
# Every value is taken from the report's "Reproduced" column, looked up by
# exact claim label, so the table can never drift from the verified numbers
# again. The one exception is the full-window Prater comparison, which
# verify_article_claims.R does not cover (only the GA-matched variant is in
# the report); its row is a documented constant sourced from
# scripts/fig_rnaseq_concordance.R output and flagged as such below.
#
# Usage: Rscript scripts/update_readme_key_numbers.R

report_path <- "replication_report.md"
readme_path <- "README.md"

lines <- readLines(report_path)
rows <- grep("^\\| .+ \\| .+ \\| .+ \\| .+ \\| .+ \\| MATCH \\|$", lines, value = TRUE)
parts <- strsplit(rows, "\\s*\\|\\s*")
claims <- vapply(parts, function(x) x[3], character(1))
values <- vapply(parts, function(x) x[5], character(1))
lookup <- setNames(values, claims)

v <- function(claim) {
  if (!claim %in% names(lookup))
    stop("claim not found in replication_report.md: ", claim)
  unname(lookup[claim])
}

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
  sprintf("| qPCR benchmark genes recovered | %s/17 | replication_report.md |",
          v("qPCR benchmark genes recovered as DEGs")),
  sprintf("| GA-matched Prater r / CCC | %s / %s | replication_report.md |",
          v("GA-matched Prater Pearson r"), v("GA-matched Prater CCC")),
  "| Prater (full window) r / CCC / shared significant DEGs | 0.670 / 0.548 / 367 | `scripts/fig_rnaseq_concordance.R` output (not covered by the report) |"
)

readme <- readLines(readme_path)
start <- grep("<!-- KEY_NUMBERS_START -->", readme)
end <- grep("<!-- KEY_NUMBERS_END -->", readme)
stopifnot(length(start) == 1, length(end) == 1, start < end)
out <- c(readme[1:start], "", tbl, "", readme[end:length(readme)])
writeLines(out, readme_path)
cat("Key numbers table regenerated:", length(tbl) - 2, "rows\n")

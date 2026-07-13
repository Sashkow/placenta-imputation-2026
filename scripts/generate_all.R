#!/usr/bin/env Rscript
# Regenerate all article figures.
# Run from repo root: Rscript scripts/generate_all.R

script_dir <- "scripts"
scripts <- list.files(script_dir, pattern = "^fig_.*\\.R$", full.names = TRUE)
scripts <- sort(scripts)

cat(sprintf("=== Regenerating all article figures (%d scripts) ===\n\n", length(scripts)))

for (i in seq_along(scripts)) {
  cat(sprintf("--- %d/%d: %s ---\n", i, length(scripts), basename(scripts[i])))
  source(scripts[i])
  cat("\n")
}

cat("=== All figures regenerated ===\n")

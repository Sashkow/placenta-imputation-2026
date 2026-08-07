## Context

The staircase visualization code exists in three copies:

1. **companion `scripts/pipeline/plot_na_staircase.R`** (458 lines) — full version with study-specific wrappers (`build_merged_from_config`, `plot_na_staircase`, `plot_na_staircase_from_config`)
2. **companion `scripts/lib/plot_na_staircase.R`** (333 lines) — improved core-only version with better axis rendering, extra `render_staircase()` parameters (`cex_axis`, `cex_lab`, `cex_ds`, `line_lab`, `line_ds`)
3. **igea-r `phase2b_direct_merge/plot_na_staircase.R`** — identical to #1

Both repos keep flat structure (no `scripts/lib/`). The improved version from #2 becomes the canonical `plot_na_staircase.R` in both repos.

## Goals / Non-Goals

**Goals:**
- Single canonical `plot_na_staircase.R` in both repos (the improved version)
- All callers updated to use `prepare_staircase()` + `render_staircase()` directly
- Delete `scripts/lib/` from companion (no more duplicates)
- Clean up vestigial external `logging_utils.R` sourcing in `run_phase2b.R`

**Non-Goals:**
- Extracting code into an R package
- Refactoring `imputation.R`, `normalization.R`, or `subsampling_helpers.R`
- Changing any pipeline behavior or output

## Decisions

### 1. Replace `plot_na_staircase()` calls with inline prepare + render

The old wrapper did three things: iterated over FDR methods, opened/closed PNG devices, and called prepare+render. Callers will do this directly.

**In `run_phase2b.R:1107`** (both repos), the current call:
```r
plot_na_staircase(staircase_matrix, staircase_sample_ds,
                  staircase_title, file.path(output_dir, "na_staircase.png"),
                  sample_group = staircase_group,
                  gene_fdr_list = gene_fdr_list,
                  baseline = config$phenotype$baseline,
                  contrast = config$phenotype$contrast,
                  coverage_threshold = selected_threshold / n_datasets_total)
```

Becomes a loop over `gene_fdr_list` with explicit `png()`/`dev.off()`, calling `prepare_staircase()` then `render_staircase()`. Roughly 12 lines replacing 1 call — acceptable for clarity.

### 2. Inline `build_merged_from_config()` into `plot_all_runs_staircase.R`

This 86-line function is called only from `plot_all_runs_staircase.R` in igea-r. It loads phenodata, expression files, applies filters, builds a merged matrix, and collects FDR vectors. Since it's the only caller, the logic moves directly into that script. The function is deleted from `plot_na_staircase.R`.

### 3. Keep fallback logging definitions, drop external sourcing attempt

`run_phase2b.R:39-68` tries two external paths for `logging_utils.R`, then falls back to inline definitions. The external paths never resolve from the companion repo and are vestigial in igea-r too. Keep only the inline definitions (lines 47-67), remove the sourcing attempt (lines 39-44).

### 4. Delete `scripts/lib/` entirely

After the staircase file is replaced in `scripts/pipeline/`, the `scripts/lib/` directory has no purpose:
- `plot_na_staircase.R` — superseded by the updated `scripts/pipeline/` version
- `subsampling_helpers.R` — a 25-line stub; the full version is in `scripts/pipeline/`

Update `fig_staircase.R:11` to source from `scripts/pipeline/plot_na_staircase.R`.

## Risks / Trade-offs

- **[Risk]** `plot_all_runs_staircase.R` becomes 86 lines longer after inlining → Acceptable; the function was single-use study-specific glue. Moving it to the caller is the right altitude.
- **[Risk]** `run_phase2b.R` staircase section grows from 1 line to ~12 → Acceptable; the explicit loop is clearer than a wrapper that hides PNG creation logic.
- **[Risk]** igea-r changes could conflict with in-flight work → Low; files being touched (`plot_na_staircase.R`, `run_phase2b.R`, `plot_all_runs_staircase.R`) are stable infrastructure, not active development targets.

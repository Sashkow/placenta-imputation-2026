## Why

The staircase visualization code exists in three copies across two repos (companion `scripts/pipeline/`, companion `scripts/lib/`, and igea-r `phase2b_direct_merge/`). The `scripts/lib/` version has an improved API (`prepare_staircase()` + `render_staircase()`) with better axis rendering and finer-grained parameters, but the pipeline version retains study-specific wrappers that callers still use. Consolidating now — while all copies are still identical at the core — prevents silent divergence and establishes one canonical staircase API.

## What Changes

- **BREAKING**: Remove `plot_na_staircase()`, `plot_na_staircase_from_config()`, and `build_merged_from_config()` from `plot_na_staircase.R` in both repos. Replace with the improved lib version's API.
- Update `run_phase2b.R` in both repos: replace `plot_na_staircase()` call with `prepare_staircase()` + `render_staircase()`.
- Update igea-r's `plot_all_runs_staircase.R`: inline the data-loading logic that was in `build_merged_from_config()` directly into the script.
- Delete `scripts/lib/` directory in companion repo (consolidate into `scripts/pipeline/`).
- Update companion's `fig_staircase.R` source path from `scripts/lib/` to `scripts/pipeline/`.
- Clean up vestigial `logging_utils.R` external sourcing in `run_phase2b.R` (both repos).

## Capabilities

### New Capabilities
- `staircase-api`: The canonical staircase visualization API — `prepare_staircase()`, `render_staircase()`, and supporting functions (`downsample_rows`, `fdr_to_color`, `load_staircase_colors`). Defines the public surface that all callers use.

### Modified Capabilities

(none — no existing specs)

## Impact

- **Companion repo**: `scripts/pipeline/plot_na_staircase.R`, `scripts/pipeline/run_phase2b.R`, `scripts/fig_staircase.R`, `scripts/lib/` (deleted)
- **igea-r**: `phase2b_direct_merge/plot_na_staircase.R`, `phase2b_direct_merge/run_phase2b.R`, `phase2b_direct_merge/plot_all_runs_staircase.R`
- **No R package dependencies change** — the staircase code uses only base R graphics

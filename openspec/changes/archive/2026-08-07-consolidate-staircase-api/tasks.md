## 1. Replace staircase file in companion repo

- [x] 1.1 Copy `scripts/lib/plot_na_staircase.R` (improved version) to `scripts/pipeline/plot_na_staircase.R`, replacing the old version
- [x] 1.2 Delete `scripts/lib/` directory entirely (both `plot_na_staircase.R` and `subsampling_helpers.R` stub)

## 2. Update companion callers

- [x] 2.1 Update `scripts/fig_staircase.R:11` source path from `scripts/lib/plot_na_staircase.R` to `scripts/pipeline/plot_na_staircase.R`
- [x] 2.2 Update `scripts/fig_validation.R:7` source path from `scripts/lib/subsampling_helpers.R` to `scripts/pipeline/subsampling_helpers.R`
- [x] 2.3 Update `scripts/pipeline/run_phase2b.R:1107-1113` — replace `plot_na_staircase()` call with a loop over `gene_fdr_list` using `prepare_staircase()` + `render_staircase()` with explicit `png()`/`dev.off()`

## 3. Clean up logging in companion run_phase2b.R

- [x] 3.1 In `scripts/pipeline/run_phase2b.R:39-68` — remove the external `logging_utils.R` sourcing attempt (lines 39-44), keep only the inline fallback definitions (lines 47-67)

## 4. Update igea-r staircase file

- [x] 4.1 Copy the new canonical `scripts/pipeline/plot_na_staircase.R` from companion to igea-r at `scripts/integrative_analysis/phase2b_direct_merge/plot_na_staircase.R`

## 5. Update igea-r callers

- [x] 5.1 Update igea-r `run_phase2b.R:1107-1113` — same `plot_na_staircase()` → `prepare_staircase()` + `render_staircase()` replacement as task 2.2
- [x] 5.2 Update igea-r `plot_all_runs_staircase.R:26` — inline the `build_merged_from_config()` logic (86 lines) directly into the script, replacing the single function call

## 6. Clean up logging in igea-r run_phase2b.R

- [x] 6.1 In igea-r `run_phase2b.R:39-68` — same logging cleanup as task 3.1

## 7. Verify

- [x] 7.1 Confirm `diff` between companion and igea-r `plot_na_staircase.R` produces no output
- [x] 7.2 Confirm `diff` between companion and igea-r `run_phase2b.R` produces no output (they should remain identical)
- [x] 7.3 Run companion pipeline to verify staircase PNGs are generated correctly
- [x] 7.4 Confirm `scripts/lib/` does not exist in companion repo

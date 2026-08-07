## ADDED Requirements

### Requirement: Canonical staircase API surface
The `plot_na_staircase.R` file SHALL export exactly these public functions: `prepare_staircase()`, `render_staircase()`, `downsample_rows()`, `fdr_to_color()`, `load_staircase_colors()`. No convenience wrappers (`plot_na_staircase`, `plot_na_staircase_from_config`, `build_merged_from_config`) SHALL exist in the file.

#### Scenario: File contains only canonical functions
- **WHEN** `plot_na_staircase.R` is sourced
- **THEN** the environment SHALL contain `prepare_staircase`, `render_staircase`, `downsample_rows`, `fdr_to_color`, `load_staircase_colors` and SHALL NOT contain `plot_na_staircase`, `plot_na_staircase_from_config`, or `build_merged_from_config`

### Requirement: Identical staircase file in both repos
The `plot_na_staircase.R` file SHALL be byte-identical between the companion repo (`scripts/pipeline/plot_na_staircase.R`) and igea-r (`scripts/integrative_analysis/phase2b_direct_merge/plot_na_staircase.R`).

#### Scenario: Files match after change
- **WHEN** both repos have been updated
- **THEN** `diff` between the two files SHALL produce no output

### Requirement: render_staircase supports fine-grained sizing parameters
`render_staircase()` SHALL accept optional parameters `cex_axis`, `cex_lab`, `cex_ds`, `line_lab`, `line_ds` for controlling axis text size and label positioning, in addition to the existing `cex_main` and `cex_legend`.

#### Scenario: Custom axis sizing
- **WHEN** `render_staircase()` is called with `cex_axis = 0.9, cex_ds = 0.55`
- **THEN** the staircase plot SHALL render with those sizing parameters applied to axis labels and dataset labels respectively

### Requirement: No duplicate library files in companion
The `scripts/lib/` directory SHALL NOT exist in the companion repo. All reusable pipeline code SHALL live in `scripts/pipeline/`.

#### Scenario: lib directory removed
- **WHEN** the change is complete
- **THEN** `scripts/lib/` SHALL not exist and `scripts/fig_staircase.R` SHALL source `scripts/pipeline/plot_na_staircase.R`

### Requirement: No external source attempts in run_phase2b.R
`run_phase2b.R` SHALL NOT attempt to source files from outside the project. The logging utility functions (`setup_logging`, `close_logging`, `archive_previous_results`) SHALL be defined inline.

#### Scenario: Logging functions defined without external sourcing
- **WHEN** `run_phase2b.R` is sourced from either repo
- **THEN** no `source()` call SHALL reference `../../utils/` or `scripts/utils/logging_utils.R`
- **THEN** `setup_logging`, `close_logging`, and `archive_previous_results` SHALL be defined in the script body

### Requirement: All callers use prepare_staircase + render_staircase
Every script that produces a staircase visualization SHALL call `prepare_staircase()` followed by `render_staircase()` directly, managing its own PNG device. No caller SHALL use the removed `plot_na_staircase()` wrapper.

#### Scenario: run_phase2b.R staircase output
- **WHEN** `run_phase2b.R` runs the staircase section
- **THEN** it SHALL iterate over `gene_fdr_list`, call `prepare_staircase()` and `render_staircase()` for each FDR method, and write one PNG per method

#### Scenario: plot_all_runs_staircase.R in igea-r
- **WHEN** `plot_all_runs_staircase.R` runs
- **THEN** it SHALL load data inline (previously done by `build_merged_from_config`) and call `prepare_staircase()` + `render_staircase()` directly

#### Scenario: fig_staircase.R in companion
- **WHEN** `fig_staircase.R` runs
- **THEN** it SHALL source from `scripts/pipeline/plot_na_staircase.R` and continue using `prepare_staircase()` + `render_staircase()` (no change to call pattern, only source path)

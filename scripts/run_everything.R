#!/usr/bin/env Rscript
#
# Master orchestration script.
#
# Two modes:
#
#   (default) FIGURES-ONLY -- regenerates every article figure and table from
#   the pre-computed outputs shipped in data/pipeline/, then verifies the
#   article's numeric claims against them. Minutes.
#
#     Rscript scripts/run_everything.R
#
#   --full -- regenerates data/pipeline/ itself: every pipeline config, every
#   validation test, the ported analyses, enrichment, then figures and
#   verification. Hours. Stages whose outputs already exist are SKIPPED unless
#   --force is given; every skip is logged, never silent.
#
#     Rscript scripts/run_everything.R --full
#     Rscript scripts/run_everything.R --full --force
#     Rscript scripts/run_everything.R --full --force --no-archive   # no archive/ copies of replaced outputs
#     Rscript scripts/run_everything.R --full --only=pipeline,validation
#     Rscript scripts/run_everything.R --full --only=holdout,analysis   # groups: selection, pipeline, validation, holdout, analysis, figures
#
# A note on reproducing --full exactly: the committed data/pipeline/ outputs
# predate the addition of a fixed RNG seed to the softImpute call, so a fresh
# --full run can differ from the committed tables for a handful of genes whose
# |logFC| sits within rounding distance of the 1.0 cutoff. Verification applies
# the documented tolerance rather than demanding an exact match. See README.

args <- commandArgs(TRUE)
has_flag <- function(f) any(args == paste0("--", f))
get_arg <- function(n, d = NULL) { h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (length(h)) sub(paste0("^--", n, "="), "", h[1]) else d }

FULL  <- has_flag("full")
FORCE <- has_flag("force")
ONLY  <- if (!is.null(get_arg("only"))) strsplit(get_arg("only"), ",")[[1]] else NULL
DRY   <- has_flag("dry-run")
## --no-archive: run_phase2b.R does not copy the previous outputs to
## <dir>/archive/<timestamp>/ (the committed outputs are already in git)
NO_ARCHIVE <- has_flag("no-archive")
## --skip=stage1,stage2: leave out individual stages by name
SKIP <- if (!is.null(get_arg("skip"))) strsplit(get_arg("skip"), ",")[[1]] else NULL

if (!file.exists("scripts/_common.R"))
  stop("run from the repository root: Rscript scripts/run_everything.R")

started <- Sys.time()
log_lines <- character(0)
say <- function(...) {
  msg <- paste0(...)
  log_lines <<- c(log_lines, msg)
  cat(msg, "\n", sep = "")
}
rule <- function(t) say("\n", strrep("=", 70), "\n", t, "\n", strrep("=", 70))

run_r <- function(path, extra = character(0)) {
  if (NO_ARCHIVE && basename(path) == "run_phase2b.R") extra <- c(extra, "--no_archive")
  ## the holdout scripts keep finished per-method results and checkpoints
  ## unless told otherwise, so --force must reach them too
  if (FORCE && basename(path) %in% c("holdout_reference_runs.R", "holdout_unified.R"))
    extra <- c(extra, "--force")
  cmd <- c(path, extra)
  say("  $ Rscript ", paste(cmd, collapse = " "))
  if (DRY) return(0L)
  st <- system2("Rscript", cmd, stdout = "", stderr = "")
  if (st != 0) say("  !! exited with status ", st)
  st
}

## A stage: name, group, the outputs that prove it ran, and how to run it.
stage <- function(name, group, outputs, run) list(name = name, group = group,
                                                  outputs = outputs, run = run)

## Hold-out choice of the softImpute lambda step (path_stop_step in every
## pipeline config); runs first because everything downstream uses it.
selection_stages <- list(
  stage("selection:lambda", "selection", "data/pipeline/lambda_selection/selected_step.txt",
        function() run_r("scripts/lambda_selection_path.R")),
  ## imputePCA comparator: ncp and input chosen on the same selection masks
  ## (ncp / input in config_holdout_*.yaml); hours, dominated by the EM fits
  stage("selection:missmda_ncp", "selection", "data/pipeline/lambda_selection/missmda_selected.txt",
        function() run_r("scripts/missmda_ncp_selection.R"))
)

pipeline_stages <- list(
  stage("pipeline:main", "pipeline", "data/pipeline/main/difexp_softimpute_combat_ref.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_pipeline.yaml")),
  stage("pipeline:balanced", "pipeline", "data/pipeline/balanced/difexp_softimpute_combat_ref.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_balanced.yaml")),
  stage("pipeline:batch_in_limma", "pipeline", "data/pipeline/batch_in_limma/difexp_softimpute_batch_in_limma.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_batch_in_limma.yaml")),
  stage("pipeline:balanced_batch_in_limma", "pipeline", "data/pipeline/balanced_batch_in_limma/difexp_softimpute_batch_in_limma.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_balanced_batch_in_limma.yaml")),
  stage("pipeline:covariate_categorical", "pipeline", "data/pipeline/covariate_categorical/difexp_softimpute_combat_ref.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_covariate_categorical.yaml")),
  stage("pipeline:covariate_poly2", "pipeline", "data/pipeline/covariate_poly2/difexp_softimpute_combat_ref.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_covariate_poly2.yaml")),
  stage("pipeline:covariate_ns3", "pipeline", "data/pipeline/covariate_ns3/difexp_softimpute_combat_ref.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_covariate_ns3.yaml")),
  stage("pipeline:ga_matched_prater", "pipeline", "data/pipeline/ga_matched_prater/difexp_softimpute_combat_ref.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_ga_matched_prater.yaml")),
  ## imputed/unimputed matrices that analysis:combat_sensitivity reads
  stage("pipeline:sensitivity", "pipeline", "data/pipeline/sensitivity/exprs_imputed_none.tsv",
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_sensitivity.yaml"))
)

validation_stages <- list(
  stage("validation:test1", "validation", "data/pipeline/validation/test1_first_trim_subsample.tsv",
        function() run_r("scripts/pipeline/test1_first_trim_subsample.R", "--config=config/config_validation.yaml")),
  stage("validation:test1b", "validation", "data/pipeline/validation/test1b_vs_balanced.tsv",
        function() run_r("scripts/pipeline/test1b_vs_balanced.R", "--config=config/config_validation.yaml")),
  stage("validation:test2", "validation", "data/pipeline/validation/test2_balanced_subsample.tsv",
        function() run_r("scripts/pipeline/test2_balanced_subsample.R", "--config=config/config_validation.yaml")),
  stage("validation:test3", "validation", "data/pipeline/validation/test3_split_half.tsv",
        function() run_r("scripts/pipeline/test3_split_half.R", "--config=config/config_validation.yaml")),
  stage("validation:test4_patterns", "validation", "data/pipeline/test4/coverage_patterns_gained.csv",
        function() run_r("scripts/test4_coverage_patterns.R"))
)
## test4_block_holdout_de.R, test4b_coverage_cv.R and test4b_refbatch_split.R are
## superseded by the unified holdout below (2026-08) and no longer run.

## Unified holdout validation (Supp S1 + S5, main-text holdout table): one set of
## masks feeds both the imputation-accuracy and the DE-call readouts.
holdout_stages <- list(
  stage("holdout:reference", "holdout", "data/pipeline/holdout/main/reference_runs_summary.csv",
        function() run_r("scripts/holdout_reference_runs.R")),
  stage("holdout:unified", "holdout", "data/pipeline/holdout/main/checkpoints/progressive_tax_block_rep3_batch_mean.done",
        function() run_r("scripts/holdout_unified.R")),
  stage("holdout:w1", "holdout", "data/pipeline/holdout/w1/checkpoints/progressive_tax_block_rep3_softimpute.done",
        function() { run_r("scripts/holdout_reference_runs.R", c("--config=config/config_holdout_w1.yaml", "--methods=softimpute"))
                     run_r("scripts/holdout_unified.R", c("--config=config/config_holdout_w1.yaml", "--methods=softimpute")) }),
  stage("holdout:plain_combat", "holdout", "data/pipeline/holdout/plain/checkpoints/progressive_tax_block_rep3_softimpute.done",
        function() { run_r("scripts/holdout_reference_runs.R", c("--config=config/config_holdout_plain_combat.yaml", "--methods=softimpute"))
                     run_r("scripts/holdout_unified.R", c("--config=config/config_holdout_plain_combat.yaml", "--methods=softimpute")) })
)

## subsampling stability for the three non-primary ComBat covariate variants
for (v in c("linear", "categorical", "poly2", "ns3")) {
  local({
    vv <- v
    validation_stages[[length(validation_stages) + 1]] <<- stage(
      paste0("validation:covariate_", vv), "validation",
      paste0("data/pipeline/validation_covariate_", vv, "/test1_within_size_jaccard.tsv"),
      function() run_r("scripts/pipeline/test1_first_trim_subsample.R",
                       paste0("--config=config/config_validation_covariate_", vv, ".yaml")))
  })
}

analysis_stages <- list(
  stage("analysis:lambda_selection_table", "analysis", "data/pipeline/lambda_selection/path.tex",
        function() run_r("scripts/lambda_selection_table.R")),
  stage("analysis:cv_comparators", "analysis", "data/pipeline/cv_comparators/table_s1_cv_comparators.csv",
        function() run_r("scripts/cv_comparators.R")),
  stage("analysis:covariate_sweep", "analysis", "data/pipeline/covariate_sweep/table_s2_covariate_comparison.csv",
        function() run_r("scripts/covariate_sweep.R")),
  stage("analysis:ga_matched_prater", "analysis", "data/pipeline/ga_matched_prater/prater_concordance.csv",
        function() run_r("scripts/ga_matched_prater.R")),
  stage("analysis:qpcr_recovery", "analysis", "data/pipeline/main/qpcr_recovery.csv",
        function() run_r("scripts/qpcr_recovery.R", "--out=data/pipeline/main/qpcr_recovery.csv")),
  stage("analysis:combat_sensitivity", "analysis", "data/pipeline/main/method_comparison.csv",
        function() run_r("scripts/pipeline/combat_sensitivity.R", "--config=config/config_sensitivity.yaml")),
  ## holdout summary tables (Supp Tables S1-S6, main-text holdout table) from the shipped holdout outputs
  stage("analysis:holdout_tables", "analysis", "data/pipeline/holdout/main/confusion_main.tex",
        function() { run_r("scripts/holdout_tables.R")
                     run_r("scripts/holdout_by_dataset.R", c("--dirs=data/pipeline/holdout/main,data/pipeline/holdout/plain", shQuote("--labels=ComBat-ref,plain ComBat")))
                     run_r("scripts/holdout_confusion.R")
                     run_r("scripts/holdout_confusion.R", c("--config=config/config_holdout_w1.yaml", "--w1_ref=data/pipeline/holdout/main/reference/softimpute_difexp.csv"))
                     run_r("scripts/holdout_confusion.R", "--config=config/config_holdout_plain_combat.yaml")
                     run_r("scripts/de_sweep_numbers.R") })
)

figure_stages <- list()
for (f in sort(list.files("scripts", pattern = "^fig_.*\\.R$"))) {
  local({
    ff <- f
    figure_stages[[length(figure_stages) + 1]] <<- stage(
      paste0("figure:", sub("\\.R$", "", ff)), "figures", NULL,
      function() run_r(file.path("scripts", ff)))
  })
}
figure_stages[[length(figure_stages) + 1]] <- stage(
  "figure:table_platform_coverage", "figures", NULL,
  function() run_r("scripts/table_platform_coverage.R"))

stages <- if (FULL) {
  c(selection_stages, pipeline_stages, validation_stages, holdout_stages, analysis_stages, figure_stages)
} else {
  ## figures-only: the ported summary analyses are cheap and read shipped data,
  ## so they run too; nothing that regenerates data/pipeline/ does.
  c(analysis_stages[sapply(analysis_stages, function(s)
      s$name %in% c("analysis:covariate_sweep", "analysis:ga_matched_prater",
                    "analysis:qpcr_recovery", "analysis:holdout_tables",
                    "analysis:lambda_selection_table"))],
    figure_stages)
}

if (!is.null(ONLY)) stages <- Filter(function(s) s$group %in% ONLY, stages)
if (!is.null(SKIP)) stages <- Filter(function(s) !s$name %in% SKIP, stages)

rule(sprintf("run_everything.R  |  mode: %s  |  %d stages",
             if (FULL) "FULL" else "figures-only", length(stages)))
if (FORCE) say("--force: existing outputs will be regenerated")
if (DRY)   say("--dry-run: commands are printed, not executed")

skipped <- character(0); failed <- character(0); ran <- character(0)

for (s in stages) {
  have <- !is.null(s$outputs) && all(file.exists(s$outputs))
  if (have && !FORCE) {
    say("[skip] ", s$name, "  (output exists: ", paste(s$outputs, collapse = ", "),
        "; use --force to regenerate)")
    skipped <- c(skipped, s$name)
    next
  }
  say("[run ] ", s$name)
  st <- tryCatch(s$run(), error = function(e) { say("  !! ", conditionMessage(e)); 1L })
  if (identical(st, 0L) || identical(st, 0)) ran <- c(ran, s$name) else failed <- c(failed, s$name)
}

rule("Claim verification")
ver <- tryCatch(run_r("scripts/verify_article_claims.R"),
                error = function(e) { say("  !! ", conditionMessage(e)); 1L })

rule("Summary")
say("ran:     ", length(ran))
say("skipped: ", length(skipped), if (length(skipped)) paste0("  [", paste(skipped, collapse = ", "), "]") else "")
say("failed:  ", length(failed), if (length(failed)) paste0("  [", paste(failed, collapse = ", "), "]") else "")
say("verifier exit status: ", ver)
say("elapsed: ", round(difftime(Sys.time(), started, units = "mins"), 1), " min")

writeLines(log_lines, "run_everything.log")
say("\nLog written to run_everything.log")
if (length(failed) || !identical(ver, 0L)) quit(status = 1)

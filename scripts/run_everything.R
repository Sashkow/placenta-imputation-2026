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
#     Rscript scripts/run_everything.R --full --only=pipeline,validation
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
        function() run_r("scripts/pipeline/run_phase2b.R", "--config=config/config_ga_matched_prater.yaml"))
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
        function() run_r("scripts/test4_coverage_patterns.R")),
  stage("validation:test4_block_holdout", "validation", "data/pipeline/test4/test4_summary.csv",
        function() run_r("scripts/test4_block_holdout_de.R")),
  stage("validation:test4b_coverage_cv", "validation", "data/pipeline/test4/test4b_summary.csv",
        function() run_r("scripts/test4b_coverage_cv.R"))
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
  stage("analysis:cv_comparators", "analysis", "data/pipeline/cv_comparators/table_s1_cv_comparators.csv",
        function() run_r("scripts/cv_comparators.R")),
  stage("analysis:covariate_sweep", "analysis", "data/pipeline/covariate_sweep/table_s2_covariate_comparison.csv",
        function() run_r("scripts/covariate_sweep.R")),
  stage("analysis:ga_matched_prater", "analysis", "data/pipeline/ga_matched_prater/prater_concordance.csv",
        function() run_r("scripts/ga_matched_prater.R")),
  stage("analysis:qpcr_recovery", "analysis", "data/pipeline/main/qpcr_recovery.csv",
        function() run_r("scripts/qpcr_recovery.R", "--out=data/pipeline/main/qpcr_recovery.csv")),
  stage("analysis:combat_sensitivity", "analysis", "data/pipeline/main/method_comparison.csv",
        function() run_r("scripts/pipeline/combat_sensitivity.R", "--config=config/config_sensitivity.yaml"))
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
  c(pipeline_stages, validation_stages, analysis_stages, figure_stages)
} else {
  ## figures-only: the ported summary analyses are cheap and read shipped data,
  ## so they run too; nothing that regenerates data/pipeline/ does.
  c(analysis_stages[sapply(analysis_stages, function(s)
      s$name %in% c("analysis:covariate_sweep", "analysis:ga_matched_prater",
                    "analysis:qpcr_recovery"))],
    figure_stages)
}

if (!is.null(ONLY)) stages <- Filter(function(s) s$group %in% ONLY, stages)

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

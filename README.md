# Companion code: Imputation-based integration of placental gene expression datasets

This repository contains the data, scripts, and pre-generated figures for:

> **Filling the gaps: matrix completion imputation for broader gene coverage when merging microarray gene expression datasets**
>
> Oleksandr Lykhenko, Yehor Polyakov, Maria Obolenskaya

## Repository structure

```
.
├── article/figures/          Pre-generated article figures (PNG, 300 DPI)
├── config/                   YAML configs (see "Configs" below)
├── data/
│   ├── phenodata.tsv         Sample metadata (186 rows, 8 datasets; the 6 analysed datasets contribute the 117 samples)
│   ├── expression/           Per-dataset expression matrices (ENTREZID-keyed)
│   ├── pipeline/             Pre-computed outputs (the authoritative record of the article's numbers)
│   │   ├── main/             6-dataset integration (softImpute + ComBat-ref)
│   │   ├── balanced/         Balanced 2-dataset reference
│   │   ├── validation/       Subsampling, split-half, permutation and leave-one-dataset-out results
│   │   ├── holdout/          Masked-cell holdout validation (main, weight1, plain)
│   │   └── ...               covariate variants, GA-matched cohort, lambda selection, census
│   └── references/           Lykhenko 2021 DEG lists, Prater 2021 supplementary tables
└── scripts/
    ├── run_everything.R      One entry point (see "Reproducing everything")
    ├── pipeline.R            THE METHOD, top to bottom (143 lines):
    │                         load -> gene coverage -> impute -> ComBat -> limma DE -> save
    ├── lib/                  One implementation of each step, shared by all callers
    │   ├── config.R                 load_config() with `base:` inheritance, CLI helpers
    │   ├── load_data.R              expression/phenodata loading, sample filters
    │   ├── gene_coverage.R          per-dataset gene coverage, per-group safety filter
    │   ├── imputation.R             softImpute plus the comparator imputers used in the holdout
    │   ├── batch_correction.R       ComBat model matrix, ComBat and ComBat-ref
    │   ├── differential_expression.R  limma design and fit (the only lmFit call)
    │   ├── masking.R                holdout masking schemes and leave-out imputation validation
    │   ├── metrics.R                Jaccard, Lin's CCC, subsample-vs-reference metrics
    │   ├── staircase_plot.R         missingness staircase plot
    │   ├── pipeline_report.R        PCA/staircase plots, method comparison, summary.txt, output archiving
    │   └── subsampling.R            data preparation and lean pipeline for the subsampling tests
    ├── validation/           Tests of the method
    │   ├── subsample_first_trimester.R      first-trimester size sweep vs the full run
    │   ├── subsample_vs_balanced_reference.R  same, against the balanced 2-dataset reference
    │   ├── subsample_balanced_cohort.R      balanced-cohort subsampling
    │   ├── split_half.R                     split-half replication
    │   ├── label_permutation.R              label-permutation null
    │   ├── leave_one_dataset_out.R          drop each dataset in turn
    │   ├── holdout_reference_runs.R         full-data run per imputer (the known answer)
    │   ├── holdout_masked_runs.R            masked runs: imputed values and DE calls
    │   ├── holdout_tables.R                 Supp Tables S1, S2, S5, S6
    │   ├── holdout_confusion.R              main-text holdout table, Supp Table S3
    │   ├── holdout_by_dataset.R             Supp Table S4: by hidden dataset, ComBat-ref vs plain
    │   └── _holdout_common.R, _lodo_permutation_common.R   shared setup
    ├── analysis/             Supporting analyses
    │   ├── covariate_sweep.R                Supp Tables S3/S4: ComBat covariate variants
    │   ├── ga_matched_prater.R              GA-matched RNA-seq concordance
    │   ├── qpcr_recovery.R                  17-gene qPCR benchmark (11/17)
    │   ├── combat_sensitivity.R             ComBat gamma/delta, intersection vs imputed
    │   ├── gene_coverage_census.R           which datasets measure which genes
    │   ├── softimpute_lambda_selection.R, softimpute_lambda_table.R
    │   ├── imputepca_ncp_selection.R
    │   └── table_platform_coverage.R        Table 6
    ├── figures/              fig_*.R (one per article figure) and _common.R
    └── verify/               verify_article_claims.R, verify_holdout_claims.R,
                              collect_article_numbers.R
```

### Configs

`main.yaml` is a complete config. Most others name a `base:` file and list only
what differs from it; `load_config()` merges them recursively. `main.yaml` is the
6-dataset pipeline, `balanced.yaml` the 2-dataset reference, `validation*.yaml` the
subsampling tests (the `validation_covariate_*` variants change only the ComBat covariate
model), `covariate_*.yaml` and `sensitivity_no_covariates.yaml` the ComBat covariate
sweep, `ga_matched_prater.yaml` the GA-matched cohort, `holdout*.yaml` the holdout
validation (`holdout_weight1`: imputed cells enter limma at weight 1; `holdout_plain_combat`:
ComBat without a reference batch), `lambda_selection.yaml` the softImpute lambda search.

## Reproducing everything

One entry point covers both the fast and the full case:

```bash
Rscript scripts/run_everything.R          # minutes: regenerate every figure and
                                          # table from the shipped data/pipeline/,
                                          # then verify the article's numbers
Rscript scripts/run_everything.R --full   # hours: regenerate data/pipeline/ itself
Rscript scripts/run_everything.R --full --force        # ignore existing outputs
Rscript scripts/run_everything.R --full --only=pipeline     # groups: pipeline, validation, holdout, analysis, figures
Rscript scripts/run_everything.R --dry-run             # print the plan, run nothing
```

Stages whose outputs already exist are skipped unless `--force`, and every skip is
printed — nothing is silently omitted. A run log is written to `run_everything.log`,
and the run ends with `verify_article_claims.R`. The sections below document the
individual steps that `run_everything.R` chains together.

## Reproducing figures

All scripts run from the repository root.

**Prerequisites:** R 4.5.0 (the version pinned in `renv.lock`, which also pins every package dependency).

**Determinism notes:** the softImpute ALS solver starts from a random matrix; the pipeline seeds it (`imputation.softimpute.seed: 42` in the config), so re-runs are repeatable. The committed outputs in `data/pipeline/` predate the seeding and were produced by an unseeded run — a fresh seeded run reproduces the article's numbers up to a few borderline genes at the |logFC| = 1 cutoff, and `verify_article_claims.R` checks the claims against the committed outputs, not a re-run. KEGG enrichment (`fig_enrichment.R`) queries the live KEGG database (accessed 2026-08-17) and is not version-pinned; the committed enrichment tables are the reference.

**Reproducibility disclaimer:** re-running the pipeline is not guaranteed to produce exactly the numbers reported in the article. Beyond the seeded softImpute step, R libraries may contain internal sources of randomness that the pipeline does not control, and floating-point results can differ across BLAS builds and thread counts. Package versions also drift: `renv.lock` pins the versions used here, but if you install newer releases instead of running `renv::restore()`, algorithmic changes in those releases may shift results slightly. Genes near the significance thresholds (adj.P.Val = 0.05, |logFC| = 1) are the most sensitive to such perturbations, so DEG counts may differ by a few genes even when effect-size estimates agree closely. The committed outputs in `data/pipeline/` are the authoritative record of what the article reports.

```bash
# Install packages (first time only, ~5-15 min)
Rscript -e 'renv::restore()'

# Regenerate all figures, tables and the claim check from the shipped data/pipeline/
Rscript scripts/run_everything.R

# Or regenerate a single figure
Rscript scripts/figures/fig_staircase.R
```

Output PNGs are written to `article/figures/`.

## Reproducing pipeline data

The `data/pipeline/` directory contains pre-computed outputs. To regenerate them from expression matrices and phenodata:

```bash
# Main integration (softImpute + ComBat-ref + limma DE)
Rscript scripts/pipeline.R --config=config/main.yaml

# ComBat sensitivity analysis (same pipeline, no covariates)
Rscript scripts/pipeline.R --config=config/sensitivity_no_covariates.yaml
Rscript scripts/analysis/combat_sensitivity.R --config=config/sensitivity_no_covariates.yaml

# Balanced 2-dataset reference (GSE100051 + GSE9984)
Rscript scripts/pipeline.R --config=config/balanced.yaml

# Subsampling validation (requires main + balanced outputs first)
Rscript scripts/validation/subsample_first_trimester.R      --config=config/validation.yaml
Rscript scripts/validation/subsample_vs_balanced_reference.R --config=config/validation.yaml
Rscript scripts/validation/subsample_balanced_cohort.R      --config=config/validation.yaml
Rscript scripts/validation/split_half.R                     --config=config/validation.yaml

# Label-permutation null: trimester labels (with gestational age) shuffled across
# all samples, 200 times
Rscript scripts/validation/label_permutation.R --n=200

# Leave-one-dataset-out: drop each dataset in turn and compare with the full DEG list
Rscript scripts/validation/leave_one_dataset_out.R

# ComBat covariate sweep -- Supplementary Tables S3 and S4
Rscript scripts/pipeline.R --config=config/covariate_categorical.yaml
Rscript scripts/pipeline.R --config=config/covariate_poly2.yaml
Rscript scripts/pipeline.R --config=config/covariate_ns3.yaml
for v in linear categorical poly2 ns3; do
  Rscript scripts/validation/subsample_first_trimester.R \
    --config=config/validation_covariate_$v.yaml
done
Rscript scripts/analysis/covariate_sweep.R

# GA-matched cohort for the Prater comparison -- Supplementary Section 6
Rscript scripts/pipeline.R --config=config/ga_matched_prater.yaml
Rscript scripts/analysis/ga_matched_prater.R

# Coverage census (which datasets measure which genes)
Rscript scripts/analysis/gene_coverage_census.R

# Holdout validation -- Supp Tables S1-S6 and the main-text holdout table.
# One set of masks (random cells / gene-dataset block / worst-case block; 10% of
# observed cells, 3 seeded repeats) is imputed by every method; the hidden cells are
# scored against the truth AND the imputed matrix is carried through ComBat and limma
# so the DE calls of the masked genes are scored against each imputer's own full run.
Rscript scripts/validation/holdout_reference_runs.R                       # full-data run per imputer
Rscript scripts/validation/holdout_masked_runs.R                          # 3 schemes x 3 repeats x 9 imputers (~4 h; imputePCA dominates)
Rscript scripts/validation/holdout_reference_runs.R --config=config/holdout_weight1.yaml --methods=softimpute
Rscript scripts/validation/holdout_masked_runs.R    --config=config/holdout_weight1.yaml --methods=softimpute      # imputed cells at weight 1
Rscript scripts/validation/holdout_reference_runs.R --config=config/holdout_plain_combat.yaml --methods=softimpute
Rscript scripts/validation/holdout_masked_runs.R    --config=config/holdout_plain_combat.yaml --methods=softimpute  # plain ComBat
Rscript scripts/validation/holdout_tables.R
Rscript scripts/validation/holdout_by_dataset.R --dirs=data/pipeline/holdout/main,data/pipeline/holdout/plain --labels="ComBat-ref,plain ComBat"
Rscript scripts/validation/holdout_confusion.R
Rscript scripts/validation/holdout_confusion.R --config=config/holdout_weight1.yaml --w1_ref=data/pipeline/holdout/main/reference/softimpute_difexp.csv
Rscript scripts/validation/holdout_confusion.R --config=config/holdout_plain_combat.yaml
Rscript scripts/verify/collect_article_numbers.R
Rscript scripts/verify/verify_holdout_claims.R                            # 296 article numbers vs their source CSVs
```

Holdout outputs live in `data/pipeline/holdout/{main,w1,plain}/`: `masks/` (the
saved masks), `reference/` (full-data DE tables and corrected matrices per imputer),
`values/` and `de/` (per scheme x repeat x imputer), and the summary tables
(`table_s1.csv`, `table_s3_arms.csv`, `table_s2.csv`, `confusion_full_vs_masked.csv`,
`by_dataset.csv`, `de_sweep_numbers.csv`, `claims_article.csv`). `run_everything.R
--full` runs these as the `holdout` stage group and `analysis:holdout_tables`.

Validation tests use parallel execution (`n_cores: 10` in config) and take ~30-60 minutes total.

Note that the covariate-variant subsampling runs write to
`data/pipeline/validation_covariate_*/`, deliberately separate from
`data/pipeline/validation/`, so the canonical validation outputs the article's
numbers derive from are never overwritten.

## Figures

Main text (`main.tex`) includes three PNG figures; the remainder are
supplementary (`supplement.tex`). The pipeline-overview diagram is TikZ
drawn inside `main.tex`. The two `fig_enrichment_*.png` files are generated for completeness but are not
included in either tex file.

| Figure | Content | Generating script | Input data |
|---|---------|-------------------|------------|
| Main 1 | Missingness staircase | `scripts/figures/fig_staircase.R` | `data/pipeline/main/exprs_imputed_softimpute.tsv`, `data/pipeline/main/difexp_softimpute_combat_ref.tsv` |
| Main 2 | PCA before/after ComBat | `scripts/figures/fig_pca_before_after_combat.R` | `data/pipeline/main/exprs_imputed_softimpute.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv`, `data/phenodata.tsv` |
| Main 3 | RNA-seq concordance (Prater) | `scripts/figures/fig_rnaseq_concordance.R` | `data/pipeline/main/difexp_softimpute_combat_ref.tsv`, `data/references/prater_2021_supp_tables.xlsx` |
| Supp S1 | Imputation CV range plot | `scripts/figures/fig_imputation_range.R` | `data/pipeline/holdout/main/table_s1.csv` |
| Supp S2 | PCA intersection vs softImpute | `scripts/figures/fig_pca_intersection_softimpute.R` | `data/pipeline/main/exprs_none_combat_ref.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv` |
| Supp S3 | ComBat sensitivity (scatter + MAE) | `scripts/figures/fig_combat_sensitivity.R` | `data/pipeline/main/exprs_none_combat_ref.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv` |
| Supp S4--S6 | Subsampling convergence / retention / CCC | `scripts/figures/fig_validation.R` | `data/pipeline/validation/test1_first_trim_subsample.tsv`, `data/pipeline/validation/test1b_vs_balanced.tsv` |
| Supp S7 | Split-half validation | `scripts/figures/fig_validation.R` | `data/pipeline/validation/test3_split_half.tsv` |
| Permutation null | Label-permutation null versus the real run | `scripts/figures/fig_permutation.R` | `data/pipeline/validation/permutation_null.tsv` |
| (extra) | GO / KEGG dotplots | `scripts/figures/fig_enrichment.R` | `data/pipeline/main/difexp_*.tsv` |

## Article tables

| # | Caption | Source data | How numbers were obtained |
|---|---------|-------------|--------------------------|
| 1 | Gene recovery & missingness (k thresholds) | `data/pipeline/main/summary.txt`, `imputed_gene_group_coverage.csv` | Pipeline run (`scripts/pipeline.R` with `config/main.yaml`) |
| 2 | Per-group safety audit | same as above | Same pipeline, per-group safety filter step |
| 3 | Imputation accuracy (leave-out CV) | `data/pipeline/main/imputation_validation.csv` | `scripts/lib/masking.R` (`validate_imputation()`), called from `scripts/pipeline.R` |
| 4 | DE counts (4 method combos + Lykhenko 2021) | `data/pipeline/main/difexp_significant_*.tsv`, `data/references/lykhenko_2021_deg.csv` | Pipeline DE step; Lykhenko reference from prior publication |
| 5 | ComBat covariate comparison | `data/pipeline/main/method_comparison.csv`, sensitivity pipeline outputs | `scripts/analysis/combat_sensitivity.R` with `config/sensitivity_no_covariates.yaml` |
| 6 | Per-dataset ENTREZID coverage | `data/expression/GSE*.tsv` | `scripts/analysis/table_platform_coverage.R` |
| 7 | Software & data resources | — | Manual/curated |

## Key numbers

<!-- KEY_NUMBERS_START -->

| Claim | Value | Verified by |
|-------|-------|-------------|
| Genes at intersection (k=6) | 8,260 | replication_report.md |
| Genes after imputation | 17,531 | replication_report.md |
| Safety-drop genes removed | 298 | replication_report.md |
| DEGs (softImpute + ComBat-ref) | 541 (448 up / 93 down) | replication_report.md |
| DEGs (intersection + ComBat-ref) | 277 | replication_report.md |
| Gained DEGs | 265 (237 testable only via imputation) | replication_report.md |
| softImpute CV r (random / block / worst-case) | 0.978 / 0.82 / 0.824 | replication_report.md |
| softImpute block RMSE / MAE | 1.5 / 1.15 | replication_report.md |
| Balanced reference DEGs / retention of full-run DEGs | 484 / 86.5% | replication_report.md |
| Subsampling DEG retention N=10 / N=30 | 90.9% / 93% | replication_report.md |
| Subsampling logFC CCC at N=10 | 0.89 | replication_report.md |
| Split-half DEG retention / logFC CCC | 86.4% / 0.771 | replication_report.md |
| GO BP terms / KEGG pathways (full 541) | 678 / 44 | replication_report.md |
| Gained genes in Lykhenko 2021 limma table | 245/265 | replication_report.md |
| Same direction as Lykhenko 2021 / already FDR-significant there | 89% / 60.8% | replication_report.md |
| Holdout block masking: true DEGs lost / non-DEGs gained | 26.6% / 0.13% | replication_report.md |
| Holdout, GSE100051 hidden: true DEGs lost | 96.8% | replication_report.md |
| qPCR benchmark genes recovered | 11/17 | replication_report.md; denominator from qpcr_benchmark_genes.csv |
| GA-matched Prater r / CCC | 0.808 / 0.755 | replication_report.md |
| Prater (full window) r / CCC / shared significant DEGs | 0.656 / 0.538 / 371 | de_sweep_numbers.csv (checked by verify/verify_holdout_claims.R) |

<!-- KEY_NUMBERS_END -->

## Verifying article claims

```bash
Rscript scripts/verify/verify_article_claims.R
```

Produces `replication_report.md` with a claim-by-claim comparison against companion data.

## Data sources

Expression data is derived from six GEO datasets:

| Dataset    | Platform  | Array                              | Samples | 1st | 2nd |
|------------|-----------|------------------------------------|---------|-----|-----|
| GSE100051  | GPL10558  | Illumina HumanHT-12 V4.0           | 49      | 42  | 7   |
| GSE122214  | GPL570    | Affymetrix HG-U133 Plus 2.0        | 4       | 4   | 0   |
| GSE28551   | GPL2986   | ABI Human Genome Survey v2         | 16      | 16  | 0   |
| GSE37901   | GPL570    | Affymetrix HG-U133 Plus 2.0        | 4       | 0   | 4   |
| GSE93520   | GPL6480   | Agilent Whole Human Genome 4x44K   | 36      | 36  | 0   |
| GSE9984    | GPL570    | Affymetrix HG-U133 Plus 2.0        | 8       | 4   | 4   |

Platform IDs and per-trimester counts are taken from `data/phenodata.tsv`
(healthy samples, first and second trimester only; 117 total).

The pipeline that produced the intermediate data files is available in the [main analysis repository](https://github.com/sashkow/integrative-gene-expression-analysis).

Reference tables in `data/references/`: `lykhenko_2021_deg.csv` and `lykhenko_2021_full_protein_coding.csv` are, respectively, the significant-DEG list and the full protein-coding limma table from the prior Affymetrix-only analysis (Lykhenko et al. 2021, 4 datasets, 22 samples); the full table backs the article's external-concordance claim (Discussion 3.2), since most gained genes fall below the 2021 significance threshold. `prater_2021_supp_tables.xlsx` is the published supplementary data of Prater et al. 2021.

## License

The **code** in this repository (`scripts/`, `config/`) is released under the
MIT License — see [LICENSE](LICENSE).

**Redistributed data** is not covered by that license and carries the terms of
its original source:

| Content | Source | Terms |
|---------|--------|-------|
| `data/expression/GSE*.tsv` | NCBI GEO series GSE100051, GSE122214, GSE28551, GSE37901, GSE93520, GSE9984 | Processed derivatives (probe-to-gene collapsed, ENTREZID-keyed) of publicly available GEO records; subject to the original depositors' terms. Cite the source studies, not this repository, when reusing the measurements. |
| `data/references/prater_2021_supp_tables.xlsx` | Supplementary tables of Prater et al. 2021 | Redistributed for verification of the concordance analysis under the publisher's terms for supplementary data; cite Prater et al. 2021. |
| `data/references/lykhenko_2021_*.csv` | Lykhenko et al. 2021 | Authors' own prior published results, redistributed here for comparison. |
| `data/pipeline/**` | Generated by this repository's scripts from the above inputs | MIT, same as the code. |

If you redistribute the GEO-derived matrices, retain the accession
attributions above.

# Companion code: Imputation-based integration of placental gene expression datasets

This repository contains the data, scripts, and pre-generated figures for:

> **Filling the gaps: matrix completion imputation for broader gene coverage when merging microarray gene expression datasets**
>
> Oleksandr Lykhenko, Yehor Polyakov, Maria Obolenskaya

## Repository structure

```
.
├── article/              Pre-generated article figures
│   └── figures/          Article figures (PNG, 300 DPI)
├── scripts/
│   ├── run_everything.R  One entry point: figures-only by default, --full to
│   │                     regenerate data/pipeline/ (see "Reproducing" below)
│   ├── fig_*.R           Figure scripts (one per article figure)
│   ├── table_platform_coverage.R  Table 6: per-dataset ENTREZID coverage
│   ├── generate_all.R    Regenerate all figures in one command
│   ├── verify_article_claims.R    Check every numeric claim against the data
│   ├── covariate_sweep.R          Supp Tables S3/S4: ComBat covariate variants
│   ├── cv_comparators.R           Earlier imputation-method comparison (superseded by holdout_*)
│   ├── ga_matched_prater.R        Supp S6: GA-matched RNA-seq concordance
│   ├── qpcr_recovery.R            17-gene benchmark recovery (11/17)
│   ├── test4_coverage_patterns.R  Which datasets measure which genes (coverage census)
│   │   (batch-in-limma variants: config_batch_in_limma.yaml, config_balanced_batch_in_limma.yaml,
│   │    config_validation_batch_in_limma.yaml -- exploratory one-step model, NOT reported in the article)
│   ├── holdout_common.R           Unified holdout: shared setup (masks, ComBat, limma)
│   ├── holdout_reference_runs.R   Full-data run per imputer (the known answer)
│   ├── holdout_unified.R          One mask, two readouts: imputed values + DE calls (M/F/O arms)
│   ├── holdout_tables.R           Supp Tables S1, S2, S5, S6 from the holdout outputs
│   ├── holdout_confusion.R        Main-text holdout table, Supp Table S3, weight-1 summary
│   ├── holdout_by_dataset.R       Supp Table S4: split by hidden dataset, ComBat-ref vs plain
│   ├── de_sweep_numbers.R         Every DE-derived article number as name/value rows
│   ├── verify_holdout_claims.R    Audit of the article's holdout/DE numbers (claims_article.csv)
│   ├── test4_block_holdout_de.R   SUPERSEDED by holdout_*.R (kept for the record)
│   ├── test4b_coverage_cv.R       SUPERSEDED by holdout_*.R (kept for the record)
│   ├── test4b_refbatch_split.R    SUPERSEDED by holdout_*.R (kept for the record)
│   ├── fig_imputation_range.R     Supp Fig S1: CV range plot (reads holdout table_s1.csv)
│   ├── _common.R         Shared constants (paths, dimensions, phenodata)
│   └── pipeline/         Data integration pipeline (reproduces data/pipeline/)
│       ├── run_phase2b.R             Main pipeline: merge → impute → ComBat → limma DE
│       ├── imputation.R              softImpute and other imputation methods
│       ├── normalization.R           ComBat and other batch correction methods
│       ├── combat_sensitivity.R      ComBat parameter comparison analysis
│       ├── subsampling_helpers.R     Shared validation functions
│       ├── test1_first_trim_subsample.R  Validation test 1
│       ├── test1b_vs_balanced.R          Validation test 1b
│       ├── test2_balanced_subsample.R    Validation test 2
│       └── test3_split_half.R            Validation test 3
├── data/
│   ├── phenodata.tsv     Sample metadata (186 rows, 8 datasets; the 6 analysed datasets contribute the 117 samples)
│   ├── expression/       Per-dataset expression matrices (ENTREZID-keyed)
│   ├── pipeline/         Pre-computed pipeline outputs
│   │   ├── main/         6-dataset integration (softImpute + ComBat-ref)
│   │   ├── balanced/     Balanced 2-dataset reference
│   │   └── validation/   Subsampling validation results
│   └── references/       External reference data
│       ├── lykhenko_2021_deg.csv           Lykhenko 2021 DEG list (310 significant genes)
│       ├── lykhenko_2021_full_protein_coding.csv  Lykhenko 2021 full limma table (16,889 genes)
│       └── prater_2021_supp_tables.xlsx    Prater 2021 RNA-seq DEGs
└── config/               Configuration files
    ├── config_pipeline.yaml     Main 6-dataset integration pipeline
    ├── config_sensitivity.yaml  ComBat sensitivity (no covariates)
    ├── config_balanced.yaml     2-dataset balanced reference
    ├── config_validation.yaml   Subsampling validation tests
    └── staircase_colors.yaml    Color scheme for staircase plot
```

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

# Regenerate all figures
Rscript scripts/generate_all.R

# Or regenerate a single figure
Rscript scripts/fig_staircase.R
```

Output PNGs are written to `article/figures/`.

## Reproducing pipeline data

The `data/pipeline/` directory contains pre-computed outputs. To regenerate them from expression matrices and phenodata:

```bash
# Main integration (softImpute + ComBat-ref + limma DE)
Rscript scripts/pipeline/run_phase2b.R --config=config/config_pipeline.yaml

# ComBat sensitivity analysis (same pipeline, no covariates)
Rscript scripts/pipeline/run_phase2b.R --config=config/config_sensitivity.yaml

# Balanced 2-dataset reference (GSE100051 + GSE9984)
Rscript scripts/pipeline/run_phase2b.R --config=config/config_balanced.yaml

# Subsampling validation (requires main + balanced outputs first)
Rscript scripts/pipeline/test1_first_trim_subsample.R --config=config/config_validation.yaml
Rscript scripts/pipeline/test1b_vs_balanced.R --config=config/config_validation.yaml
Rscript scripts/pipeline/test2_balanced_subsample.R --config=config/config_validation.yaml
Rscript scripts/pipeline/test3_split_half.R --config=config/config_validation.yaml

# ComBat covariate sweep -- Supplementary Tables S3 and S4
Rscript scripts/pipeline/run_phase2b.R --config=config/config_covariate_categorical.yaml
Rscript scripts/pipeline/run_phase2b.R --config=config/config_covariate_poly2.yaml
Rscript scripts/pipeline/run_phase2b.R --config=config/config_covariate_ns3.yaml
for v in linear categorical poly2 ns3; do
  Rscript scripts/pipeline/test1_first_trim_subsample.R \
    --config=config/config_validation_covariate_$v.yaml
done
Rscript scripts/covariate_sweep.R

# GA-matched cohort for the Prater comparison -- Supplementary Section 6
Rscript scripts/pipeline/run_phase2b.R --config=config/config_ga_matched_prater.yaml
Rscript scripts/ga_matched_prater.R

# Coverage census (which datasets measure which genes)
Rscript scripts/test4_coverage_patterns.R

# Unified holdout validation -- Supp Tables S1-S6 and the main-text holdout table.
# One set of masks (random cells / gene-dataset block / worst-case block; 10% of
# observed cells, 3 seeded repeats) is imputed by every method; the hidden cells are
# scored against the truth AND the imputed matrix is carried through ComBat and limma
# so the DE calls of the masked genes are scored against each imputer's own full run.
Rscript scripts/holdout_reference_runs.R                       # full-data run per imputer
Rscript scripts/holdout_unified.R                              # 3 schemes x 3 repeats x 5 imputers (~3 h; imputePCA dominates)
Rscript scripts/holdout_reference_runs.R --config=config/config_holdout_w1.yaml --methods=softimpute
Rscript scripts/holdout_unified.R        --config=config/config_holdout_w1.yaml --methods=softimpute      # imputed cells at weight 1
Rscript scripts/holdout_reference_runs.R --config=config/config_holdout_plain_combat.yaml --methods=softimpute
Rscript scripts/holdout_unified.R        --config=config/config_holdout_plain_combat.yaml --methods=softimpute  # plain ComBat
Rscript scripts/holdout_tables.R
Rscript scripts/holdout_by_dataset.R --dirs=data/pipeline/holdout/main,data/pipeline/holdout/plain --labels="ComBat-ref,plain ComBat"
Rscript scripts/holdout_confusion.R
Rscript scripts/holdout_confusion.R --config=config/config_holdout_w1.yaml --w1_ref=data/pipeline/holdout/main/reference/softimpute_difexp.csv
Rscript scripts/holdout_confusion.R --config=config/config_holdout_plain_combat.yaml
Rscript scripts/de_sweep_numbers.R
Rscript scripts/verify_holdout_claims.R                        # 272 article numbers vs their source CSVs
```

Holdout outputs live in `data/pipeline/holdout/{main,w1,plain}/`: `masks/` (the
saved masks), `reference/` (full-data DE tables and corrected matrices per imputer),
`values/` and `de/` (per scheme x repeat x imputer), and the summary tables
(`table_s1.csv`, `table_s3_arms.csv`, `table_s2.csv`, `confusion_full_vs_masked.csv`,
`by_dataset.csv`, `de_sweep_numbers.csv`, `claims_article.csv`). `run_everything.R
--full` runs these as the `holdout` stage group and `analysis:holdout_tables`;
the earlier `test4_block_holdout_de.R` / `test4b_*` stages were retired in 2026-08
when the unified design replaced them.

Validation tests use parallel execution (`n_cores: 10` in config) and take ~30-60 minutes total.

Note that the covariate-variant subsampling runs write to
`data/pipeline/validation_covariate_*/`, deliberately separate from
`data/pipeline/validation/`, so the canonical validation outputs the article's
numbers derive from are never overwritten.

## Figures

Main text (`main.tex`) includes three PNG figures; the remainder are
supplementary (`supplement.tex`). The pipeline-overview diagram is TikZ
drawn inside `main.tex`. `fig_venn_three_way.png` and the two
`fig_enrichment_*.png` files are generated for completeness but are not
included in either tex file.

| Figure | Content | Generating script | Input data |
|---|---------|-------------------|------------|
| Main 1 | Missingness staircase | `scripts/fig_staircase.R` | `data/pipeline/main/exprs_imputed_softimpute.tsv`, `data/pipeline/main/difexp_softimpute_combat_ref.tsv` |
| Main 2 | PCA before/after ComBat | `scripts/fig_pca_before_after_combat.R` | `data/pipeline/main/exprs_imputed_softimpute.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv`, `data/phenodata.tsv` |
| Main 3 | RNA-seq concordance (Prater) | `scripts/fig_rnaseq_concordance.R` | `data/pipeline/main/difexp_softimpute_combat_ref.tsv`, `data/references/prater_2021_supp_tables.xlsx` |
| Supp S1 | Imputation CV range plot | `scripts/fig_imputation_range.R` | `data/pipeline/holdout/main/table_s1.csv` |
| Supp S2 | PCA intersection vs softImpute | `scripts/fig_pca_intersection_softimpute.R` | `data/pipeline/main/exprs_none_combat_ref.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv` |
| Supp S3 | ComBat sensitivity (scatter + MAE) | `scripts/fig_combat_sensitivity.R` | `data/pipeline/main/exprs_none_combat_ref.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv` |
| Supp S4--S6 | Subsampling convergence / retention / CCC | `scripts/fig_validation.R` | `data/pipeline/validation/test1_first_trim_subsample.tsv`, `data/pipeline/validation/test1b_vs_balanced.tsv` |
| Supp S7 | Split-half validation | `scripts/fig_validation.R` | `data/pipeline/validation/test3_split_half.tsv` |
| (extra) | Three-way Venn | `scripts/fig_venn.R` | `data/pipeline/main/difexp_significant_*.tsv`, `data/references/lykhenko_2021_deg.csv` |
| (extra) | GO / KEGG dotplots | `scripts/fig_enrichment.R` | `data/pipeline/main/difexp_*.tsv` |

## Article tables

| # | Caption | Source data | How numbers were obtained |
|---|---------|-------------|--------------------------|
| 1 | Gene recovery & missingness (k thresholds) | `data/pipeline/main/summary.txt`, `imputed_gene_group_coverage.csv` | Pipeline run (`scripts/pipeline/run_phase2b.R` with `config/config_pipeline.yaml`) |
| 2 | Per-group safety audit | same as above | Same pipeline, per-group safety filter step |
| 3 | Imputation accuracy (leave-out CV) | `data/pipeline/main/imputation_validation.csv` | `scripts/pipeline/imputation.R` cross-validation |
| 4 | DE counts (4 method combos + Lykhenko 2021) | `data/pipeline/main/difexp_significant_*.tsv`, `data/references/lykhenko_2021_deg.csv` | Pipeline DE step; Lykhenko reference from prior publication |
| 5 | ComBat covariate comparison | `data/pipeline/main/method_comparison.csv`, sensitivity pipeline outputs | `scripts/pipeline/combat_sensitivity.R` with `config/config_sensitivity.yaml` |
| 6 | Per-dataset ENTREZID coverage | `data/expression/GSE*.tsv` | `scripts/table_platform_coverage.R` |
| 7 | Software & data resources | — | Manual/curated |

## Key numbers

<!-- KEY_NUMBERS_START -->

| Claim | Value | Verified by |
|-------|-------|-------------|
| Genes at intersection (k=6) | 8,260 | replication_report.md |
| Genes after imputation | 17,531 | replication_report.md |
| Safety-drop genes removed | 298 | replication_report.md |
| DEGs (softImpute + ComBat-ref) | 530 (436 up / 94 down) | replication_report.md |
| DEGs (intersection + ComBat-ref) | 277 | replication_report.md |
| Gained DEGs | 253 (243 testable only via imputation) | replication_report.md |
| softImpute CV r (random / block / worst-case) | 0.993 / 0.817 / 0.824 | replication_report.md |
| softImpute block RMSE / MAE | 1.52 / 1.16 | replication_report.md |
| Balanced reference DEGs / retention of full-run DEGs | 474 / 86.4% | replication_report.md |
| Subsampling DEG retention N=10 / N=30 | 90.6% / 92.7% | replication_report.md |
| Subsampling logFC CCC at N=10 | 0.89 | replication_report.md |
| Split-half DEG retention / logFC CCC | 86.1% / 0.763 | replication_report.md |
| GO BP terms / KEGG pathways (full 530) | 659 / 39 | replication_report.md |
| Gained genes in Lykhenko 2021 limma table | 233/253 | replication_report.md |
| Same direction as Lykhenko 2021 / already FDR-significant there | 89.3% / 62.7% | replication_report.md |
| Holdout block masking: true DEGs lost / non-DEGs gained | 24.8% / 0.21% | replication_report.md |
| Holdout, GSE100051 hidden: true DEGs lost | 87.1% | replication_report.md |
| qPCR benchmark genes recovered | 11/17 | replication_report.md; denominator from qpcr_benchmark_genes.csv |
| GA-matched Prater r / CCC | 0.808 / 0.755 | replication_report.md |
| Prater (full window) r / CCC / shared significant DEGs | 0.670 / 0.548 / 367 | de_sweep_numbers.csv (checked by verify_holdout_claims.R) |

<!-- KEY_NUMBERS_END -->

## Verifying article claims

```bash
Rscript scripts/verify_article_claims.R
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

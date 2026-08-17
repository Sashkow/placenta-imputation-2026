# Companion code: Imputation-based integration of placental gene expression datasets

This repository contains the data, scripts, and article source for:

> **Filling the gaps: matrix completion imputation for broader gene coverage in cross-dataset direct-merge expression analysis**
>
> Oleksandr Lykhenko, Yehor Polyakov, Maria Obolenskaya

## Repository structure

```
.
├── article/              Pre-generated article figures
│   └── figures/          Article figures (PNG, 300 DPI)
├── scripts/
│   ├── fig_*.R           Figure scripts (one per article figure)
│   ├── table_platform_coverage.R  Table 6: per-dataset ENTREZID coverage
│   ├── generate_all.R    Regenerate all figures in one command
│   ├── _common.R         Shared constants (paths, dimensions, phenodata)
│   ├── lib/              Helper functions sourced by figure scripts
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
│   ├── phenodata.tsv     Sample metadata (117 samples, 6 datasets)
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

## Reproducing figures

All scripts run from the repository root.

**Prerequisites:** R (>= 4.1). All R package dependencies are pinned in `renv.lock`.

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
```

Validation tests use parallel execution (`n_cores: 10` in config) and take ~30-60 minutes total.

## Article figures

| # | Caption | Generating script | Input data |
|---|---------|-------------------|------------|
| 1 | Missingness staircase | `scripts/fig_staircase.R` | `data/pipeline/main/exprs_imputed_softimpute.tsv`, `data/pipeline/main/difexp_softimpute_combat_ref.tsv` |
| 2 | PCA before/after ComBat | `scripts/fig_pca_before_after_combat.R` | `data/pipeline/main/exprs_imputed_softimpute.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv`, `data/phenodata.tsv` |
| 3 | PCA intersection vs softImpute | `scripts/fig_pca_intersection_softimpute.R` | `data/pipeline/main/exprs_none_combat_ref.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv` |
| 4 | ComBat sensitivity (scatter + MAE) | `scripts/fig_combat_sensitivity.R` | `data/pipeline/main/exprs_none_combat_ref.tsv`, `data/pipeline/main/exprs_softimpute_combat_ref.tsv` |
| 5 | Subsampling convergence | `scripts/fig_validation.R` | `data/pipeline/validation/test1_first_trim_subsample.tsv`, `data/pipeline/validation/test1b_vs_balanced.tsv` |
| 6 | DEG retention across sizes | `scripts/fig_validation.R` | same as above |
| 7 | Lin's CCC agreement | `scripts/fig_validation.R` | same as above |
| 8 | Split-half validation | `scripts/fig_validation.R` | `data/pipeline/validation/test3_split_half.tsv` |
| 9 | Three-way Venn | `scripts/fig_venn.R` | `data/pipeline/main/difexp_significant_*.tsv`, `data/references/lykhenko_2021_deg.csv` |
| 10 | RNA-seq concordance (Prater) | `scripts/fig_rnaseq_concordance.R` | `data/pipeline/main/difexp_softimpute_combat_ref.tsv`, `data/references/prater_2021_supp_tables.xlsx` |
| 11 | Pipeline overview diagram | manually created | — |

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

| Claim | Value | Derived from |
|-------|-------|-------------|
| Genes at intersection (k=6) | 8,260 | `data/pipeline/main/exprs_imputed_none.tsv` row count |
| Genes after imputation | 17,531 | `data/pipeline/main/exprs_imputed_softimpute.tsv` row count |
| Samples | 117 | `data/phenodata.tsv` row count |
| softImpute CV accuracy (random) | r=0.994, RMSE=0.330 | `data/pipeline/main/imputation_validation.csv` |
| softImpute CV accuracy (block) | r=0.816, RMSE=1.52 | same |
| DEGs (softImpute + ComBat-ref) | 538 | `data/pipeline/main/difexp_significant_softimpute_combat_ref.tsv` |
| DEGs (intersection + ComBat-ref) | 277 | `data/pipeline/main/difexp_significant_none_combat_ref.tsv` |
| Gained DEGs | 262 | difference of above two sets |
| ComBat sensitivity r | 0.9994 | `scripts/fig_combat_sensitivity.R` computation |
| Prater concordance | r=0.646, CCC=0.529, 370 shared DEGs | `scripts/fig_rnaseq_concordance.R` |
| Balanced reference DEGs | 484 | `data/pipeline/balanced/difexp_significant_softimpute_combat_ref.tsv` |
| Balanced overlap | 86.8% (467/538) | `data/pipeline/validation/test1b_vs_balanced.tsv` |
| Subsampling Jaccard @N=30 | 0.689 | `data/pipeline/validation/test1_first_trim_subsample.tsv` |
| Split-half Jaccard | 0.45 | `data/pipeline/validation/test3_split_half.tsv` |
| GO BP terms (full 538) | 678 | `scripts/fig_enrichment.R` output |
| KEGG pathways (full 538) | 43 | same |
| Gained genes in Lykhenko 2021 limma table | 242/262 | `scripts/verify_article_claims.R` vs `data/references/lykhenko_2021_full_protein_coding.csv` |
| Gained-gene direction concordance (Lykhenko 2021) | 88.8% (215/242), r=0.600 | same |
| Gained genes already FDR<0.05 in Lykhenko 2021 | 61.2% (148/242) | same |
| qPCR-validated genes recovered | 11/17 | `scripts/pipeline/combat_sensitivity.R` |

## Verifying article claims

```bash
Rscript scripts/verify_article_claims.R
```

Produces `replication_report.md` with a claim-by-claim comparison against companion data.

## Data sources

Expression data is derived from six GEO datasets:

| Dataset    | Platform    | Samples | Trimester |
|------------|-------------|---------|-----------|
| GSE100051  | GPL6244     | 49      | 1st + 2nd |
| GSE122214  | GPL6244     | 4       | 1st       |
| GSE28551   | GPL6947     | 16      | 1st       |
| GSE37901   | GPL6947     | 4       | 2nd       |
| GSE93520   | GPL6244     | 36      | 1st       |
| GSE9984    | GPL570      | 8       | 1st + 2nd |

The pipeline that produced the intermediate data files is available in the [main analysis repository](https://github.com/sashkow/integrative-gene-expression-analysis).

Reference tables in `data/references/`: `lykhenko_2021_deg.csv` and `lykhenko_2021_full_protein_coding.csv` are, respectively, the significant-DEG list and the full protein-coding limma table from the prior Affymetrix-only analysis (Lykhenko et al. 2021, 4 datasets, 22 samples); the full table backs the article's external-concordance claim (Discussion 3.2), since most gained genes fall below the 2021 significance threshold. `prater_2021_supp_tables.xlsx` is the published supplementary data of Prater et al. 2021.

## License

This work is provided for research reproducibility. Expression data originates from NCBI GEO and is subject to the original depositors' terms.

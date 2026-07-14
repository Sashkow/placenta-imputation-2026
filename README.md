# Companion code: Imputation-based integration of placental gene expression datasets

This repository contains the data, scripts, and article source for:

> **Imputation-based integration of placental gene expression datasets across pregnancy trimesters**
>
> Olexandr Lykhenko, Yehor Polyakov, Maria Obolenskaya

## Repository structure

```
.
├── article/              Pre-generated article figures
│   └── figures/          Article figures (PNG, 300 DPI)
├── scripts/
│   ├── fig_*.R           Figure scripts (one per article figure)
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
│   │   ├── sensitivity/  ComBat sensitivity analysis
│   │   ├── balanced/     Balanced 2-dataset reference
│   │   └── validation/   Subsampling validation results
│   └── references/       External reference data
│       ├── lykhenko_2021_deg.csv           Lykhenko 2021 DEG list
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

## License

This work is provided for research reproducibility. Expression data originates from NCBI GEO and is subject to the original depositors' terms.

# Companion code: Imputation-based integration of placental gene expression datasets

This repository contains the data, scripts, and article source for:

> **Imputation-based integration of placental gene expression datasets across pregnancy trimesters**
>
> Olexandr Lykhenko, Yehor Chernysh

## Repository structure

```
.
├── article/              LaTeX source and pre-generated figures
│   ├── main.tex
│   ├── main.bbl
│   └── figures/          Article figures (PNG, 300 DPI)
├── scripts/              R scripts that produce each figure
│   ├── _common.R         Shared constants (paths, dimensions, phenodata)
│   ├── fig_staircase.R               Figure 1: NA staircase
│   ├── fig_pca_intersection_softimpute.R  Figure 3: PCA comparison
│   ├── fig_pca_before_after_combat.R      Figure: PCA before/after ComBat
│   ├── fig_combat_sensitivity.R           Figure 4: ComBat sensitivity
│   ├── fig_validation.R                   Figure 5: Subsampling validation
│   ├── fig_venn.R                         Figure 6: DEG Venn diagram
│   ├── fig_rnaseq_concordance.R           Figure 7: RNA-seq concordance
│   ├── generate_all.R    Regenerate all figures in one command
│   └── lib/              Pipeline library functions (sourced by figure scripts)
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
    ├── staircase_colors.yaml    Color scheme for staircase plot
    ├── config_pipeline.yaml     Pipeline configuration
    └── config_validation.yaml   Validation configuration
```

## Reproducing figures

All scripts run from the repository root.

**Prerequisites:** R (>= 4.1) with packages: `limma`, `ggplot2`, `gridExtra`, `VennDiagram`, `openxlsx`, `yaml`, `org.Hs.eg.db`, `softImpute`, `sva`.

```bash
# Regenerate all figures
Rscript scripts/generate_all.R

# Or regenerate a single figure
Rscript scripts/fig_staircase.R
```

Output PNGs are written to `article/figures/`.

## Data sources

Expression data is derived from six GEO datasets:

| Dataset    | Platform    | Samples | Trimester |
|------------|-------------|---------|-----------|
| GSE100051  | GPL6244     | 49      | 1st + 2nd |
| GSE122214  | GPL6244     | 4       | 2nd       |
| GSE28551   | GPL6947     | 16      | 1st + 2nd |
| GSE37901   | GPL6947     | 4       | 1st       |
| GSE93520   | GPL6244     | 36      | 1st + 2nd |
| GSE9984    | GPL570      | 8       | 1st + 2nd |

The pipeline that produced the intermediate data files is available in the [main analysis repository](https://github.com/sashkow/integrative-gene-expression-analysis).

## License

This work is provided for research reproducibility. Expression data originates from NCBI GEO and is subject to the original depositors' terms.

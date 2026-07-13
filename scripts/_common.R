# Shared constants for article figure generation.
# Run from repo root: source("scripts/_common.R")

options(bitmapType = "cairo")

TW  <- 6.27   # A4 textwidth, 1in margins
DPI <- 300
PT  <- 11

fig_dir      <- "article/figures"
base_dir     <- "data/pipeline/main"
combat_dir   <- "data/pipeline/sensitivity"
val_dir      <- "data/pipeline/validation"
pheno_path   <- "data/phenodata.tsv"

datasets_6ds <- c("GSE100051", "GSE122214", "GSE28551",
                   "GSE37901", "GSE93520", "GSE9984")

pheno_full <- read.delim(pheno_path, stringsAsFactors = FALSE)
pheno <- pheno_full[pheno_full$dataset_id %in% datasets_6ds &
                    pheno_full$condition == "healthy" &
                    pheno_full$trimester %in% c("First trimester",
                                                "Second trimester"), ]

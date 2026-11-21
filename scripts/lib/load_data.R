# Loading the per-dataset expression matrices and the sample table.

SAMPLE_ID  <- "arraydatafile_exprscolumnnames"   # phenodata column: sample name
DATASET_ID <- "secondaryaccession"               # phenodata column: GEO dataset

# Sample table, with the columns renamed to the names the code uses
# (config: paths.column_map).
load_phenodata <- function(config) {
  path <- config$paths$phenodata
  pheno <- if (grepl("\\.tsv$", path)) read.delim(path, stringsAsFactors = FALSE)
           else read.csv(path, stringsAsFactors = FALSE)
  for (to in names(config$paths$column_map)) {
    from <- config$paths$column_map[[to]]
    if (from %in% colnames(pheno)) pheno[[to]] <- pheno[[from]]
  }
  pheno
}

# One genes x samples matrix (ENTREZID rows) per dataset in config$files$datasets
load_expression <- function(config) {
  exprs_list <- list()
  for (ds in config$files$datasets) {
    file <- config$files$file_map[[ds]] %||% paste0(ds, config$files$suffix %||% "", ".tsv")
    path <- file.path(config$paths$mapped_data, file)
    if (!file.exists(path)) { cat("  WARNING: file not found:", path, "\n"); next }
    exprs_list[[ds]] <- read.delim(path, row.names = 1, check.names = FALSE)
    cat("  Loaded", ds, ":", nrow(exprs_list[[ds]]), "genes x", ncol(exprs_list[[ds]]), "samples\n")
  }
  exprs_list
}

# Gestational age in weeks: the exact week if reported, else the mean of the
# reported range, else a trimester midpoint.
compute_ga_weeks <- function(pdata, group_column = "trimester") {
  midpoint <- c("first trimester" = 8, "second trimester" = 16, "third trimester" = 39)
  ga <- numeric(nrow(pdata))
  for (i in seq_len(nrow(pdata))) {
    week  <- suppressWarnings(as.numeric(pdata$fetus_week[i]))
    range <- as.character(pdata$fetus_range_week[i])
    ga[i] <- if (!is.na(week) && week > 0) week
      else if (!is.na(range) && nchar(range) > 0 && range != "0")
        mean(as.numeric(strsplit(range, ",")[[1]]), na.rm = TRUE)
      else unname(midpoint[tolower(as.character(pdata[[group_column]][i]))]) %|na|% 12
  }
  ga
}
`%|na|%` <- function(a, b) if (length(a) != 1 || is.na(a)) b else a

# Keep in every dataset only the given samples; datasets left empty are dropped.
keep_samples <- function(exprs_list, sample_ids) {
  exprs_list <- lapply(exprs_list, function(e)
    e[, intersect(colnames(e), sample_ids), drop = FALSE])
  exprs_list[sapply(exprs_list, ncol) > 0]
}

# Samples allowed by config$sample_filter (column -> allowed values) and, if
# set, by config$ga_filter (group -> [min week, max week]).
select_samples <- function(phenodata, config) {
  ok <- rep(TRUE, nrow(phenodata))
  for (col in names(config$sample_filter)) {
    if (!col %in% colnames(phenodata)) stop("sample_filter column not in phenodata: ", col)
    ok <- ok & phenodata[[col]] %in% config$sample_filter[[col]]
  }
  if (length(config$ga_filter) > 0) {
    ga <- compute_ga_weeks(phenodata, config$phenotype$group_column)
    group <- as.character(phenodata[[config$phenotype$group_column]])
    in_window <- rep(FALSE, nrow(phenodata))
    for (g in names(config$ga_filter)) {
      window <- as.numeric(unlist(config$ga_filter[[g]]))
      if (length(window) != 2) stop("ga_filter['", g, "'] must be [min, max]")
      in_window <- in_window | (!is.na(ga) & group == g & ga >= window[1] & ga <= window[2])
    }
    ok <- ok & in_window
  }
  phenodata[[SAMPLE_ID]][ok]
}

# ENTREZIDs annotated as protein-coding in org.Hs.eg.db
protein_coding_genes <- function() {
  db <- org.Hs.eg.db::org.Hs.eg.db
  types <- AnnotationDbi::select(db, keys = AnnotationDbi::keys(db, "ENTREZID"),
                                 columns = c("ENTREZID", "GENETYPE"), keytype = "ENTREZID")
  types$ENTREZID[types$GENETYPE == "protein-coding"]
}

# Expression matrices and sample table, filtered as the config asks.
load_datasets <- function(config) {
  phenodata  <- load_phenodata(config)
  exprs_list <- load_expression(config)
  exprs_list <- keep_samples(exprs_list, select_samples(phenodata, config))
  if (isTRUE(config$gene_filter$protein_coding_only)) {
    coding <- protein_coding_genes()
    exprs_list <- lapply(exprs_list, function(e) e[rownames(e) %in% coding, , drop = FALSE])
  }
  cat("Datasets:", length(exprs_list), " samples:", sum(sapply(exprs_list, ncol)), "\n")
  list(exprs_list = exprs_list, phenodata = phenodata)
}

# Sample table rows for the given samples, in that order, with empty covariate
# values turned into NA.
phenodata_for <- function(phenodata, sample_ids, covariates = character(0)) {
  pd <- phenodata[phenodata[[SAMPLE_ID]] %in% sample_ids, ]
  rownames(pd) <- pd[[SAMPLE_ID]]
  pd <- pd[sample_ids, ]
  for (cov in covariates) {
    if (!cov %in% colnames(pd)) stop("phenotype.covariate not in phenodata: ", cov)
    pd[[cov]][pd[[cov]] %in% c("", "_", "NA")] <- NA
  }
  pd
}

#' NA staircase diagram for merged expression matrices.
#'
#' Entry points: prepare_staircase() + render_staircase()

suppressPackageStartupMessages(library(yaml))

downsample_rows <- function(mat_bool, max_rows = 1600) {
  nr <- nrow(mat_bool)
  if (nr <= max_rows) return(mat_bool)
  idx <- round(seq(1, nr, length.out = max_rows))
  mat_bool[idx, , drop = FALSE]
}

load_staircase_colors <- function(path = NULL) {
  defaults <- list(
    background = "white", text_primary = "black", text_dim = "black",
    fdr_gradient = list(significant = "#004D40", nonsig = "#B2DFDB", min_alpha = 0.05),
    categories = list(
      no_fdr = "#E0F2F1", imputable_na = "#FFF176",
      group_excluded = "#B71C1C", coverage_excluded = "#F9A825",
      padding = "#9E9E9E"
    ),
    legacy = list(na = "#D93025", excluded = "#D93025")
  )
  if (is.null(path) || !file.exists(path)) return(defaults)
  cfg <- yaml::read_yaml(path)
  for (n in names(defaults)) {
    if (!is.null(cfg[[n]])) {
      if (is.list(defaults[[n]])) {
        for (m in names(cfg[[n]])) defaults[[n]][[m]] <- cfg[[n]][[m]]
      } else {
        defaults[[n]] <- cfg[[n]]
      }
    }
  }
  defaults
}

STAIRCASE_COLORS <- load_staircase_colors(Sys.getenv("STAIRCASE_COLORS", unset = NA))

fdr_to_color <- function(fdr,
                         sig_color = STAIRCASE_COLORS$fdr_gradient$significant,
                         nonsig_color = STAIRCASE_COLORS$fdr_gradient$nonsig) {
  min_a <- STAIRCASE_COLORS$fdr_gradient$min_alpha
  alpha <- ifelse(fdr < 0.001, 1.0,
           ifelse(fdr < 0.01, 1.0 - 0.5 * ((fdr - 0.001) / 0.009),
           ifelse(fdr < 0.05, 0.5 - 0.3 * ((fdr - 0.01) / 0.04),
                  0.2 * exp(-3 * (fdr - 0.05)))))
  alpha <- pmax(alpha, min_a)
  rgb_sig    <- col2rgb(sig_color)[, 1] / 255
  rgb_nonsig <- col2rgb(nonsig_color)[, 1] / 255
  r <- rgb_sig[1] * alpha + rgb_nonsig[1] * (1 - alpha)
  g <- rgb_sig[2] * alpha + rgb_nonsig[2] * (1 - alpha)
  b <- rgb_sig[3] * alpha + rgb_nonsig[3] * (1 - alpha)
  rgb(r, g, b)
}

prepare_staircase <- function(mat, sample_group = NULL, gene_fdr = NULL,
                              baseline = NULL, contrast = NULL,
                              coverage_threshold = NULL) {
  is_na <- is.na(mat)
  gene_na_count <- rowSums(is_na)
  nc <- ncol(mat)

  # Padding genes (added for uniform height across runs)
  is_padding <- grepl("^_pad_", rownames(mat))

  # Group-excluded: one comparison group has zero real values
  is_group_excluded <- rep(FALSE, nrow(mat))
  n_baseline_missing <- 0L
  n_contrast_missing <- 0L
  if (!is.null(sample_group) && !is.null(baseline) && !is.null(contrast)) {
    bl_idx <- which(sample_group == baseline)
    ct_idx <- which(sample_group == contrast)
    if (length(bl_idx) > 0 && length(ct_idx) > 0) {
      bl_real <- rowSums(!is_na[, bl_idx, drop = FALSE])
      ct_real <- rowSums(!is_na[, ct_idx, drop = FALSE])
      is_group_excluded <- !is_padding & (bl_real == 0 | ct_real == 0)
      n_baseline_missing <- sum(!is_padding & bl_real == 0 & ct_real > 0)
      n_contrast_missing <- sum(!is_padding & ct_real == 0 & bl_real > 0)
    }
  }

  # Coverage-excluded: below threshold (not padding or group-excluded)
  is_coverage_excluded <- rep(FALSE, nrow(mat))
  if (!is.null(coverage_threshold) && coverage_threshold > 0) {
    gene_coverage <- 1 - gene_na_count / nc
    is_coverage_excluded <- !is_padding & !is_group_excluded &
                            gene_coverage < coverage_threshold
  }

  gene_excluded <- is_padding | is_group_excluded | is_coverage_excluded

  # Primary: descending NA count. Secondary: ascending FDR (dark genes first).
  fdr_for_sort <- rep(Inf, nrow(mat))
  if (!is.null(gene_fdr)) {
    gn <- rownames(mat)
    matched <- gene_fdr[gn]
    has <- !is.na(matched)
    fdr_for_sort[has] <- matched[has]
  }
  gene_ord <- order(-gene_na_count, fdr_for_sort)
  sorted_na  <- is_na[gene_ord, , drop = FALSE]
  gene_names <- rownames(mat)[gene_ord]
  excluded   <- gene_excluded[gene_ord]
  padding_sorted     <- is_padding[gene_ord]
  grp_excl_sorted    <- is_group_excluded[gene_ord]
  cov_excl_sorted    <- is_coverage_excluded[gene_ord]

  # Rearrange: within each row push NAs left
  col_perm <- matrix(0L, nrow(sorted_na), nc)
  for (i in seq_len(nrow(sorted_na))) {
    col_perm[i, ] <- order(-sorted_na[i, ])
    sorted_na[i, ] <- sorted_na[i, col_perm[i, ]]
  }

  list(
    sorted_na  = sorted_na,
    col_perm   = col_perm,
    gene_names = gene_names,
    excluded   = excluded,
    padding          = padding_sorted,
    group_excluded   = grp_excl_sorted,
    coverage_excluded = cov_excl_sorted,
    nr_orig    = nrow(mat),
    nc         = nc,
    na_pct     = round(100 * sum(is_na) / length(is_na), 1),
    genes_full    = sum(gene_na_count == 0),
    genes_partial = sum(gene_na_count > 0 & gene_na_count < nc),
    genes_excluded = sum(gene_excluded),
    n_padding         = sum(is_padding),
    n_group_excluded  = sum(is_group_excluded),
    n_baseline_missing = n_baseline_missing,
    n_contrast_missing = n_contrast_missing,
    n_coverage_excluded = sum(is_coverage_excluded),
    sample_group  = sample_group,
    gene_fdr      = gene_fdr,
    baseline      = baseline,
    contrast      = contrast
  )
}

COL_NO_FDR            <- STAIRCASE_COLORS$categories$no_fdr
COL_NA                <- STAIRCASE_COLORS$legacy$na
COL_EXCLUDED          <- STAIRCASE_COLORS$legacy$excluded
COL_IMPUTABLE_NA      <- STAIRCASE_COLORS$categories$imputable_na
COL_GROUP_EXCLUDED    <- STAIRCASE_COLORS$categories$group_excluded
COL_COVERAGE_EXCLUDED <- STAIRCASE_COLORS$categories$coverage_excluded
COL_PADDING           <- STAIRCASE_COLORS$categories$padding

render_staircase <- function(staircase, title, sample_ds = NULL,
                             xlim = NULL, ylim = NULL,
                             show_legend = TRUE, cex_main = 1.1,
                             cex_legend = 0.5,
                             cex_axis = 0.6, cex_lab = 0.7,
                             cex_ds = 0.38, line_lab = 3.5,
                             line_ds = 0.8) {
  sorted_na <- downsample_rows(staircase$sorted_na, 1600)
  excluded  <- staircase$excluded
  padding   <- if (!is.null(staircase$padding)) staircase$padding else rep(FALSE, length(excluded))
  grp_excl  <- if (!is.null(staircase$group_excluded)) staircase$group_excluded else rep(FALSE, length(excluded))
  cov_excl  <- if (!is.null(staircase$coverage_excluded)) staircase$coverage_excluded else rep(FALSE, length(excluded))
  if (nrow(sorted_na) < length(excluded)) {
    idx <- round(seq(1, length(excluded), length.out = nrow(sorted_na)))
    excluded <- excluded[idx]
    padding  <- padding[idx]
    grp_excl <- grp_excl[idx]
    cov_excl <- cov_excl[idx]
  }
  gene_names_ds <- staircase$gene_names
  if (length(gene_names_ds) > nrow(sorted_na)) {
    idx <- round(seq(1, length(gene_names_ds), length.out = nrow(sorted_na)))
    gene_names_ds <- gene_names_ds[idx]
  }

  nr <- nrow(sorted_na)
  nc <- staircase$nc
  nr_orig <- staircase$nr_orig
  na_pct  <- staircase$na_pct
  has_categories <- any(padding) || any(grp_excl) || any(cov_excl)

  eff_xlim <- if (!is.null(xlim)) xlim else c(0, nc)
  eff_ylim <- if (!is.null(ylim)) ylim else c(0, nr)

  gf <- staircase$gene_fdr

  # Build per-gene observed color: FDR gradient or white (no FDR data)
  gene_color <- rep(COL_NO_FDR, nr)
  if (!is.null(gf) && length(gf) > 0) {
    fdr_vec <- gf[gene_names_ds]
    has_fdr <- !is.na(fdr_vec) & !padding
    if (any(has_fdr)) {
      gene_color[has_fdr] <- fdr_to_color(fdr_vec[has_fdr])
    }
  }

  # Build color matrix with category-specific NA colors
  if (has_categories) {
    col_mat <- matrix(COL_IMPUTABLE_NA, nr, nc)
    if (any(padding))  col_mat[padding, ] <- COL_PADDING
    for (i in which(grp_excl)) col_mat[i, sorted_na[i, ]] <- COL_GROUP_EXCLUDED
    for (i in which(cov_excl)) col_mat[i, sorted_na[i, ]] <- COL_COVERAGE_EXCLUDED
    for (i in which(!padding & rowSums(!sorted_na) > 0)) {
      col_mat[i, !sorted_na[i, ]] <- gene_color[i]
    }
  } else {
    col_mat <- matrix(COL_NA, nr, nc)
    if (any(excluded)) col_mat[excluded, ] <- COL_EXCLUDED
    observed <- !sorted_na & !excluded
    for (i in which(rowSums(observed) > 0)) {
      col_mat[i, observed[i, ]] <- gene_color[i]
    }
  }

  txt_col <- STAIRCASE_COLORS$text_primary
  dim_col <- STAIRCASE_COLORS$text_dim
  bg_col  <- STAIRCASE_COLORS$background

  plot(NA, xlim = eff_xlim, ylim = eff_ylim,
       xlab = "", ylab = "", main = title,
       xaxs = "i", yaxs = "i", axes = FALSE,
       cex.main = cex_main, col.main = txt_col)

  rimg <- as.raster(col_mat)
  rasterImage(rimg, 0, 0, nc, nr, interpolate = FALSE)
  box(col = dim_col)

  # Y-axis: round tick marks
  y_breaks <- pretty(c(0, nr_orig), n = 5)
  y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= nr_orig]
  y_scaled <- y_breaks * (nr / nr_orig)
  axis(2, at = y_scaled, labels = y_breaks, las = 1, cex.axis = cex_axis, col = dim_col, col.axis = txt_col)
  mtext("Genes", side = 2, line = line_lab, cex = cex_lab, col = dim_col)

  # X-axis: dataset boundary ticks and labels
  if (!is.null(sample_ds)) {
    ds_order <- unique(sample_ds)
    ds_counts <- table(factor(sample_ds, levels = ds_order))
    cumpos <- cumsum(as.numeric(ds_counts))
    midpos <- c(0, cumpos[-length(cumpos)]) + as.numeric(ds_counts) / 2
    axis(1, at = cumpos, labels = FALSE, col = dim_col, tcl = -0.3)
    ds_labels <- paste0(ds_order, "\n(n=", ds_counts, ")")
    mtext(ds_labels, side = 1, at = midpos, line = line_ds, cex = cex_ds, col = dim_col)
  } else {
    mtext(sprintf("Samples (%d)", nc), side = 1, line = 1.5, cex = 0.7, col = dim_col)
  }

  if (show_legend) {
    n_imp <- staircase$genes_partial -
      staircase$n_group_excluded -
      staircase$n_coverage_excluded
    final_size <- staircase$genes_full + n_imp

    usr <- par("usr")
    lx <- usr[2] - (usr[2] - usr[1]) * 0.01   # right edge
    ly <- usr[3] + (usr[4] - usr[3]) * 0.01   # bottom edge
    line_h <- (usr[4] - usr[3]) * 0.033
    swatch_w <- (usr[2] - usr[1]) * 0.025
    swatch_h <- line_h * 0.7
    pad_x <- (usr[2] - usr[1]) * 0.01
    num_col_x <- lx - pad_x  # right-aligned numbers column

    # Build legend entries: list of (label, fill, number_or_NULL)
    entries <- list()
    add_entry <- function(label, fill, number = NULL)
      entries[[length(entries) + 1L]] <<- list(l = label, f = fill, n = number)

    if (has_categories) {
      fdr_mid <- fdr_to_color(0.01)
      add_entry("Observed (FDR shading)", fdr_mid)
      add_entry("Observed (no DE result)", COL_NO_FDR)
      if (n_imp > 0)
        add_entry("Missing, imputable", COL_IMPUTABLE_NA)
      if (!is.null(staircase$n_group_excluded) && staircase$n_group_excluded > 0)
        add_entry("One group all-NA (dropped)", COL_GROUP_EXCLUDED)
      if (!is.null(staircase$n_coverage_excluded) && staircase$n_coverage_excluded > 0)
        add_entry("Coverage excluded", COL_COVERAGE_EXCLUDED)
      if (!is.null(staircase$n_padding) && staircase$n_padding > 0)
        add_entry("Not in run", COL_PADDING)
    } else {
      add_entry(sprintf("NA (%.1f%%)", na_pct), COL_NA)
      add_entry("Observed (no DE result)", COL_NO_FDR)
      if (staircase$genes_excluded > 0)
        add_entry("Excluded", COL_EXCLUDED,
                  formatC(staircase$genes_excluded, format = "d", big.mark = ","))
    }

    # Separator + summary rows (no swatch)
    add_entry("---", NA)
    add_entry("Fully observed genes",  NA,
              formatC(staircase$genes_full, format = "d", big.mark = ","))
    add_entry("Genes needing imputation", NA,
              formatC(staircase$genes_partial, format = "d", big.mark = ","))
    add_entry("  kept (both groups real)", NA,
              formatC(n_imp, format = "d", big.mark = ","))

    n_bl <- staircase$n_baseline_missing %||% 0L
    n_ct <- staircase$n_contrast_missing %||% 0L
    bl_short <- sub(" [Tt]rimester$", " Trim", staircase$baseline %||% "group 1")
    ct_short <- sub(" [Tt]rimester$", " Trim", staircase$contrast %||% "group 2")
    if (n_ct > 0 || n_bl > 0) {
      add_entry(sprintf("  dropped (all-NA in %s)", ct_short), NA,
                formatC(n_ct, format = "d", big.mark = ","))
      add_entry(sprintf("  dropped (all-NA in %s)", bl_short), NA,
                formatC(n_bl, format = "d", big.mark = ","))
    } else if (staircase$n_group_excluded > 0) {
      add_entry("  dropped (one group all-NA)", NA,
                formatC(staircase$n_group_excluded, format = "d", big.mark = ","))
    }
    if (staircase$n_coverage_excluded > 0)
      add_entry("  dropped (low coverage)", NA,
                formatC(staircase$n_coverage_excluded, format = "d", big.mark = ","))
    add_entry("Final merged matrix", NA,
              paste0(formatC(final_size, format = "d", big.mark = ","), " genes"))

    n_entries <- length(entries)
    # count non-separator entries for height
    n_visible <- sum(sapply(entries, function(e) e$l != "---"))
    box_h <- (n_visible + 1.2) * line_h
    box_w <- (usr[2] - usr[1]) * 0.46
    box_x0 <- lx - box_w
    box_y0 <- ly

    rect(box_x0, box_y0, lx, box_y0 + box_h,
         col = adjustcolor(bg_col, alpha.f = 0.93), border = dim_col, lwd = 0.5)

    row_i <- 0
    for (e in entries) {
      if (e$l == "---") {
        # draw separator line
        sep_y <- box_y0 + box_h - (row_i + 0.7) * line_h
        segments(box_x0 + pad_x, sep_y, lx - pad_x, sep_y,
                 col = dim_col, lwd = 0.5)
        row_i <- row_i + 0.4
        next
      }
      row_i <- row_i + 1
      y_pos <- box_y0 + box_h - row_i * line_h

      # swatch
      if (!is.na(e$f)) {
        rect(box_x0 + pad_x, y_pos - swatch_h / 2,
             box_x0 + pad_x + swatch_w, y_pos + swatch_h / 2,
             col = e$f, border = "black", lwd = 0.5)
        text(box_x0 + pad_x + swatch_w + pad_x * 0.5, y_pos,
             e$l, adj = 0, cex = cex_legend, col = txt_col)
      } else {
        text(box_x0 + pad_x, y_pos,
             e$l, adj = 0, cex = cex_legend, col = txt_col)
      }

      # right-aligned number
      if (!is.null(e$n)) {
        text(num_col_x, y_pos, e$n, adj = 1, cex = cex_legend,
             col = txt_col, family = "mono")
      }
    }

    # FDR gradient bar (bottom-left)
    if (!is.null(gf)) {
      bar_x0 <- usr[1] + (usr[2] - usr[1]) * 0.01
      bar_x1 <- usr[1] + (usr[2] - usr[1]) * 0.12
      bar_y1 <- usr[3] + (usr[4] - usr[3]) * 0.20
      bar_y0 <- usr[3] + (usr[4] - usr[3]) * 0.03
      bar_h  <- bar_y1 - bar_y0

      fdr_min_log <- 4
      fdr_max_log <- 0
      fdr_to_barfrac <- function(f) {
        nlog <- -log10(pmax(f, 1e-4))
        (nlog - fdr_max_log) / (fdr_min_log - fdr_max_log)
      }

      n_steps <- 64
      fdr_seq <- 10^(-seq(fdr_min_log, 0, length.out = n_steps))
      col_seq <- fdr_to_color(fdr_seq)
      rect(bar_x0, bar_y0, bar_x1, bar_y1, col = NA, border = dim_col, lwd = 0.5)
      for (k in seq_len(n_steps)) {
        frac0 <- (k - 1) / n_steps
        frac1 <- k / n_steps
        rect(bar_x0, bar_y1 - frac1 * bar_h, bar_x1, bar_y1 - frac0 * bar_h,
             col = col_seq[k], border = NA)
      }
      rect(bar_x0, bar_y0, bar_x1, bar_y1, col = NA, border = dim_col, lwd = 0.5)

      tck_fdr <- c(0.001, 0.01, 0.05, 0.25, 1.0)
      tck_y   <- bar_y0 + fdr_to_barfrac(tck_fdr) * bar_h
      segments(bar_x1, tck_y, bar_x1 + (bar_x1 - bar_x0) * 0.12, tck_y, lwd = 0.5, col = dim_col)
      text(bar_x1 + (bar_x1 - bar_x0) * 0.18, tck_y,
           labels = tck_fdr, adj = 0, cex = cex_legend * 0.85, col = txt_col)
      text((bar_x0 + bar_x1) / 2, bar_y1 + bar_h * 0.06,
           "FDR", cex = cex_legend, font = 2, col = txt_col)
    }
  }

}

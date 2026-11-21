#!/usr/bin/env Rscript
##
## Supplementary table of the softImpute lambda path: hold-out error per path
## step on the selection masks, from data/pipeline/lambda_selection/path.csv
## (scripts/analysis/softimpute_lambda_selection.R). Writes path.tex next to it.
##
## Usage: Rscript scripts/analysis/softimpute_lambda_table.R

DIR <- "data/pipeline/lambda_selection"
p <- read.csv(file.path(DIR, "path.csv"), stringsAsFactors = FALSE)
blk <- p[p$scheme == "gene_dataset_block", ]
rnd <- p[p$scheme == "random_cells", ]
rnd <- rnd[match(blk$step, rnd$step), ]

fmt <- function(m, s) sprintf("$%.3f \\pm %.3f$", m, s)
rows <- sprintf("%d & %s & %.3f & %s & %s & %s & %s \\\\",
  blk$step, format(round(blk$lambda_full, 1), nsmall = 1), blk$lambda_ratio,
  ifelse(blk$selected, sprintf("\\textbf{%d}", blk$rank_full), as.character(blk$rank_full)),
  fmt(blk$rmse_mean, blk$rmse_se), sprintf("$%.3f$", blk$r_mean), fmt(rnd$rmse_mean, rnd$rmse_se))
rows <- sub("\\\\\\\\$", "\\\\\\\\", rows)
tex <- c("\\begin{tabular}{@{}rrrrccc@{}}", "\\toprule",
  "\\textbf{Step} & $\\boldsymbol{\\lambda}$ & $\\boldsymbol{\\lambda/\\lambda_0}$ & \\textbf{Factors} & \\textbf{Block RMSE} & \\textbf{Block $r$} & \\textbf{Random-cell RMSE} \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(DIR, "path.tex"))
cat(tex, sep = "\n")

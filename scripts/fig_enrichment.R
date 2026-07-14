#!/usr/bin/env Rscript
# GO Biological Process and KEGG pathway enrichment for DEG sets.
# Reads pre-computed DE results; outputs enrichment tables and
# compareCluster dotplots to article/figures/.

source("scripts/_common.R")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

enrichment_dir <- file.path(base_dir, "enrichment")
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

sig_softimpute <- read.delim(file.path(base_dir,
  "difexp_significant_softimpute_combat_ref.tsv"), stringsAsFactors = FALSE)
sig_intersection <- read.delim(file.path(base_dir,
  "difexp_significant_none_combat_ref.tsv"), stringsAsFactors = FALSE)

full_genes <- as.character(sig_softimpute$gene)
intersection_genes <- as.character(sig_intersection$gene)
gained_genes <- setdiff(full_genes, intersection_genes)

cat("Full DEGs:", length(full_genes), "\n")
cat("Intersection DEGs:", length(intersection_genes), "\n")
cat("Gained DEGs:", length(gained_genes), "\n\n")

full_limma <- read.delim(file.path(base_dir,
  "difexp_softimpute_combat_ref.tsv"), stringsAsFactors = FALSE)
universe <- as.character(full_limma$gene)
cat("Background universe:", length(universe), "genes\n\n")

run_enrichment <- function(gene_list, name, universe) {
  cat("Running GO BP for:", name, "(", length(gene_list), "genes)\n")
  ego <- enrichGO(gene          = gene_list,
                  universe      = universe,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  cat("  Significant GO BP terms:", nrow(as.data.frame(ego)), "\n")

  cat("Running KEGG for:", name, "\n")
  ekegg <- enrichKEGG(gene         = gene_list,
                      universe     = universe,
                      organism     = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  cat("  Significant KEGG pathways:", nrow(as.data.frame(ekegg)), "\n\n")

  list(go = ego, kegg = ekegg)
}

enrich_full <- run_enrichment(full_genes,
  paste0("full_", length(full_genes)), universe)
enrich_intersection <- run_enrichment(intersection_genes,
  paste0("intersection_", length(intersection_genes)), universe)
enrich_gained <- run_enrichment(gained_genes,
  paste0("gained_", length(gained_genes)), universe)

save_enrichment <- function(result, prefix) {
  go_df <- as.data.frame(result$go)
  kegg_df <- as.data.frame(result$kegg)
  if (nrow(go_df) > 0)
    write.csv(go_df, file.path(enrichment_dir,
      paste0(prefix, "_GO_BP.csv")), row.names = FALSE)
  if (nrow(kegg_df) > 0)
    write.csv(kegg_df, file.path(enrichment_dir,
      paste0(prefix, "_KEGG.csv")), row.names = FALSE)
}

save_enrichment(enrich_full,
  paste0("enrichment_full_", length(full_genes)))
save_enrichment(enrich_intersection,
  paste0("enrichment_intersection_", length(intersection_genes)))
save_enrichment(enrich_gained,
  paste0("enrichment_gained_", length(gained_genes)))

gene_clusters <- list(
  "Full" = full_genes,
  "Intersection" = intersection_genes,
  "Gained" = gained_genes
)

cat("Running compareCluster GO BP...\n")
cc_go <- compareCluster(geneCluster = gene_clusters,
                        fun = "enrichGO",
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        universe = universe,
                        readable = TRUE)

cat("Running compareCluster KEGG...\n")
cc_kegg <- compareCluster(geneCluster = gene_clusters,
                          fun = "enrichKEGG",
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          universe = universe)

if (nrow(as.data.frame(cc_go)) > 0) {
  p_go <- dotplot(cc_go, showCategory = 15) +
    ggtitle("GO Biological Process: DEG set comparison") +
    theme(axis.text.y = element_text(size = 8))
  png(file.path(fig_dir, "fig_enrichment_GO_BP.png"),
      width = TW * DPI * 2, height = TW * DPI * 1.6, res = DPI)
  print(p_go)
  invisible(dev.off())
  cat("Saved GO BP compareCluster dot plot\n")
}

if (nrow(as.data.frame(cc_kegg)) > 0) {
  p_kegg <- dotplot(cc_kegg, showCategory = 15) +
    ggtitle("KEGG pathways: DEG set comparison") +
    theme(axis.text.y = element_text(size = 8))
  png(file.path(fig_dir, "fig_enrichment_KEGG.png"),
      width = TW * DPI * 2, height = TW * DPI * 1.3, res = DPI)
  print(p_kegg)
  invisible(dev.off())
  cat("Saved KEGG compareCluster dot plot\n")
}

cat("\n=== Enrichment summary ===\n")
cat(sprintf("Full %d: %d GO BP terms, %d KEGG pathways\n",
  length(full_genes), nrow(as.data.frame(enrich_full$go)),
  nrow(as.data.frame(enrich_full$kegg))))
cat(sprintf("Intersection %d: %d GO BP terms, %d KEGG pathways\n",
  length(intersection_genes), nrow(as.data.frame(enrich_intersection$go)),
  nrow(as.data.frame(enrich_intersection$kegg))))
cat(sprintf("Gained %d: %d GO BP terms, %d KEGG pathways\n",
  length(gained_genes), nrow(as.data.frame(enrich_gained$go)),
  nrow(as.data.frame(enrich_gained$kegg))))

cat("\nEnrichment tables saved to:", enrichment_dir, "\n")
cat("Figures saved to:", fig_dir, "\n")

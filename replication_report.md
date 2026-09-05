# Replication Report: Companion Repository vs Article Claims

**Generated:** 2026-09-05
**R version:** R version 4.5.0 (2025-04-11)
**Platform:** x86_64-pc-linux-gnu

## Environment

**renv R version:** 4.5.0

## Summary

- **Total claims verified:** 125
- **Matches:** 125
- **Mismatches:** 0
- **Not verified:** 0

## Comparison Table

| Section | Claim | Article | Reproduced | Tolerance | Status |
|---------|-------|---------|------------|-----------|--------|
| Table 1 | Intersection gene count (k=6) | 8260 | 8260 | exact | MATCH |
| Table 1 | Union gene count (k=1) | 17531 | 17531 | exact | MATCH |
| Table 1 | Safety-drop genes removed | 298 | 298 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref total DEGs | 530 | 530 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref up | 436 | 436 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref down | 94 | 94 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) total DEGs | 277 | 277 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) up | 224 | 224 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) down | 53 | 53 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) total DEGs | 525 | 525 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) up | 433 | 433 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) down | 92 | 92 | exact | MATCH |
| Supp S1 | softImpute random-cell r | 0.993 | 0.993 | 0.0005 | MATCH |
| Supp S1 | softImpute gene-dataset block r | 0.817 | 0.817 | 0.0005 | MATCH |
| Supp S1 | softImpute block RMSE | 1.52 | 1.52 | 0.005 | MATCH |
| Supp S1 | softImpute block MAE | 1.16 | 1.16 | 0.005 | MATCH |
| Supp S1 | softImpute worst-case block r | 0.824 | 0.824 | 0.0005 | MATCH |
| Supp S1 | KNN block r | 0.683 | 0.683 | 0.0005 | MATCH |
| Supp S1 | imputePCA worst-case r | 0.436 | 0.436 | 0.0005 | MATCH |
| Section 3.4 | Jaccard at N=10 | 0.52 | 0.52 | 0.01 | MATCH |
| Section 3.4 | Split-half median Jaccard between halves | 0.445 | 0.445 | 0.01 | MATCH |
| Supp S3 | Split-half DEG retention | 86.1 | 86.1 | 0.1 | MATCH |
| Supp S3 | Split-half logFC CCC between halves | 0.763 | 0.763 | 0.01 | MATCH |
| Supp S3 | Subsampling DEG retention at N=10 | 90.6 | 90.6 | 0.1 | MATCH |
| Supp S3 | Subsampling DEG retention at N=30 | 92.7 | 92.7 | 0.1 | MATCH |
| Supp S3 | Subsampling logFC CCC at N=10 | 0.89 | 0.89 | 0.01 | MATCH |
| Supp S3 | Subsampling logFC CCC at N=100 | 1.00 | 1 | 0.01 | MATCH |
| Section 3.4 | Jaccard at N=100 | 0.98 | 0.98 | 0.01 | MATCH |
| Section 3.4 | Balanced reference DEG count | 474 | 474 | exact | MATCH |
| Section 3.4 | Balanced vs full Jaccard | 0.839 | 0.839 | 0.001 | MATCH |
| Section 3.4 | Balanced retention of full DEGs | 86.4% | 86.4% | exact | MATCH |
| Results 2.3 | Shared genes compared | 8260 | 8260 | exact | MATCH |
| Results 2.3 | Post-ComBat Pearson r (all cells) | 0.9999 | 0.9999 | 0.0005 | MATCH |
| Results 2.3 | Mean abs diff | 0.01 | 0.01 | 0.005 | MATCH |
| Results 2.3 | Median per-gene MAE | 0.009 | 0.009 | 0.005 | MATCH |
| Enrichment | 530-DEG GO BP terms (q<0.05) | 659 | 659 | exact | MATCH |
| Enrichment | 530-DEG KEGG pathways (q<0.05) | 39 | 39 | exact | MATCH |
| Enrichment | 277-DEG GO BP terms | 292 | 292 | exact | MATCH |
| Enrichment | 277-DEG KEGG pathways | 25 | 25 | exact | MATCH |
| Enrichment | 253-DEG GO BP terms | 250 | 250 | exact | MATCH |
| Enrichment | 253-DEG KEGG pathways | 17 | 17 | exact | MATCH |
| Discussion 3.2 | Gained DEGs (530 set minus 277 set) | 253 | 253 | exact | MATCH |
| Discussion 3.2 | Gained genes present in Lykhenko 2021 limma table | 233 | 233 | exact | MATCH |
| Discussion 3.2 | Same direction of change (%) | 89.3 | 89.3 | 0.05 | MATCH |
| Discussion 3.2 | logFC Pearson r vs Lykhenko 2021 | 0.654 | 0.654 | 0.0051 | MATCH |
| Discussion 3.2 | Already FDR-significant in Lykhenko 2021 (%) | 62.7 | 62.7 | 0.05 | MATCH |
| Supp S4 | Direction-discordant gained genes | 25 | 25 | exact | MATCH |
| Supp S4 | Discordant genes FDR-significant in 2021 | 0 | 0 | exact | MATCH |
| Supp S4 | Median 2021 |logFC| of discordant genes | 0.13 | 0.13 | 0.005 | MATCH |
| Supp S4 | Median 2021 |logFC| of concordant genes | 0.60 | 0.6 | 0.005 | MATCH |
| Supp S4 | 2021-significant matched genes agreeing in direction | 146/146 | 146/146 | exact | MATCH |
| Section 3.3 | GO BP terms (full 530) | 659 | 659 | exact | MATCH |
| Section 3.3 | KEGG pathways (full 530) | 39 | 39 | exact | MATCH |
| Section 3.3 | GO BP terms (intersection 277) | 292 | 292 | exact | MATCH |
| Section 3.3 | KEGG pathways (intersection 277) | 25 | 25 | exact | MATCH |
| Section 3.3 | GO BP terms (gained 253) | 250 | 250 | exact | MATCH |
| Section 3.3 | KEGG pathways (gained 253) | 17 | 17 | exact | MATCH |
| Section 3.3 | q: positive regulation of cytokine production (full) | 7e-14 | 7e-14 | 0.000000000000001 | MATCH |
| Section 3.3 | q: leukocyte migration (full) | 6.7e-12 | 6.7e-12 | 0.000000000001 | MATCH |
| Section 3.3 | q: chemotaxis (full) | 2.8e-11 | 2.8e-11 | 0.000000000001 | MATCH |
| Section 3.3 | q: humoral immune response (full) | 5.3e-11 | 5.3e-11 | 0.00000000001 | MATCH |
| Section 3.3 | q: hsa05150 S. aureus infection (full) | 2.2e-10 | 2.2e-10 | 0.00000000001 | MATCH |
| Section 3.3 | q: hsa04610 complement and coagulation (full) | 1e-09 | 1e-09 | 0.0000000001 | MATCH |
| Section 3.3 | q: hsa04514 cell adhesion molecules (full) | 7.7e-06 | 7.7e-06 | 0.0000001 | MATCH |
| Section 3.3 | q: chemotaxis (gained) | 2.4e-07 | 2.4e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: leukocyte migration (gained) | 2.4e-07 | 2.4e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: hsa04145 phagocytosis (gained) | 0.0023 | 0.0023 | 0.0001 | MATCH |
| Section 3.4 | Gained GO BP terms shared with intersection | 83 | 83 | exact | MATCH |
| Section 3.4 | Gained KEGG pathways shared with intersection | 7 | 7 | exact | MATCH |
| Section 3.3 | Gained DEGs that are intersection genes | 10 | 10 | exact | MATCH |
| Section 3.3 | Gained DEGs testable only via imputation | 243 | 243 | exact | MATCH |
| Supp S6 | qPCR benchmark genes recovered as DEGs | 11 | 11 | exact | MATCH |
| Supp S5 | GA-matched testable genes | 11761 | 11761 | exact | MATCH |
| Supp S5 | GA-matched DEGs | 602 | 602 | exact | MATCH |
| Supp S5 | GA-matched shared genes with Prater | 1992 | 1992 | exact | MATCH |
| Supp S5 | GA-matched Prater Pearson r | 0.808 | 0.808 | 0.002 | MATCH |
| Supp S5 | GA-matched Prater CCC | 0.755 | 0.755 | 0.002 | MATCH |
| Supp S6 | Covariate sweep DEGs (categorical) | 436 | 436 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (linear) | 530 | 530 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (poly2) | 486 | 486 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (ns3) | 513 | 513 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (categorical) | 9 | 9 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (linear) | 11 | 11 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (poly2) | 8 | 8 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (ns3) | 10 | 10 | exact | MATCH |
| Supp S6 | Within-size Jaccard N=10 (categorical) | 0.453 | 0.453 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (linear) | 0.481 | 0.481 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (poly2) | 0.429 | 0.429 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (ns3) | 0.443 | 0.443 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (categorical) | 0.625 | 0.625 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (linear) | 0.647 | 0.647 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (poly2) | 0.565 | 0.565 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (ns3) | 0.573 | 0.573 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (categorical) | 0.967 | 0.967 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (linear) | 0.972 | 0.972 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (poly2) | 0.961 | 0.961 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (ns3) | 0.957 | 0.957 | 0.001 | MATCH |
| Supp S6 | Linear GA most reproducible at every subsample size | TRUE | TRUE | exact | MATCH |
| Section 3.3 | Gained DEGs carrying the ComBat reference batch | 245 | 245 | exact | MATCH |
| Section 3.3 | Gained DEGs lacking the ComBat reference batch | 8 | 8 | exact | MATCH |
| Holdout | Block masking: masked genes per run | 5157 | 5157 | 0.5 | MATCH |
| Holdout | Block masking: true DEGs lost (%) | 24.8 | 24.8 | 0.05 | MATCH |
| Holdout | Block masking: non-DEGs gained (%) | 0.21 | 0.21 | 0.005 | MATCH |
| Holdout | Block masking: genes gaining DEG status per run | 10 | 10 | 0.5 | MATCH |
| Holdout | Block masking, GSE100051 visible: lost (%) | 6.4 | 6.4 | 0.05 | MATCH |
| Holdout | Block masking, GSE100051 hidden: lost (%) | 87.1 | 87.1 | 0.05 | MATCH |
| Holdout | Block masking, GSE100051 hidden: true DEGs lost per run | 36 | 36 | 0.5 | MATCH |
| Holdout | Block masking, GSE100051 hidden: |logFC| ratio | 0.48 | 0.48 | 0.005 | MATCH |
| Holdout | Random-cell masking: lost (%) | 4.2 | 4.2 | 0.05 | MATCH |
| Holdout | Random-cell masking: gained (%) | 0.16 | 0.16 | 0.005 | MATCH |
| Holdout | Worst-case masking: lost (%) | 63.9 | 63.9 | 0.05 | MATCH |
| Holdout | Worst-case masking, GSE100051 hidden: lost (%) | 96.4 | 96.4 | 0.05 | MATCH |
| Holdout | Hidden GSE100051 (ComBat-ref): lost (%) | 81 | 81 | 0.5 | MATCH |
| Holdout | Hidden GSE100051 (ComBat-ref): |logFC| ratio | 0.52 | 0.52 | 0.005 | MATCH |
| Holdout | Hidden GSE100051 (ComBat-ref): oracle sensitivity (%) | 95 | 95 | 0.5 | MATCH |
| Holdout | Hidden GSE37901 (ComBat-ref): lost (%) | 23 | 23 | 0.5 | MATCH |
| Holdout | Hidden GSE37901 (ComBat-ref): |logFC| ratio | 0.88 | 0.88 | 0.005 | MATCH |
| Holdout | Hidden GSE100051 (plain ComBat): lost (%) | 78 | 78 | 0.5 | MATCH |
| Holdout | Hidden GSE100051 (plain ComBat): |logFC| ratio | 0.52 | 0.52 | 0.005 | MATCH |
| Holdout | Imputed cells at weight 1: full-run DEGs | 527 | 527 | exact | MATCH |
| Holdout | Imputed cells at weight 1: DEGs shared with weight 0 | 505 | 505 | exact | MATCH |
| Holdout | Imputed cells at weight 1: FDR-only genes | 9109 | 9109 | exact | MATCH |
| Holdout | Imputed cells at weight 0: FDR-only genes | 6794 | 6794 | exact | MATCH |
| Holdout | Plain ComBat full-run DEGs | 525 | 525 | exact | MATCH |
| Holdout | Article number audit (verify_holdout_claims.R): claims needing attention | 0 | 0 | exact | MATCH |

## Reproduction Outcome

All verified quantitative claims in the article are reproduced by
the companion repository pipeline.

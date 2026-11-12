# Replication Report: Companion Repository vs Article Claims

**Generated:** 2026-11-12
**R version:** R version 4.5.0 (2025-04-11)
**Platform:** x86_64-pc-linux-gnu

## Environment

**renv R version:** 4.5.0

## Summary

- **Total claims verified:** 126
- **Matches:** 126
- **Mismatches:** 0
- **Not verified:** 0

## Comparison Table

| Section | Claim | Article | Reproduced | Tolerance | Status |
|---------|-------|---------|------------|-----------|--------|
| Table 1 | Intersection gene count (k=6) | 8260 | 8260 | exact | MATCH |
| Table 1 | Union gene count (k=1) | 17531 | 17531 | exact | MATCH |
| Table 1 | Safety-drop genes removed | 298 | 298 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref total DEGs | 541 | 541 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref up | 448 | 448 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref down | 93 | 93 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) total DEGs | 277 | 277 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) up | 224 | 224 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) down | 53 | 53 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) total DEGs | 508 | 508 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) up | 425 | 425 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) down | 83 | 83 | exact | MATCH |
| Supp S1 | softImpute random-cell r | 0.978 | 0.978 | 0.0005 | MATCH |
| Supp S1 | softImpute gene-dataset block r | 0.820 | 0.82 | 0.0005 | MATCH |
| Supp S1 | softImpute block RMSE | 1.50 | 1.5 | 0.005 | MATCH |
| Supp S1 | softImpute block MAE | 1.15 | 1.15 | 0.005 | MATCH |
| Supp S1 | softImpute worst-case block r | 0.824 | 0.824 | 0.0005 | MATCH |
| Supp S1 | KNN (centred) block r | 0.815 | 0.815 | 0.0005 | MATCH |
| Supp S1 | biScale fill block RMSE | 1.56 | 1.56 | 0.005 | MATCH |
| Supp S1 | imputePCA worst-case r | 0.842 | 0.842 | 0.0005 | MATCH |
| Section 3.4 | Jaccard at N=10 | 0.54 | 0.54 | 0.01 | MATCH |
| Section 3.4 | Split-half median Jaccard between halves | 0.453 | 0.453 | 0.01 | MATCH |
| Supp S3 | Split-half DEG retention | 86.4 | 86.4 | 0.1 | MATCH |
| Supp S3 | Split-half logFC CCC between halves | 0.771 | 0.771 | 0.01 | MATCH |
| Supp S3 | Subsampling DEG retention at N=10 | 90.9 | 90.9 | 0.1 | MATCH |
| Supp S3 | Subsampling DEG retention at N=30 | 93.0 | 93 | 0.1 | MATCH |
| Supp S3 | Subsampling logFC CCC at N=10 | 0.89 | 0.89 | 0.01 | MATCH |
| Supp S3 | Subsampling logFC CCC at N=100 | 1.00 | 1 | 0.01 | MATCH |
| Section 3.4 | Jaccard at N=100 | 0.98 | 0.98 | 0.01 | MATCH |
| Section 3.4 | Balanced reference DEG count | 484 | 484 | exact | MATCH |
| Section 3.4 | Balanced vs full Jaccard | 0.840 | 0.84 | 0.001 | MATCH |
| Section 3.4 | Balanced retention of full DEGs | 86.5% | 86.5% | exact | MATCH |
| Results 2.3 | Shared genes compared | 8260 | 8260 | exact | MATCH |
| Results 2.3 | Post-ComBat Pearson r (all cells) | 0.9994 | 0.9994 | 0.0005 | MATCH |
| Results 2.3 | Mean abs diff | 0.031 | 0.031 | 0.005 | MATCH |
| Results 2.3 | Median per-gene MAE | 0.027 | 0.027 | 0.005 | MATCH |
| Enrichment | 541-DEG GO BP terms (q<0.05) | 678 | 678 | exact | MATCH |
| Enrichment | 541-DEG KEGG pathways (q<0.05) | 44 | 44 | exact | MATCH |
| Enrichment | 277-DEG GO BP terms | 292 | 292 | exact | MATCH |
| Enrichment | 277-DEG KEGG pathways | 26 | 26 | exact | MATCH |
| Enrichment | 265-DEG GO BP terms | 243 | 243 | exact | MATCH |
| Enrichment | 265-DEG KEGG pathways | 25 | 25 | exact | MATCH |
| Discussion 3.2 | Gained DEGs (541 set minus 277 set) | 265 | 265 | exact | MATCH |
| Discussion 3.2 | Gained genes present in Lykhenko 2021 limma table | 245 | 245 | exact | MATCH |
| Discussion 3.2 | Same direction of change (%) | 89.0 | 89 | 0.05 | MATCH |
| Discussion 3.2 | logFC Pearson r vs Lykhenko 2021 | 0.61 | 0.608 | 0.0051 | MATCH |
| Discussion 3.2 | Already FDR-significant in Lykhenko 2021 (%) | 60.8 | 60.8 | 0.05 | MATCH |
| Supp S4 | Direction-discordant gained genes | 27 | 27 | exact | MATCH |
| Supp S4 | Discordant genes FDR-significant in 2021 | 0 | 0 | exact | MATCH |
| Supp S4 | Median 2021 |logFC| of discordant genes | 0.13 | 0.13 | 0.005 | MATCH |
| Supp S4 | Median 2021 |logFC| of concordant genes | 0.56 | 0.56 | 0.005 | MATCH |
| Supp S4 | 2021-significant matched genes agreeing in direction | 149/149 | 149/149 | exact | MATCH |
| Section 3.3 | GO BP terms (full 541) | 678 | 678 | exact | MATCH |
| Section 3.3 | KEGG pathways (full 541) | 44 | 44 | exact | MATCH |
| Section 3.3 | GO BP terms (intersection 277) | 292 | 292 | exact | MATCH |
| Section 3.3 | KEGG pathways (intersection 277) | 26 | 26 | exact | MATCH |
| Section 3.3 | GO BP terms (gained 265) | 243 | 243 | exact | MATCH |
| Section 3.3 | KEGG pathways (gained 265) | 25 | 25 | exact | MATCH |
| Section 3.3 | q: positive regulation of cytokine production (full) | 9.5e-15 | 9.5e-15 | 0.000000000000001 | MATCH |
| Section 3.3 | q: leukocyte migration (full) | 3.3e-12 | 3.3e-12 | 0.000000000001 | MATCH |
| Section 3.3 | q: chemotaxis (full) | 2e-11 | 2e-11 | 0.000000000001 | MATCH |
| Section 3.3 | q: humoral immune response (full) | 1e-10 | 1e-10 | 0.00000000001 | MATCH |
| Section 3.3 | q: hsa05150 S. aureus infection (full) | 3.2e-10 | 3.2e-10 | 0.00000000001 | MATCH |
| Section 3.3 | q: hsa04610 complement and coagulation (full) | 1.5e-09 | 1.5e-09 | 0.0000000001 | MATCH |
| Section 3.3 | q: hsa04514 cell adhesion molecules (full) | 6.3e-07 | 6.3e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: chemotaxis (gained) | 1.4e-07 | 1.4e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: leukocyte migration (gained) | 1.4e-07 | 1.4e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: hsa04145 phagocytosis (gained) | 0.0034 | 0.0034 | 0.0001 | MATCH |
| Section 3.4 | Gained GO BP terms shared with intersection | 85 | 85 | exact | MATCH |
| Section 3.4 | Gained KEGG pathways shared with intersection | 9 | 9 | exact | MATCH |
| Section 3.3 | Gained DEGs that are intersection genes | 28 | 28 | exact | MATCH |
| Section 3.3 | Gained DEGs testable only via imputation | 237 | 237 | exact | MATCH |
| Supp S6 | qPCR benchmark genes recovered as DEGs | 11 | 11 | exact | MATCH |
| Supp S5 | GA-matched testable genes | 11761 | 11761 | exact | MATCH |
| Supp S5 | GA-matched DEGs | 603 | 603 | exact | MATCH |
| Supp S5 | GA-matched shared genes with Prater | 1992 | 1992 | exact | MATCH |
| Supp S5 | GA-matched Prater Pearson r | 0.808 | 0.808 | 0.002 | MATCH |
| Supp S5 | GA-matched Prater CCC | 0.755 | 0.755 | 0.002 | MATCH |
| Supp S6 | Covariate sweep DEGs (categorical) | 444 | 444 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (linear) | 541 | 541 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (poly2) | 474 | 474 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (ns3) | 511 | 511 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (categorical) | 9 | 9 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (linear) | 11 | 11 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (poly2) | 8 | 8 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (ns3) | 10 | 10 | exact | MATCH |
| Supp S6 | Within-size Jaccard N=10 (categorical) | 0.452 | 0.452 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (linear) | 0.483 | 0.483 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (poly2) | 0.436 | 0.436 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (ns3) | 0.45 | 0.45 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (categorical) | 0.638 | 0.638 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (linear) | 0.662 | 0.662 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (poly2) | 0.578 | 0.578 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (ns3) | 0.593 | 0.593 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (categorical) | 0.969 | 0.969 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (linear) | 0.976 | 0.976 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (poly2) | 0.959 | 0.959 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (ns3) | 0.95 | 0.95 | 0.001 | MATCH |
| Supp S6 | Linear GA most reproducible at every subsample size | TRUE | TRUE | exact | MATCH |
| Section 3.3 | Gained DEGs carrying the ComBat reference batch | 265 | 265 | exact | MATCH |
| Section 3.3 | Gained DEGs lacking the ComBat reference batch | 0 | 0 | exact | MATCH |
| Holdout | Block masking: masked genes per run | 5157 | 5157 | 0.5 | MATCH |
| Holdout | Block masking: true DEGs lost (%) | 26.6 | 26.6 | 0.05 | MATCH |
| Holdout | Block masking: non-DEGs gained (%) | 0.13 | 0.13 | 0.005 | MATCH |
| Holdout | Block masking: genes gaining DEG status per run | 7 | 7 | 0.5 | MATCH |
| Holdout | Block masking, GSE100051 visible: lost (%) | 7.0 | 7 | 0.05 | MATCH |
| Holdout | Block masking, GSE100051 hidden: lost (%) | 96.8 | 96.8 | 0.05 | MATCH |
| Holdout | Block masking, GSE100051 hidden: true DEGs lost per run | 40 | 40 | 0.5 | MATCH |
| Holdout | Block masking, GSE100051 hidden: |logFC| ratio | 0.28 | 0.28 | 0.005 | MATCH |
| Holdout | Random-cell masking: lost (%) | 7.6 | 7.6 | 0.05 | MATCH |
| Holdout | Random-cell masking: gained (%) | 0.13 | 0.13 | 0.005 | MATCH |
| Holdout | Worst-case masking: lost (%) | 66.7 | 66.7 | 0.05 | MATCH |
| Holdout | Worst-case masking, GSE100051 hidden: lost (%) | 100.0 | 100 | 0.05 | MATCH |
| Holdout | Hidden GSE100051 (ComBat-ref): lost (%) | 95 | 95 | 0.5 | MATCH |
| Holdout | Hidden GSE100051 (ComBat-ref): |logFC| ratio | 0.33 | 0.33 | 0.005 | MATCH |
| Holdout | Hidden GSE100051 (ComBat-ref): oracle sensitivity (%) | 95 | 95 | 0.5 | MATCH |
| Holdout | Hidden GSE37901 (ComBat-ref): lost (%) | 18 | 18 | 0.5 | MATCH |
| Holdout | Hidden GSE37901 (ComBat-ref): |logFC| ratio | 0.89 | 0.89 | 0.005 | MATCH |
| Holdout | Hidden GSE100051 (plain ComBat): lost (%) | 93 | 93 | 0.5 | MATCH |
| Holdout | Hidden GSE100051 (plain ComBat): |logFC| ratio | 0.38 | 0.38 | 0.005 | MATCH |
| Holdout | Imputed cells at weight 1: full-run DEGs | 540 | 540 | exact | MATCH |
| Holdout | Imputed cells at weight 1: DEGs shared with weight 0 | 525 | 525 | exact | MATCH |
| Holdout | Imputed cells at weight 1: FDR-only genes | 10664 | 10664 | exact | MATCH |
| Holdout | Imputed cells at weight 0: FDR-only genes | 7779 | 7779 | exact | MATCH |
| Holdout | Plain ComBat full-run DEGs | 508 | 508 | exact | MATCH |
| Holdout | Article number audit (verify_holdout_claims.R): claims needing attention | 0 | 0 | exact | MATCH |

## Reproduction Outcome

All verified quantitative claims in the article are reproduced by
the companion repository pipeline.

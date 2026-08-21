# Replication Report: Companion Repository vs Article Claims

**Generated:** 2026-08-20
**R version:** R version 4.5.0 (2025-04-11)
**Platform:** x86_64-pc-linux-gnu

## Environment

**renv R version:** 4.5.0

## Summary

- **Total claims verified:** 91
- **Matches:** 91
- **Mismatches:** 0
- **Not verified:** 0

## Comparison Table

| Section | Claim | Article | Reproduced | Tolerance | Status |
|---------|-------|---------|------------|-----------|--------|
| Table 1 | Intersection gene count (k=6) | 8260 | 8260 | exact | MATCH |
| Table 1 | Union gene count (k=1) | 17531 | 17531 | exact | MATCH |
| Table 1 | Safety-drop genes removed | 298 | 298 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref total DEGs | 538 | 538 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref up | 446 | 446 | exact | MATCH |
| Table 4 | SoftImpute+ComBat-ref down | 92 | 92 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) total DEGs | 277 | 277 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) up | 224 | 224 | exact | MATCH |
| Table 4 | Intersection-only (none+ref) down | 53 | 53 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) total DEGs | 504 | 504 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) up | 421 | 421 | exact | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) down | 83 | 83 | exact | MATCH |
| Table 3 | Block-mask softImpute r | 0.816 | 0.816 | 0.001 | MATCH |
| Table 3 | Block-mask RMSE | 1.52 | 1.52 | 0.015 | MATCH |
| Table 3 | Block-mask MAE | 1.17 | 1.16 | 0.015 | MATCH |
| Section 3.4 | Jaccard at N=10 | 0.54 | 0.54 | 0.01 | MATCH |
| Section 3.4 | Split-half median Jaccard between halves | 0.402 | 0.402 | 0.01 | MATCH |
| Supp S3 | Split-half DEG retention | 84.1 | 84.1 | 0.1 | MATCH |
| Supp S3 | Split-half logFC CCC between halves | 0.705 | 0.705 | 0.01 | MATCH |
| Supp S3 | Subsampling DEG retention at N=10 | 90.3 | 90.3 | 0.1 | MATCH |
| Supp S3 | Subsampling DEG retention at N=30 | 91.8 | 91.8 | 0.1 | MATCH |
| Supp S3 | Subsampling logFC CCC at N=10 | 0.88 | 0.88 | 0.01 | MATCH |
| Supp S3 | Subsampling logFC CCC at N=100 | 1.00 | 1 | 0.01 | MATCH |
| Section 3.4 | Jaccard at N=100 | 0.98 | 0.98 | 0.01 | MATCH |
| Section 3.4 | Balanced reference DEG count | 484 | 484 | exact | MATCH |
| Section 3.4 | Balanced vs full Jaccard | 0.841 | 0.841 | 0.001 | MATCH |
| Section 3.4 | Balanced retention of full DEGs | 86.8% | 86.8% | exact | MATCH |
| Results 2.3 | Shared genes compared | 8260 | 8260 | exact | MATCH |
| Results 2.3 | Post-ComBat Pearson r (all cells) | 0.9994 | 0.9994 | 0.0005 | MATCH |
| Results 2.3 | Mean abs diff | 0.031 | 0.031 | 0.005 | MATCH |
| Results 2.3 | Median per-gene MAE | 0.027 | 0.027 | 0.005 | MATCH |
| Enrichment | 538-DEG GO BP terms (q<0.05) | 678 | 678 | exact | MATCH |
| Enrichment | 538-DEG KEGG pathways (q<0.05) | 43 | 43 | exact | MATCH |
| Enrichment | 277-DEG GO BP terms | 292 | 292 | exact | MATCH |
| Enrichment | 277-DEG KEGG pathways | 25 | 25 | exact | MATCH |
| Enrichment | 262-DEG GO BP terms | 240 | 240 | exact | MATCH |
| Enrichment | 262-DEG KEGG pathways | 24 | 24 | exact | MATCH |
| Discussion 3.2 | Gained DEGs (538 set minus 277 set) | 262 | 262 | exact | MATCH |
| Discussion 3.2 | Gained genes present in Lykhenko 2021 limma table | 242 | 242 | exact | MATCH |
| Discussion 3.2 | Same direction of change (%) | 88.8 | 88.8 | 0.05 | MATCH |
| Discussion 3.2 | logFC Pearson r vs Lykhenko 2021 | 0.600 | 0.6 | 0.0051 | MATCH |
| Discussion 3.2 | Already FDR-significant in Lykhenko 2021 (%) | 61.2 | 61.2 | 0.05 | MATCH |
| Section 3.3 | GO BP terms (full 538) | 678 | 678 | exact | MATCH |
| Section 3.3 | KEGG pathways (full 538) | 43 | 43 | exact | MATCH |
| Section 3.3 | GO BP terms (intersection 277) | 292 | 292 | exact | MATCH |
| Section 3.3 | KEGG pathways (intersection 277) | 25 | 25 | exact | MATCH |
| Section 3.3 | GO BP terms (gained 262) | 240 | 240 | exact | MATCH |
| Section 3.3 | KEGG pathways (gained 262) | 24 | 24 | exact | MATCH |
| Section 3.3 | q: positive regulation of cytokine production (full) | 7.9e-15 | 7.9e-15 | 0.000000000000001 | MATCH |
| Section 3.3 | q: leukocyte migration (full) | 1.3e-11 | 1.3e-11 | 0.000000000001 | MATCH |
| Section 3.3 | q: chemotaxis (full) | 6.9e-11 | 6.9e-11 | 0.000000000001 | MATCH |
| Section 3.3 | q: humoral immune response (full) | 1.1e-10 | 1.1e-10 | 0.00000000001 | MATCH |
| Section 3.3 | q: hsa05150 S. aureus infection (full) | 3.2e-10 | 3.2e-10 | 0.00000000001 | MATCH |
| Section 3.3 | q: hsa04610 complement and coagulation (full) | 1.5e-09 | 1.5e-09 | 0.0000000001 | MATCH |
| Section 3.3 | q: hsa04514 cell adhesion molecules (full) | 3.1e-06 | 3.1e-06 | 0.0000001 | MATCH |
| Section 3.3 | q: chemotaxis (gained) | 6.5e-07 | 6.5e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: leukocyte migration (gained) | 6.5e-07 | 6.5e-07 | 0.0000001 | MATCH |
| Section 3.3 | q: hsa04145 phagocytosis (gained) | 3.2e-03 | 0.0032 | 0.0001 | MATCH |
| Section 3.4 | Gained GO BP terms shared with intersection | 83 | 83 | exact | MATCH |
| Section 3.4 | Gained KEGG pathways shared with intersection | 8 | 8 | exact | MATCH |
| Section 3.3 | Gained DEGs that are intersection genes | 28 | 28 | exact | MATCH |
| Section 3.3 | Gained DEGs testable only via imputation | 234 | 234 | exact | MATCH |
| Supp S6 | qPCR benchmark genes recovered as DEGs | 11 | 11 | exact | MATCH |
| Supp S5 | GA-matched testable genes | 11761 | 11761 | exact | MATCH |
| Supp S5 | GA-matched DEGs | 603 | 603 | exact | MATCH |
| Supp S5 | GA-matched shared genes with Prater | 1992 | 1992 | exact | MATCH |
| Supp S5 | GA-matched Prater Pearson r | 0.808 | 0.808 | 0.002 | MATCH |
| Supp S5 | GA-matched Prater CCC | 0.755 | 0.755 | 0.002 | MATCH |
| Supp S6 | Covariate sweep DEGs (categorical) | 447 | 447 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (linear) | 538 | 538 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (poly2) | 472 | 472 | exact | MATCH |
| Supp S6 | Covariate sweep DEGs (ns3) | 518 | 518 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (categorical) | 9 | 9 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (linear) | 11 | 11 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (poly2) | 8 | 8 | exact | MATCH |
| Supp S6 | Covariate sweep qPCR recovered (ns3) | 10 | 10 | exact | MATCH |
| Supp S6 | Within-size Jaccard N=10 (categorical) | 0.453 | 0.453 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (linear) | 0.484 | 0.484 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (poly2) | 0.437 | 0.437 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=10 (ns3) | 0.451 | 0.451 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (categorical) | 0.639 | 0.639 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (linear) | 0.662 | 0.662 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (poly2) | 0.579 | 0.579 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=30 (ns3) | 0.594 | 0.594 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (categorical) | 0.97 | 0.97 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (linear) | 0.978 | 0.978 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (poly2) | 0.959 | 0.959 | 0.001 | MATCH |
| Supp S6 | Within-size Jaccard N=100 (ns3) | 0.951 | 0.951 | 0.001 | MATCH |
| Supp S6 | Linear GA most reproducible at every subsample size | TRUE | TRUE | exact | MATCH |
| Section 3.3 | Gained DEGs carrying the ComBat reference batch | 262 | 262 | exact | MATCH |
| Section 3.3 | Gained DEGs lacking the ComBat reference batch | 0 | 0 | exact | MATCH |

## Reproduction Outcome

All verified quantitative claims in the article are reproduced by
the companion repository pipeline.

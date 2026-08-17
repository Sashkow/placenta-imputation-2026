# Replication Report: Companion Repository vs Article Claims

**Generated:** 2026-08-17
**R version:** R version 4.5.0 (2025-04-11)
**Platform:** x86_64-pc-linux-gnu

## Environment

**renv R version:** 4.5.0

## Summary

- **Total claims verified:** 35
- **Matches:** 35
- **Mismatches:** 0
- **Not verified:** 0

## Comparison Table

| Section | Claim | Article | Reproduced | Status |
|---------|-------|---------|------------|--------|
| Table 1 | Intersection gene count (k=6) | 8260 | 8260 | MATCH |
| Table 1 | Union gene count (k=1) | 17531 | 17531 | MATCH |
| Table 1 | Safety-drop genes removed | 298 | 298 | MATCH |
| Table 4 | SoftImpute+ComBat-ref total DEGs | 538 | 538 | MATCH |
| Table 4 | SoftImpute+ComBat-ref up | 446 | 446 | MATCH |
| Table 4 | SoftImpute+ComBat-ref down | 92 | 92 | MATCH |
| Table 4 | Intersection-only (none+ref) total DEGs | 277 | 277 | MATCH |
| Table 4 | Intersection-only (none+ref) up | 224 | 224 | MATCH |
| Table 4 | Intersection-only (none+ref) down | 53 | 53 | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) total DEGs | 504 | 504 | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) up | 421 | 421 | MATCH |
| Table 4 | SoftImpute+ComBat (no ref) down | 83 | 83 | MATCH |
| Table 3 | Block-mask softImpute r | 0.816 | 0.816 | MATCH |
| Table 3 | Block-mask RMSE | 1.52 | 1.52 | MATCH |
| Table 3 | Block-mask MAE | 1.17 | 1.16 | MATCH |
| Section 3.4 | Jaccard at N=10 | 0.54 | 0.54 | MATCH |
| Section 3.4 | Split-half median Jaccard between halves | 0.45 | 0.4 | MATCH |
| Section 3.4 | Balanced reference DEG count | 484 | 484 | MATCH |
| Section 3.4 | Balanced vs full Jaccard | 0.841 | 0.841 | MATCH |
| Section 3.4 | Balanced retention of full DEGs | 86.8% | 86.8% | MATCH |
| Results 2.3 | Shared genes compared | 8260 | 8260 | MATCH |
| Results 2.3 | Post-ComBat Pearson r (all cells) | 0.9994 | 0.9994 | MATCH |
| Results 2.3 | Mean abs diff | 0.031 | 0.031 | MATCH |
| Results 2.3 | Median per-gene MAE | 0.027 | 0.027 | MATCH |
| Enrichment | 538-DEG GO BP terms (q<0.05) | 678 | 678 | MATCH |
| Enrichment | 538-DEG KEGG pathways (q<0.05) | 43 | 43 | MATCH |
| Enrichment | 277-DEG GO BP terms | 292 | 292 | MATCH |
| Enrichment | 277-DEG KEGG pathways | 25 | 25 | MATCH |
| Enrichment | 262-DEG GO BP terms | 240 | 240 | MATCH |
| Enrichment | 262-DEG KEGG pathways | 24 | 24 | MATCH |
| Discussion 3.2 | Gained DEGs (538 set minus 277 set) | 262 | 262 | MATCH |
| Discussion 3.2 | Gained genes present in Lykhenko 2021 limma table | 242 | 242 | MATCH |
| Discussion 3.2 | Same direction of change (%) | 88.8 | 88.8 | MATCH |
| Discussion 3.2 | logFC Pearson r vs Lykhenko 2021 | 0.600 | 0.6 | MATCH |
| Discussion 3.2 | Already FDR-significant in Lykhenko 2021 (%) | 61.2 | 61.2 | MATCH |

## Reproduction Outcome

All verified quantitative claims in the article are reproduced by
the companion repository pipeline.

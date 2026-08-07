# Replication Report: Companion Repository vs Article Claims

**Generated:** 2026-08-07
**R version:** R version 4.5.0 (2025-04-11)
**Platform:** x86_64-pc-linux-gnu

## Environment

**renv R version:** 4.5.0

## Summary

- **Total claims verified:** 24
- **Matches:** 23
- **Mismatches:** 0
- **Not verified:** 1

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
| Section 3.2 | Sensitivity Pearson r (8260 genes) | 0.9994 | 0.9991 | MATCH |
| Section 3.2 | Sensitivity mean abs diff | 0.031 | 0.033 | MATCH |
| Section 3.2 | Sensitivity median per-gene MAE | 0.027 | 0.024 | MATCH |
| Enrichment | All enrichment counts | — | — | NOT VERIFIED (enrichment dir missing) |

## Reproduction Outcome

All 23 verified claims match the article. 1 claim(s) could not be verified due to missing data.

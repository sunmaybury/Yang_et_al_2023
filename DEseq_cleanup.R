
## ========================================================
## 22.02 WT and I-SceI fibroblasts H3K27ac CUT&Tag
## data processed and DEseq2 analysis done on O2
## R version 3.6.1 on O2
## CNT_DEseq2.R run on O2
## Cleaning up results locally for gene assignment etc
## ========================================================

setwd("~/Dropbox (HMS)/SunMayburyLewis/_Jae_ICE_paper/22.02/Nextseq/DEseq2/")


## =========================================
## Data cleanup and write out files
## =========================================

## Load data
res1 <- read.csv("results_WT_ISceI_H3K27ac.csv")
## remove first column
res1[1] <- NULL

## Reorder by padj
res2 <- res1[order(res1$padj), ]
## Filter by padj < 0.05
res2 <- res2[res2$padj < 0.05, ]

## Clean up data frame
res2clean <- res2[, c(1,2,3,12,13,14,15,16,17)]

# write.table(res2, file="results_WT_ISceI_H3K27ac_filtered.csv", sep=",", row.names=F, col.names=T, quote=F)
# write.table(res2clean, file="results_WT_ISceI_H3K27ac_clean.csv", sep=",", row.names=F, col.names=T, quote=F)

## Filter further by abs[log2FC] > 1.5
res3 <- res2[abs(res2$log2FoldChange) >=1.5, ]
res3clean <- res3[, c(1,2,3,12,13,14,15,16,17)]

# write.table(res3, file="results_WT_ISceI_H3K27ac_filtered.csv", sep=",", row.names=F, col.names=T, quote=F)
# write.table(res3clean, file="results_WT_ISceI_H3K27ac_clean.csv", sep=",", row.names=F, col.names=T, quote=F)

## Also get txt files for gene assignment
# write.table(res3clean, file="results_WT_ISceI_H3K27ac_clean.txt", sep="\t", row.names=F, col.names=F, quote=F)

## ======================================================
## Separate increased and decreased H3K27ac peak list
## Want to run GREAT separately on these two lists
## ======================================================

incre <- res3clean[res3clean$log2FoldChange>0, ]
decre <- res3clean[res3clean$log2FoldChange<0, ]

write.table(incre, file="results_WT_ISceI_H3K27ac_increasing.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(decre, file="results_WT_ISceI_H3K27ac_decreasing.txt", sep="\t", row.names=F, col.names=T, quote=F)



## =============================================================================
## Write out res2clean without filtering by padj, without filtering by log2FC
## For Jae to run GO analysis with various cutoffs
## =============================================================================

incre <- res2clean[res2clean$log2FoldChange>0, ]
decre <- res2clean[res2clean$log2FoldChange<0, ]

write.table(incre, file="results_WT_ISceI_H3K27ac_increasing_nofilt.txt", sep="\t", row.names=F, col.names=T, quote=F)
write.table(decre, file="results_WT_ISceI_H3K27ac_decreasing_nofilt.txt", sep="\t", row.names=F, col.names=T, quote=F)

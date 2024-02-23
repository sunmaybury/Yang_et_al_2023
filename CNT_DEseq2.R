## ==========================================================
## Experiment 22.02 H3K27ac CUT&Tag for WT and I-SceI cells
## For Jae's Cell manuscript (ICE) revision
## Experiment done, data processed
## ==========================================================

library('dplyr')
library('readr')
library('stringr')
library('ggplot2')
library('viridis')
library('GenomicRanges')
library('chromVAR') ## For FRiP analysis and differential analysis
library('DESeq2') ## For differential analysis section
library('Rsubread')
library('ggrepel')
library('pheatmap')
library('RColorBrewer')

setwd("/n/data2/hms/genetics/sinclair/Sun/Jae_H3K27ac_CNT/")

## ===================================================
## Differential analysis
## ===================================================

sampleList <- c("A3_S1_pr1","A3_S1_pr2","A5_S5","A6_S6")
exp="H3K27ac"

## Create the peak x sample matrix
## Create a master peak list merging all peaks called for each sample
mPeak <- GRanges()
## overlap with bam file to get count
for(samp in sampleList){
  peakRes <- read.table(paste0("peakCalling/SEACR-stringent/", samp, "_seacr_IgG.peaks.stringent.bed"), header=FALSE, fill=TRUE)
  mPeak <- GRanges(seqnames=peakRes$V1, IRanges(start=peakRes$V2, end=peakRes$V3), strand="*") %>% append(mPeak, .)
}
masterPeak <- reduce(mPeak)

## Get the seq data separately from masterPeak list
## This should be appended to counts and results csv files later
masterSeq <- as.data.frame(granges(masterPeak))
# write.table(masterSeq, file="WT_ISceI_H3K27ac_master_peaklist.txt", sep="\t", quote=F, row.names=F, col.names=F)
masterSeq <- masterSeq[, -c(4,5)]

## Get the fragment counts for each peak in the master peak list
countMat <- matrix(NA, length(masterPeak), length(sampleList))
i=1
for(samp in sampleList){
  bamFile <- paste0("mapped/bam/", samp, ".mapped.bam")
  fragment_counts <- getCounts(bamFile, masterPeak, paired=TRUE, by_rg=FALSE, format="bam")
  countMat[, i] <- counts(fragment_counts)[,1]
  i=i+1
}
colnames(countMat) <- sampleList
head(countMat)

cts <- (countMat)
## Appending seq info for csv files
ctsApp <- cbind(masterSeq, countMat)
## writing raw feature counts (not normalized)
# write.csv(ctsApp, file=paste("raw-featurecounts_DESEQ2_", exp, ".csv", sep=""))

## if need to read feature counts after they are written
#cts <- read.csv(paste("raw-featurecounts_DESEQ2_",exp,".csv", sep="")) 
#rownames(cts) <- cts[, 1] ## set rownames
#cts <- cts[, -1]          ## remove the first variable

group <- c("WT","WT","ISceI","ISceI")
coldata <- data.frame(sampleList, group)

## remove low count genes
selectR <- which(rowSums(cts) > 3)
dataS <- cts[selectR,]
dds <- DESeqDataSetFromMatrix(countData = dataS,
                             colData = coldata,
                             design = ~ group)
## DEseq analysis
DDS <- DESeq(dds)

## normalization with respect to the sequencing depth
normDDS <- counts(DDS, normalized = TRUE)
colnames(normDDS) <- paste0(colnames(normDDS), "_norm")

## write normalized feature counts
normDDSApp <- as.data.frame(normDDS)
normDDSApp <- cbind(masterSeq, normDDSApp)

# write.csv(normDDSApp, file="norm_feature_counts_H3K27ac.csv") 

## ===================================================
## Results and plots
## ===================================================

## Plot dispersion estimate
plotDispEsts(DDS)

## WT vs I-SceI
res <- results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs", contrast=c("group","ISceI","WT"))
countMatDiff <- cbind(masterSeq, dataS, normDDS, res)
head(countMatDiff)
# write.csv(as.data.frame(countMatDiff), file="results_WT_ISceI_H3K27ac.csv")

## MAplot
jpeg("MA_plot_WT_ISceI_H3K27ac.jpeg")
plotMA(res)
dev.off()

## Summary
summary <- summary(res)
# out of 138520 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 7424, 5.4%
# LFC < 0 (down)     : 5726, 4.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

## PCA
## Transform counts for data visualization
# rld <- rlog(dds, blind=FALSE)
# jpeg("PCA_plot_H3K27ac_rld.jpeg")
# plotPCA(vsd, intgroup="group")
# dev.off()

vsd <- vst(dds, blind=FALSE)
p <- plotPCA(vsd, intgroup="group")

pdf.options(reset=T, onefile=T)
pdf(file="PCA_vsd.pdf")
p + geom_text_repel(aes(label=name), color="black", size=4) +
  labs(title="H3K27ac") +
  coord_fixed(ratio=1, ylim=c(-3,3), xlim=c(-8,8)) +
  theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
dev.off()

## Heatmap of sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
          # A3_S1_pr1 A3_S1_pr2    A5_S5    A6_S6
# A3_S1_pr1   0.00000  25.53785 70.80134 63.16305
# A3_S1_pr2  25.53785   0.00000 70.72273 63.13301
# A5_S5      70.80134  70.72273  0.00000 50.67277
# A6_S6      63.16305  63.13301 50.67277  0.00000

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf.options(reset=T, onefile=T)
pdf(file="DistMatrix_H3K27ac.pdf")
pheatmap(sampleDistMatrix,
         cellwidth=20, cellheight=20,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

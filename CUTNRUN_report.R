
## =======================================================
## 22.02 H3K27ac CUT&Tag experiment
## For Jae's WT and I-SceI cells
## Report
## April 2022
## =======================================================

setwd("~/Dropbox (HMS)/SunMayburyLewis/_Jae_ICE_paper/22.02/Nextseq/")

library("dplyr")
library("readr")
library("stringr")
library("ggplot2")
library("directlabels")
library("viridis")
library("GenomicRanges")
library("chromVAR") ## For FRiP analysis and differential analysis
library("DESeq2") ## For differential analysis section
library("ggpubr") ## For customizing figures
library("corrplot") ## For correlation plot

sampleList <- c("A3_S1","B2_S2","A4_S3","B3_S4","A5_S5","A6_S6","A7_S7")

# ===========================================================
## Make graphs of mapping summaries of mouse and spike in
# ===========================================================

## changed row 6 to row 10 since bowtie2 gave slightly different result file
alignResult <- c()
for(samp in sampleList){
  alignRes <- read.table(paste0("mapped/", samp, ".txt"), header = FALSE, fill = TRUE)
  alignRate <- substr(alignRes$V1[10], 1, nchar(as.character(alignRes$V1[10]))-1)
  alignResult <- data.frame(Sample = samp,
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric,
                           MappedFragNum_mm10 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric,
                           AlignmentRate_mm10 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}

alignResult %>% mutate(AlignmentRate_mm10 = paste0(AlignmentRate_mm10, "%"))

  # Sample SequencingDepth MappedFragNum_mm10 AlignmentRate_mm10
# 1  A3_S1        57188958           53939776             94.56%
# 2  B2_S2           52393               1413              5.22%
# 3  A4_S3        15156926           14339909             94.89%
# 4  B3_S4         6120871            5736544             94.34%
# 5  A5_S5        30074389           28650208             95.42%
# 6  A6_S6        26261132           25050292              95.6%
# 7  A7_S7          358772             319670             89.74%

# ======================================
## Generate sequencing depth boxplot
# ======================================

fig1A <- alignResult %>% ggplot(aes(x = Sample, y = SequencingDepth/1000000)) +
  geom_boxplot() +
  #geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  #scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 80) +
  ylab("Sequencing Depth per Million") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")

# ggsave(fig1A, filename = "Fig1A.pdf", width = 4, height = 4)

fig1B <- alignResult %>% ggplot(aes(x = Sample, y = MappedFragNum_mm10/1000000)) +
  geom_boxplot() +
  #geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  #scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 80) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (mm10)")

# ggsave(fig1B, filename = "Fig1B.pdf", width = 4, height = 4)

fig1C <- alignResult %>% ggplot(aes(x = Sample, y = AlignmentRate_mm10)) +
  geom_boxplot() +
  #geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  #scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 80) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (mm10)")

# ggsave(fig1C, filename = "Fig1C.pdf", width = 4, height = 4)

# ggarrange(fig1A, fig1B, fig1C, ncol=2, nrow=2, common.legend=TRUE, legend="bottom")
# ggsave("Fig1.pdf")

# ===============================================
## graph of fragment lengths
# ===============================================

## Collect the fragment size information
fragLen <- c()

for(samp in sampleList){
  fragLen <- read.table(paste0("mapped/fragmentLen/", samp, ".txt"), header = FALSE) %>% 
  mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Sample = samp) %>% rbind(fragLen, .) 
}
#fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
#fragLen$Histone = factor(fragLen$Histone, levels = histList)
head(fragLen)

## Generate the fragment size density plot (violin plot)
fig2A <- fragLen %>% ggplot(aes(x = Sample, y = fragLen, weight = Weight)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  #scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  #scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 60) +
  ylab("Fragment Length") +
  xlab("") +
  ggtitle("A. Fragment lengths")

# ggsave(fig2A, filename = "Fig2A.pdf", width = 4, height = 4)

fig2B <- fragLen %>% ggplot(aes(x = fragLen, y = fragCount, fill = Sample)) +
  geom_line(aes(color = Sample), size = 0.5) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "turbo") +
  theme_bw(base_size = 12) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 350)) +
  ggtitle("B. Fragment lengths")

# ggsave(fig2B, filename = "Fig2B.pdf", width = 8, height = 6)

# ggarrange(fig2A, fig2B, ncol = 1, nrow = 2, legend = T)
# ggsave("Fig2.pdf")

# ==============================================
## visualisation of peak width etc
# ==============================================

## SEACR setting stringent
sampleList2 <- c("A3_S1_pr1","A3_S1_pr2","A4_B3_merge","A5_S5","A6_S6")

## number of peaks
peakN = c()
peakWidth = c()
for(samp in sampleList2){
    peakInfo = read.table(paste0("peakCalling/SEACR-stringent/", samp, "_seacr_IgG.peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
    peakN = data.frame(peakN = nrow(peakInfo), Sample = samp) %>% rbind(peakN, .)
    peakWidth = data.frame(width = peakInfo$width, Sample = samp)  %>% rbind(peakWidth, .)
  }
peakNumOption1A <- peakN %>% select(Sample, peakN)
peakNumOption1A

        # Sample  peakN
# 1   A3_S1_pr1 103714
# 2   A3_S1_pr2 114176
# 3 A4_B3_merge 105519
# 4       A5_S5  93622
# 5       A6_S6 115895

## SEACR setting stringent, top 1 percent peaks
peakN = c()
peakWidth = c()
for(samp in sampleList2){
  peakInfo = read.table(paste0("peakCalling/SEACR-stringent/", samp, "_seacr_top0.01.peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
  peakN = data.frame(peakN = nrow(peakInfo), Sample = samp) %>% rbind(peakN, .)
  peakWidth = data.frame(width = peakInfo$width, Sample = samp)  %>% rbind(peakWidth, .)
}
peakNumOption1B <- peakN %>% select(Sample, peakN)
peakNumOption1B

        # Sample peakN
# 1   A3_S1_pr1 19740
# 2   A3_S1_pr2 19796
# 3 A4_B3_merge 17967
# 4       A5_S5 18515
# 5       A6_S6 21458



## ==============================================================================
## Hold off for now
## ==============================================================================






## fragment proportion in peaks regions (FRiPs)
inPeakData = c()

## overlap with bam file to get count
for(samp in sampleList2){
  peakRes <- read.table(paste0("peakCalling/SEACR-stringent/", samp, "_seacr_IgG.peaks.stringent.bed"), header = FALSE, fill = TRUE)
  peak.gr <- GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
  bamFile <- paste0("mapped/bam/", samp, ".mapped.bam")
  fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
  inPeakN = counts(fragment_counts)[,1] %>% sum
  inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Sample = samp))
}

frip = left_join(inPeakData, alignResult, by = c("Sample")) %>% mutate(frip = inPeakN/MappedFragNum_mm10 * 100)
frip %>% select(Sample, SequencingDepth, MappedFragNum_mm10, AlignmentRate_mm10, FragInPeakNum = inPeakN, FRiPs = frip)

       
## visualisation

fig4A = peakN %>% ggplot(aes(x = Sample, y = peakN, fill = Sample)) +
  geom_boxplot() +
  #geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  #facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 75) +
  ylab("Number of Peaks") +
  xlab("") +
  ggtitle("A. Number of peaks")

# ggsave(fig4A, filename = "Fig4A.pdf", width = 6, height = 4)

fig4B = peakWidth %>% ggplot(aes(x = Sample, y = width, fill = Sample)) +
  geom_violin() +
  #facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 75) +
  ylab("Width of Peaks") +
  xlab("") +
  ggtitle("B. Width of peaks")

# ggsave(fig4B, filename = "Fig4B.pdf", width = 6, height = 4)

fig4C = frip %>% ggplot(aes(x = Sample, y = frip, fill = Sample, label = round(frip, 2))) +
  geom_boxplot() +
  #geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 12) +
  ggpubr::rotate_x_text(angle = 75) +
  ylab("% of Fragments in Peaks") +
  xlab("") +
  ggtitle("C. FRiP scores")

# ggsave(fig4C, filename = "Fig4C.pdf", width = 6, height = 4)

ggarrange(fig4A, fig4B, fig4C, ncol = 2, nrow=2)
ggsave("Fig4.pdf")

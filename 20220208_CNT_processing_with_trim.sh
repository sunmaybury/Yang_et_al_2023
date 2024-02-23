#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-12:00						 
#SBATCH -p priority
#SBATCH --mem 12G
#SBATCH -o CUTNTAG-2202_%j.out

SAMP=(
A3_S1  ## // split into 2 pseudoreplicates
#### B2_S2 ## // do not use
A4_S3 ## // combine A4_S3 and B3_S4 
B3_S4 
A5_S5
A6_S6
A7_S7)

module load gcc fastqc bowtie2 samtools R/3.6.1 bedtools python/3.7.4 picard/2.8.0 cutadapt java

cd /n/data2/hms/genetics/sinclair/Sun/Jae_H3K27ac_CNT/

# =================================
## Part 1 FASTQC and trimming
# =================================

# cd ./fastq/I-SceI_CUT_and_Tag/
# 
# mkdir -p fastqc
# mkdir -p trimmed

# for i in "${!SAMP[@]}"; do

# echo fastqc ${SAMP[i]}
# fastqc -o fastqc ${SAMP[i]}_R1_001.fastq.gz ${SAMP[i]}_R2_001.fastq.gz

# echo trim ${SAMP[i]}
# cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -j 4 --nextseq-trim=20 -u -50 -U -50 -o ./trimmed/${SAMP[i]}_R1_001.fastq.gz -p ./trimmed/${SAMP[i]}_R2_001.fastq.gz ./${SAMP[i]}_R1_001.fastq.gz ./${SAMP[i]}_R2_001.fastq.gz
# ### Use the command below
# cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -j 4 --nextseq-trim=20 -o ./trimmed/${SAMP[i]}_R1_001.fastq.gz -p ./trimmed/${SAMP[i]}_R2_001.fastq.gz ./${SAMP[i]}_R1_001.fastq.gz ./${SAMP[i]}_R2_001.fastq.gz

# done


# ======================================================================
## PART 2 MAPPING (this takes a long time, approx 15 mins for one file)
# ======================================================================

# mkdir -p mapped
# mkdir -p mapped/fragmentLen
# 
# for i in "${!SAMP[@]}"; do
#  
# ## mm10 alignment
# echo align mm10 reads ${SAMP[i]}
# bowtie2 -p 6 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /home/sum879/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 ./fastq/I-SceI_CUT_and_Tag/trimmed/${SAMP[i]}_R1_001.fastq.gz -2 ./fastq/I-SceI_CUT_and_Tag/trimmed/${SAMP[i]}_R2_001.fastq.gz -S ./mapped/${SAMP[i]}.sam &> ./mapped/${SAMP[i]}.txt
# 
# ## spike in alignment
# echo aligned mm10 reads ${SAMP[i]}
# bowtie2 -p 6 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x /home/sum879/genomes/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index/genome -1 ./fastq/I-SceI_CUT_and_Tag/trimmed/${SAMP[i]}_R1_001.fastq.gz -2 ./fastq/I-SceI_CUT_and_Tag/trimmed/${SAMP[i]}_R2_001.fastq.gz -S ./mapped/${SAMP[i]}.spikeIn.sam &> ./mapped/${SAMP[i]}.spikeIn.txt
# 
# ## Extract fragment length
# echo extracted fragment length ${SAMP[i]}
# samtools view -F 0x04 ./mapped/${SAMP[i]}.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > ./mapped/fragmentLen/${SAMP[i]}.txt
# 
# done

# ======================================================================
## PART 3 FILTER BAM AND SPLIT WT TO PSEUDOREPLICATES AND COMBINE SC1
# ======================================================================

# mkdir -p mapped/bam
# mkdir -p mapped/bed

# for i in "${!SAMP[@]}"; do
# 
# # echo Filter and keep mapped pairs ${SAMP[i]}
# # samtools view -bS -F 0x04 ./mapped/${SAMP[i]}.sam > ./mapped/bam/${SAMP[i]}.mapped.bam
# 
# echo Sort samples by name
# samtools sort -n ./mapped/bam/${SAMP[i]}.mapped.bam > ./mapped/bam/${SAMP[i]}.mapped.sorted.bam
# 
# done


# ## To do:
# ## [1] Split A3_S1 into 2 pseudoreplicates (sample WT3) using random_split_bed.py
# ## [2] Combine A4_S3 and B3_S4 (Sc1)



# ## ================================ [1]
# ## Start from ./mapped/bam/A3_S1.mapped.sorted.bam
# cd mapped/bam
# 
# ## // for paired end bam file split
# 
# ## Count mapped reads
# PAIRS=$(samtools view -c -f 2 A3_S1.mapped.sorted.bam)
# MATES=$((PAIRS /2))
# echo PAIRS IS $PAIRS
# echo MATES is $MATES
# 
# ## If mates is odd, round down to even number
# if [[ $((MATES % 2)) -eq 0 ]];
#    then LINES=$MATES;
#    else LINES=$((MATES - 1));
# fi
# echo LINE is $LINES
# 
# ## Shuffle and split file
# echo Shuffle and Split bam
# samtools view -f 2  A3_S1.mapped.sorted.bam | awk '{printf("%s%s",$0,(NR%2==0)?"\n":"\0")}' | shuf | tr "\0" "\n" | split -d -l "${LINES}" - "A3_S1_pr"
# 
# ## Get header
# echo Get header 
# samtools view -H A3_S1.mapped.sorted.bam > A3_S1.header.sam
# 
# ## Replace header
# echo Replace header
# cat A3_S1.header.sam  A3_S1_pr00 | samtools sort -@ 6 -T "A3_S1_pr1" -o "A3_S1_pr1.mapped.bam" -
# cat A3_S1.header.sam A3_S1_pr01 | samtools sort -@ 6 -T "A3_S1_pr2" -o "A3_S1_pr2.mapped.bam" -

 
# ## ================================ [2]
# ## Combine A4_S3 and B3_S4 by the following command
# samtools merge A4_B3_merge.mapped.bam A4_S3.mapped.sorted.bam B3_S4.mapped.sorted.bam

## Now new sample list:

SAMP=(
A3_S1_pr1
A3_S1_pr2
# B2_S2
A4_B3_merge)
# A5_S5
# A6_S6
# A7_S7)

# ===============================================================
## PART 4 FILE CONVERSION FOR PEAKCALLING ETC
# ===============================================================

# for i in "${!SAMP[@]}"; do
# 
# echo Sort new samples by name ${SAMP[i]}
# samtools sort -n ./mapped/bam/${SAMP[i]}.mapped.bam > ./mapped/bam/${SAMP[i]}.mapped.sorted.bam
# 
# echo Convert to bed ${SAMP[i]}
# bedtools bamtobed -i ./mapped/bam/${SAMP[i]}.mapped.sorted.bam -bedpe > ./mapped/bed/${SAMP[i]}.bed
# 
# echo Keep read pairs on same chromosomes and fragment length less than 1000bp ${SAMP[i]}
# awk '$1==$4 && $6-$2 < 1000 {print $0}' ./mapped/bed/${SAMP[i]}.bed > ./mapped/bed/${SAMP[i]}.clean.bed
# 
# echo Extract fragment related columns only ${SAMP[i]}
# cut -f 1,2,6 ./mapped/bed/${SAMP[i]}.clean.bed | sort -k1,1 -k2,2n -k3,3n  >./mapped/bed/${SAMP[i]}.fragments.bed 
# 
# done

# =============================================================
## PART 5 Creation of bedgraph files 
# =============================================================
  
# chromSize="/home/sum879/genomes/Mus_musculus/mm10.chrom.sizes"

# mkdir -p ./mapped/bedgraph 

# for i in "${!SAMP[@]}"; do
# 
# echo Make bedgraph ${SAMP[i]}
# bedtools genomecov -bg -i ./mapped/bed/${SAMP[i]}.fragments.bed -g $chromSize > ./mapped/bedgraph/${SAMP[i]}.fragments.bedgraph
# 
# done

# ==================================================================
## PART 6 call peaks (use IgG control as reference)
# ==================================================================
 
## For peakcalling, excluding IgG sample A7

SAMP=(
A3_S1_pr1
A3_S1_pr2
A4_B3_merge
A5_S5
A6_S6)

# seacr="/home/sum879/SEACR-master/SEACR_1.3.sh"
# 
# mkdir -p ./peakCalling/SEACR-stringent
# mkdir -p ./peakCalling/SEACR-relaxed

# for i in "${!SAMP[@]}"; do

# echo Call peaks against IgG stringent ${SAMP[i]}
# bash $seacr ./mapped/bedgraph/${SAMP[i]}.fragments.bedgraph ./mapped/bedgraph/A7_S7.fragments.bedgraph norm stringent ./peakCalling/SEACR-stringent/${SAMP[i]}_seacr_IgG.peaks
# 
# echo Call peaks against IgG relaxed ${SAMP[i]}
# bash $seacr ./mapped/bedgraph/${SAMP[i]}.fragments.bedgraph ./mapped/bedgraph/A7_S7.fragments.bedgraph norm relaxed ./peakCalling/SEACR-relaxed/${SAMP[i]}_seacr_IgG.peaks

# echo Call top 1% peaks stringent ${SAMP[i]}
# bash $seacr ./mapped/bedgraph/${SAMP[i]}.fragments.bedgraph 0.01 norm stringent ./peakCalling/SEACR-stringent/${SAMP[i]}_seacr_top0.01.peaks
# 
# echo Call top 1% peaks relaxed ${SAMP[i]}
# bash $seacr ./mapped/bedgraph/${SAMP[i]}.fragments.bedgraph 0.01 norm relaxed ./peakCalling/SEACR-relaxed/${SAMP[i]}_seacr_top0.01.peaks
# 
# done


###### Skip this option since we're not using spike-in calibration for now
# ## Calibrated to spike-in

# bash $seacr ./mapped/bedgraph/${SAMP[i]}.fragments.spikeIn.normalized.bedgraph ./mapped/bedgraph/2104-B2_S5.fragments.spikeIn.normalized.bedgraph norm stringent ./peakCalling/SEACR-stringent/${SAMP[i]}_seacr_control.spikeIn.peaks
# echo called peaks with IgG control ${SAMP[i]} spikeIn

# bash $seacr ./mapped/bedgraph/${SAMP[i]}.fragments.spikeIn.normalized.bedgraph 0.01 stringent ./peakCalling/SEACR-stringent/${SAMP[i]}_seacr_top0.01.spikeIn.peaks
# echo called top 1 percent peaks ${SAMP[i]} spikeIn


# =================================================================
# ## PART 7 make bigwig files to check on igv
# # ===============================================================

SAMP=(
A3_S1_pr1
A3_S1_pr2
# B2_S2
A4_B3_merge
A5_S5
A6_S6
A7_S7)


## Make merged bigwig files for WT and I-SceI

mkdir -p ./mapped/bigwig

for i in "${!SAMP[@]}"; do

echo Sort mapped bam ${SAMP[i]}
samtools sort ./mapped/bam/${SAMP[i]}.mapped.bam -o ./mapped/bam/${SAMP[i]}.sorted.bam
echo Index bam ${SAMP[i]}                                                   
samtools index ./mapped/bam/${SAMP[i]}.sorted.bam

## merge bam files
## For merged WT just use A3_S1 before splitting to pr1 and pr2
# echo Sort merged bam A3_S1
# samtools sort ./mapped/bam/A3_S1.mapped.bam -o ./mapped/bam/A3_S1.sorted.bam
# echo Index merged bam A3_S1
# samtools index ./mapped/bam/A3_S1.sorted.bam
# echo Merge sc1 sc2 samples
# samtools merge ./mapped/bam/I-SceI.merged.bam ./mapped/bam/A5_S5.mapped.bam ./mapped/bam/A6_S6.mapped.bam
# echo Sort merged bam I-SceI
# samtools sort ./mapped/bam/I-SceI.merged.bam -o ./mapped/bam/I-SceI.merged.sorted.bam
# echo Index merged bam I-SceI
# samtools index ./mapped/bam/I-SceI.merged.sorted.bam

## need to be in python virtual environment to use deeptools
# source /home/sum879/python3.7.4/bin/activate

# echo Generate bigWig ${SAMP[i]}
# bamCoverage -b ./mapped/bam/${SAMP[i]}.sorted.bam -o ./mapped/bigwig/${SAMP[i]}_raw.bw
# echo Generate bigwig A3_S1
# bamCoverage -b ./mapped/bam/A3_S1.sorted.bam -o ./mapped/bigwig/A3_S1_raw.bw
# echo Generate bigwig I-SceI merged
# bamCoverage -b ./mapped/bam/I-SceI.merged.sorted.bam -o ./mapped/bigwig/I-SceI_merge_raw.bw
	
## after done with python environment
# deactivate

# done



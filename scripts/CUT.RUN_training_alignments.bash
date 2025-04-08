#!/bin/bash
#
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --verbose
#SBATCH --error=/ix1/yarbely/mam835/training/CR_PDNC4/job.%J.err
#SBATCH --output=/ix1/yarbely/mam835/training/CR_PDNC4/job.%J.out

module load cutadapt/2.10
module load gcc/8.2.0
module load bowtie2/2.4.1
module load samtools/1.14

cd /ix1/yarbely/mam835/training/CR_PDNC4/

#bowtie2 \
#	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 8 \
#	-x /ix1/yarbely/mam835/training/CR_PDNC4/h38.p14 \
#	-1 /bgfs/yarbely/mam835/CUT.RUN/Pool5/C2_12m_sc_D4/C2_12m_sc_D4.trimmed.R1.fastq \
#	-2 /bgfs/yarbely/mam835/CUT.RUN/Pool5/C2_12m_sc_D4/C2_12m_sc_D4.trimmed.R2.fastq \
#	-S E2_12m_sc_D4.sam

# If you are wondering why the fastq files are called C2_12m but I named the output as E2_12m, 
## it's because once it was sequenced, we discovered there was a mixup between C2 12m and E2 12m cells

#bowtie2 \
#	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 8 \
#	-x /ix1/yarbely/mam835/training/CR_PDNC4/h38.p14 \
#	-1 /bgfs/yarbely/mam835/CUT.RUN/Pool2/PDNC4_CA-HJ-LAP_cC4_Y/PDNC4_CA-HJ-LAP_cC4_Y.trimmed.R2.fastq \
#	-2 /bgfs/yarbely/mam835/CUT.RUN/Pool2/PDNC4_CA-HJ-LAP_cC4_Y/PDNC4_CA-HJ-LAP_cC4_Y.trimmed.R1.fastq \
#	-S PDNC4_CA-HJ-LAP_cC4_Y.sam

### For the E.coli alignments, I've used my E.coli genome that I downlaoded and indexed
## Before running, download the E.coli genome and make an index with bowtie2
## I used https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001308125.1/
## Look at the previous training module for information on how to do that

#bowtie2 \
#	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 8 \
#	-x /bgfs/yarbely/mam835/CUT.RUN/Ecoli \
#	-1 /bgfs/yarbely/mam835/CUT.RUN/Pool5/C2_12m_sc_D4/C2_12m_sc_D4.trimmed.R1.fastq \
#	-2 /bgfs/yarbely/mam835/CUT.RUN/Pool5/C2_12m_sc_D4/C2_12m_sc_D4.trimmed.R2.fastq \
#	-S E2_12m_sc_D4_ecoli.sam

#bowtie2 \
#	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 8 \
#	-x /bgfs/yarbely/mam835/CUT.RUN/Ecoli \
#	-1 /bgfs/yarbely/mam835/CUT.RUN/Pool2/PDNC4_CA-HJ-LAP_cC4_Y/PDNC4_CA-HJ-LAP_cC4_Y.trimmed.R2.fastq \
#	-2 /bgfs/yarbely/mam835/CUT.RUN/Pool2/PDNC4_CA-HJ-LAP_cC4_Y/PDNC4_CA-HJ-LAP_cC4_Y.trimmed.R1.fastq \
#	-S PDNC4_CA-HJ-LAP_cC4_Y_ecoli.sam

#bowtie2 \
#	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 8 \
#	-x /bgfs/yarbely/mam835/CUT.RUN/Ecoli \
#	-1 PDNC4_1.trimmed.fq \
#	-2 PDNC4_2.trimmed.fq \
#	-S PDNC4_test_ecoli.sam

samtools stats PDNC4_test.sam > PDNC4_test_stats.txt
samtools stats PDNC4_CA-HJ-LAP_cC4_Y.sam > PDNC4_CA-HJ-LAP_cC4_Y_stats.txt
samtools stats E2_12m_sc_D4.sam > E2_12m_sc_D4_stats.txt

samtools stats PDNC4_test_ecoli.sam > PDNC4_test_ecoli_stats.txt
samtools stats PDNC4_CA-HJ-LAP_cC4_Y_ecoli.sam > PDNC4_CA-HJ-LAP_cC4_Y_ecoli_stats.txt
samtools stats E2_12m_sc_D4_ecoli.sam > E2_12m_sc_D4_ecoli_stats.txt

samtools view -b -h -F 3852 PDNC4_test.sam > PDNC4_test.bam
samtools sort -o PDNC4_test_sorted.bam PDNC4_test.bam
samtools index PDNC4_test_sorted.bam

samtools view -b -h -F 3852 PDNC4_CA-HJ-LAP_cC4_Y.sam > PDNC4_CA-HJ-LAP_cC4_Y.bam
samtools sort -o PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam PDNC4_CA-HJ-LAP_cC4_Y.bam
samtools index PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam

samtools view -b -h -F 3852 E2_12m_sc_D4.sam > E2_12m_sc_D4.bam
samtools sort -o E2_12m_sc_D4_sorted.bam E2_12m_sc_D4.bam
samtools index E2_12m_sc_D4_sorted.bam

#Make a - control for this example by downsampling
#Get a sample (this is an H3K27me3 sample)
module load sratoolkit/3.0.2
fasterq-dump ERR14789999

cutadapt \
	-m 20 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o control.trimmed.R1.fastq -p control.trimmed.R2.fastq \
    ERR14789999_1.fastq ERR14789999_2.fastq
	
	
bowtie2 \
	--end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -p 8 \
	-x /ix1/yarbely/mam835/training/CR_PDNC4/h38.p14 \
	-1 ERR14789999_1.fastq \
	-2 ERR14789999_2.fastq \
	-S control.sam
	
samtools view -b -h -F 3852 control.sam > control.bam
samtools sort -o control.bam control.bam
samtools index control.bam
samtools view -h -b -s 0.01 control.bam > Neg_control.bam
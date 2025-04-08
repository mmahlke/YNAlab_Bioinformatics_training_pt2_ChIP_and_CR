#!/bin/bash
#
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --verbose
#SBATCH --error=/ix1/yarbely/<your_user>/training/CR_PDNC4/job.%J.err
#SBATCH --output=/ix1/yarbely/<your_user>/training/CR_PDNC4/job.%J.out

module load gcc/8.2.0
module load bedtools/2.29.0
module load samtools/1.14

cd /ix1/yarbely/<your_user>/training/PD-NC4

#First, we convert the bam file to a bed file

bedtools bamtobed -bedpe -i PDNC4_test_sorted.bam > PDNC4_test_seacr.bed

#Now we are selecting the information from the bed file that we want to keep, sorting it and sending the output to a new bed file

cut -f 1,2,3,5 PDNC4_test_seacr.bed | sort -k1,1 -k2,2n -k3,3n > PDNC4_test_seacr.clean.bed

#Next we are scaling the clean bed file using our calculated scaling factor and the size of the genome we aligned to, then sending the output to a bedgraph format
## Bedgraph format is the required file format for peak calling with SEACR 
## We need to prepare an index for our alignment file to use in the genomecov command
## We need to unzip the file to build the index

gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna.gz

bedtools genomecov -bg -scale 0.458 -i PDNC4_test_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > PDNC4_test_seacr.bedgraph

#Now, repeat for the other two files
### Remember that the last file will not be scaled, so remove the -scale option

#CA-HJ-LAP_cC4
bedtools bamtobed -i PDNC4_CA-HJ-LAP_cC4_Y_sorted.bam > CA-HJ-LAP_cC4_seacr.bed
cut -f 1,2,3,5 CA-HJ-LAP_cC4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > CA-HJ-LAP_cC4_seacr.clean.bed
bedtools genomecov -bg -scale 0.399 -i CA-HJ-LAP_cC4_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > CA-HJ-LAP_cC4_seacr.bedgraph

#E2_12m_scD4
bedtools bamtobed -i E2_12m_sc_D4_sorted.bam > E2_12m_sc_D4_seacr.bed
cut -f 1,2,3,5 E2_12m_sc_D4_seacr.bed | sort -k1,1 -k2,2n -k3,3n > E2_12m_sc_D4_seacr.clean.bed
bedtools genomecov -bg -i E2_12m_sc_D4_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > E2_12m_sc_D4_seacr.bedgraph

#We also need to prepare our - control 
bedtools bamtobed -i Neg_control.bam > Neg_control_seacr.bed
cut -f 1,2,3,5 Neg_control_seacr.bed | sort -k1,1 -k2,2n -k3,3n > Neg_control_seacr.clean.bed
bedtools genomecov -bg -i Neg_control_seacr.clean.bed -g GCF_000001405.40_GRCh38.p14_genomic.fna.fai > Neg_control_seacr.bedgraph

# All done!
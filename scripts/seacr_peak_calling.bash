#!/bin/bash
#
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --verbose
#SBATCH --error=/ix1/yarbely/<your_user>/training/CR_PDNC4/job.%J.err
#SBATCH --output=/ix1/yarbely/<your_user>/training/CR_PDNC4/job.%J.out

cd /ix1/yarbely/<your_user>/training/PD-NC4

mkdir ./seacr_peaks

#Now call peaks with seacr
## Note about seacr--it requires bedtools module to be loaded
## So if you run peak calling with seacr separately, make sure to load gcc and bedtools as well
module load gcc/8.2.0
module load bedtools/2.29.0
module load seacr/1.3

##PDNC4_test
#With a - control
SEACR_1.3.sh PDNC4_test_seacr.bedgraph Neg_control_seacr.bedgraph norm stringent PDNC4_test_ctrl 
SEACR_1.3.sh PDNC4_test_seacr.bedgraph Neg_control_seacr.bedgraph norm relaxed PDNC4_test_ctrl

#Without a - control
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.00001 non stringent PDNC4_test_0.00001
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.00001 non relaxed PDNC4_test_0.00001
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.01 non stringent PDNC4_test_0.01
SEACR_1.3.sh PDNC4_test_seacr.bedgraph 0.01 non relaxed PDNC4_test_0.01


##CA-HJ-LAP_cC4
#With a - control
SEACR_1.3.sh CA-HJ-LAP_cC4_seacr.bedgraph Neg_control_seacr.bedgraph norm stringent CA-HJ-LAP_cC4 
SEACR_1.3.sh CA-HJ-LAP_cC4_seacr.bedgraph Neg_control_seacr.bedgraph norm relaxed CA-HJ-LAP_cC4

#Without a - control
SEACR_1.3.sh CA-HJ-LAP_cC4_seacr.bedgraph 0.00001 non stringent CA-HJ-LAP_cC4_0.00001
SEACR_1.3.sh CA-HJ-LAP_cC4_seacr.bedgraph 0.00001 non relaxed CA-HJ-LAP_cC4_0.00001
SEACR_1.3.sh CA-HJ-LAP_cC4_seacr.bedgraph 0.01 non stringent CA-HJ-LAP_cC4_0.01
SEACR_1.3.sh CA-HJ-LAP_cC4_seacr.bedgraph 0.01 non relaxed CA-HJ-LAP_cC4_0.01


##E2_12m_sc_D4
#With a - control
SEACR_1.3.sh E2_12m_sc_D4_seacr.bedgraph Neg_control_seacr.bedgraph norm stringent E2_12m_sc_D4 
SEACR_1.3.sh E2_12m_sc_D4_seacr.bedgraph Neg_control_seacr.bedgraph norm relaxed E2_12m_sc_D4

#Without a - control
SEACR_1.3.sh E2_12m_sc_D4_seacr.bedgraph 0.00001 non stringent E2_12m_sc_D4_0.00001
SEACR_1.3.sh E2_12m_sc_D4_seacr.bedgraph 0.00001 non relaxed E2_12m_sc_D4_0.00001
SEACR_1.3.sh E2_12m_sc_D4_seacr.bedgraph 0.01 non stringent E2_12m_sc_D4_0.01
SEACR_1.3.sh E2_12m_sc_D4_seacr.bedgraph 0.01 non relaxed E2_12m_sc_D4_0.01


#Move these new peak files to ./seacr_peaks
mv *stringent.bed ./seacr_peaks/
mv *relaxed.bed ./seacr_peaks/

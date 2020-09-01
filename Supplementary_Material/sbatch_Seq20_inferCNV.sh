#!/bin/bash

# Set job time to 3 days.
#SBATCH --time=3-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name=infercnv_Seq20

# One node.
#SBATCH --nodes=1

# One task
#SBATCH --ntasks=1

# One CPU/core per task, n of threads
#SBATCH --cpus-per-task=8

# 50GB of RAM
#SBATCH --mem=50G

# Who to send mail to.
#SBATCH --mail-user=mjprzy@stanford.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

# which account
#SBATCH --account=ccurtis2
 
# load required modules
source activate snakes

# run script
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R1_T2
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R1_T7
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R1_T12
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R2_T2
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R2_T3
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R2_T7
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing20_C1_T2_T7_T12_ECB_WT/6077_P_C1_ECB_R2_T12

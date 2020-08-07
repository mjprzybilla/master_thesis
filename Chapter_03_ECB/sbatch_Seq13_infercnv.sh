#!/bin/bash

# Set job time to 3 days.
#SBATCH --time=3-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name=infercnv_Seq13

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

# run scripts
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing13_6077_EL_cr3_WT/6077_C1
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing13_6077_EL_cr3_WT/6077_C5
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing13_6077_EL_cr3_WT/6077_D3

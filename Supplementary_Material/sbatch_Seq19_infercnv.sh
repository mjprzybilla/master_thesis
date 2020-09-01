#!/bin/bash

# Set job time to 3 days.
#SBATCH --time=3-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name=infercnv_Seq19

# One node.
#SBATCH --nodes=1

# One task
#SBATCH --ntasks=1

# One CPU/core per task, n of threads
#SBATCH --cpus-per-task=8

# 50GB of RAM
#SBATCH --mem=50GB

# Who to send mail to.
#SBATCH --mail-user=mjprzy@stanford.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

# which account
#SBATCH --account=ccurtis2
 
# load required modules/conda environment
source activate snakes

# run script
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing19_4230_EL_Base_WT/4230_APA6
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing19_4230_EL_Base_WT/4230_PB6
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing19_4230_EL_Base_WT/4230_APA4

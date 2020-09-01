#!/bin/bash

# Set job time to 3 days.
#SBATCH --time=3-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name=infercnv_Seq27

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
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing27_0891EL_WT/0891_PB2
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing27_0891EL_WT/0891_P31
Rscript /home/mjprzy/infercnv_gastric/inferCNV_analysis.R /labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing27_0891EL_WT/0891_P41

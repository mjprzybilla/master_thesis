# Tracking copy number evolution over time

**Seq8_ECB_heatmaps.R** - Script to re-plot the inferCNV output in a format which is overlapping with scDNA and CONICSmat. The representation is shown in Figure 13b.

**inferCNV_analysis.R** - R script with the commands to run inferCNV.

**inferCNV_prepare_Seq8_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the analysis script. Here, data from Sequencing 8, corresponding to D2R2 day 168 is specifically used.

**sbatch_Seq8_inferCNV.sh** - Bash script to run inferCNV on the data from D2R2 day 168 on a slurm cluster.

# Inferring fitness effects from copy number alterations through time

**ECB_01_inferCNV_prepare_Seq8_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the analysis script. Here, data from Sequencing 8, corresponding to D2C2 replicate 2 at day 168 is specifically used.

**ECB_02_inferCNV_analysis.R** - R script with the commands to run inferCNV.

**ECB_03_sbatch_Seq8_inferCNV.sh** - Bash script to run inferCNV on the data from D2R2 replicate 2 at day 168 on a slurm cluster.

**ECB_04_Seq8_ECB_heatmaps.R** - Script to re-plot the inferCNV output with cells ordered according to their cell-subclone-relationship for data from Sequencing 8. The resulting heatmap is shown in Figure 15b.

**ECB_05_Seq8_ECB_heatmap_chr.R** - Script to plot a zoom-in of specific chromosomes from D2C2 replicate 2 day 168, that show copy number alterations between the winning subclone and the base-line subclone. The zoom-in is represented in Figure 15d.

**ECB_06_Seq8_DEG_analysis.R** - Script for the performance of differential gene expression analysis between individual subclones from the expressed cellular barcode system (ECB). Here, Sequencing 8 is used. Further, the proportion of differential expressed genes, and gene set enrichment analysis is performed. The results from this script are shown in Figure 15e-g. 

**ECB_07_inferCNV_prepare_Seq18_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the analysis script. Here, data from Sequencing 18, corresponding to D3C2 replicate 1 at day 168 is specifically used.

**ECB_08_sbatch_Seq18_inferCNV.sh** - Bash script to run inferCNV on the data from D3C2 replicate 1 at day 168 on a slurm cluster.

**ECB_09_Seq18_ECB_heatmaps.R** - Script to re-plot the inferCNV output with cells ordered according to their cell-subclone-relationship for data from Sequencing 18. The resulting heatmap is shown in Figure 16b.

**ECB_10_Seq18_ECB_heatmap_chr.R** - Script to plot a zoom-in of specific chromosomes from D3C2 replicate 1 day 168, that show copy number alterations between the winning subclone and the base-line subclone. The zoom-in is represented in Figure 16d.

**ECB_11_Seq18_DEG_analysis.R** - Script for the performance of differential gene expression analysis between individual subclones from the expressed cellular barcode system (ECB). Here, Sequencing 18 is used. Further, the proportion of differential expressed genes, and gene set enrichment analysis is performed. The results from this script are shown in Figure 16e-g. 

**ECB_12_get_CNV_clones.R** - Script leveraging the hidden markov model (HMM) output from inferCNV to generate copy number clone assignments per cell. For this reason, the *genes x cells* matrix containing predicted copy number states is transformed to a *chromosome arm x cells* matrix, which is subsequently used for the assignment of copy number clones.

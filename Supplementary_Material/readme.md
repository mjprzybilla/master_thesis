# Supplementary Material providing additional information about the thesis

***

**CNV_01_inferCNV_prepare_Seq13_input.R** - R script taking the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the inferCNV analysis script. Here, data from Sequencing 13, corresponding to the *Early* and *Late* time points from donor 2, is used.

**CNV_02_sbatch_Seq13_inferCNV.sh** - Bash script to run inferCNV on the data from organoids of donor 2 at the *Early* and *Late* time point from Sequencing 13 on a slurm cluster. The output from this script is visualized in Supplementary Figure 8.

**CNV_03_inferCNV_prepare_Seq19_input.R** - R script taking the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the inferCNV analysis script. Here, data from Sequencing 19, corresponding to the *Early* and *Late* time points from donor 3, is used.

**CNV_04_sbatch_Seq19_inferCNV.sh** - Bash script to run inferCNV on the data from organoids of donor 3 at the *Early* and *Late* time point from Sequencing 19 on a slurm cluster. The output from this script is visualized in Supplementary Figure 9.

**CNV_05_inferCNV_prepare_Seq27_input.R** - R script taking the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape. Here, data from Sequencing 27, corresponding to the *Early* and *Late* time points from donor 1, is used.

**CNV_06_sbatch_Seq27_inferCNV.sh** - Bash script to run inferCNV on the data from organoids of donor 1 at the *Early* and *Late* time point from Sequencing 27 on a slurm cluster. The output from this script is visualized in Supplementary Figure 8.

***

**ECB_01_inferCNV_prepare_Seq20_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the inferCNV analysis script. Here, data from Sequencing 20, corresponding to D2C1 replicate 1 and 2 at different time points, is used.

**ECB_02_sbatch_Seq20_inferCNV.sh** - Bash script to run inferCNV on the data from D2C1 replicate 1 and 2 at different time points on a slurm cluster.

**ECB_03_Seq20_ECB_heatmaps.R** - Script to re-plot the inferCNV output with cells ordered according to their cell-subclone-relationship for data from Sequencing 20. The resulting heatmap is shown in Supplementary Figure 6b.

**ECB_04_Seq20_ECB_heatmap_chr.R** - Script to plot a zoom-in of specific chromosomes from D2C1 replicate 1 at day 245, that show copy number alterations between the winning subclone and the base-line subclone. The zoom-in is represented in Supplementary Figure 6d.

**ECB_05_Seq20_DEG_analysis.R** - Script for the performance of differential gene expression analysis between individual subclones from the expressed cellular barcode system (ECB). Here, Sequencing 20 is used. Further, the proportion of differential expressed genes, and gene set enrichment analysis is performed. The results from this script are shown in Supplementary Figure 6e-g. 

**ECB_06_inferCNV_prepare_Seq25_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the inferCNV analysis script. Here, data from Sequencing 20, corresponding to D2C2 replicate 1, 2 and 3 at different time points, is used.

**ECB_07_sbatch_Seq25_inferCNV.sh** - Bash script to run inferCNV on the data from D2C2 replicate 1, 2 and 3  at different time points on a slurm cluster.

**ECB_08_Seq25_ECB_heatmaps.R** - Script to re-plot the inferCNV output with cells ordered according to their cell-subclone-relationship for data from Sequencing 25. The resulting heatmaps are shown in Supplementary Figure 4.

***

**LSI_01_Seurat_ZhangEtal.R** - Script to re-analyze the publicly available dataset from Zhang et al., 2019, Cell Reports using Seurat (<https://www.cell.com/cell-reports/pdf/S2211-1247(19)30525-X.pdf>).

**LSI_02_COMBAT.R** - Script for the removal of the organoid-specific gene signature by using reference data from Zhang et al.

**LSI_03_Clustering_UMAP_Zhang.R** - Script for clustering of the scRNA-seq data from Zhang et al., generated with *LSI_01_Seurat_ZhangEtal.R*, which are then subsequently visualized using UMAP. The UMAP representation is then saved for the projection in *LSI_04_Project_Zhang_w_EL_Organoids.R*. This script generates Supplementary Figure 12a and Supplementary Figure 13a.

**LSI_04_Project_Zhang_w_EL_Organoids.R** - Here, the LSI projection itself is implemented. This script projects each individual organoid clone at the *Early* and *Late* time point into the gastric dataset form *LSI_03_Clustering_UMAP_Zhang_v2.R*. Subsequently, the nearest neighbors are quantified and used to quantify the NN cell type frequency. This script generates Supplementary Figure 12b and Supplementary Figure 12c, as well as Supplementary Figure 13b-d.

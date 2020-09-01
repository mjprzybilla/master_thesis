# Evaluating the performance of single-cell RNA copy number estimation methods

**CNV_01_CONICSmat.R** - Script to run CONICSmat on the data from D2C2 replicate 2 at day 168 data. This script is used to generate Supplementary Figure 3.

**CNV_02_honeybadger.R** - Script to run HoneyBADGER on the data from D2C2 replicate 2 at day 168 data.

**CNV_03_inferCNV_prepare_Seq8_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the *CNV_04_inferCNV_analysis.R* script. Here, data from Sequencing 8, corresponding to D2C2 replicate 2 at day 168, is used.

**CNV_04_inferCNV_analysis.R** - R script with the commands to run inferCNV.

**CNV_05_sbatch_Seq8_inferCNV.sh** - Bash script to run inferCNV on the data from D2C2 replicate 2 at day 168 on a slurm cluster.

**CNV_06_Method_Stats.R** - Script to compare the output from CONICSmat, HoneyBADGER and inferCNV to each other and assess the normalization quality. The results of this comparison are visualized in Figure 11c and 11d.

**CNV_07_WGS_CN.R** - Script to visualize the results from shallow WGS (sWGS) for D2C2 replicate 2 at day 168, 245 and 315 in Figure 12b.

**CNV_08_Pseudobulk_Correlation.R** - Script to assess the performance of each method in comparison to the copy number estimations from sWGS. Here, pseudo-bulk modified expression is compared to logR values per gene. The output of this script is depicted in Figure 12c as well as Supplementary Figure 2.

**CNV_09_construct_ploidyMatrix_scDNA.R** - This script generates a *bins x cells* matrix from the output of the cellranger scCNV protocol using bedtools intersect.

**CNV_10_plot_scDNA_heatmap.R** - Script to visualize the *bins x cells* matrix from *CNV_09_construct_ploidyMatrix_scDNA.R*.

**CNV_11_Seq8_ECB_heatmaps.R** - Script to re-plot the inferCNV output in a format which is overlapping with scDNA and CONICSmat. The representation is shown in Figure 13b.

**CNV_12_CONICSmat_heatmap.R** - Script to re-plot the normalized gene expression data from CONICSmat to the common format used in Figure 13c of the Results section.

**CNV_13_CellFrequency.R** - Script to calculate the frequency of altered cell per chromosome arm from CONICSmat, inferCNV and scDNA-seq data. The output of this script is represented in Figure 13d.

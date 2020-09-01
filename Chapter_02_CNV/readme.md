# Comparison of methods for copy number estimations from scRNA-seq

**CNV_01_CONICSmat.R** - Script to run CONICSmat on the data from D2R2 day 168 data. This script is used to generate Supplementary Figure 3.

**CNV_02_CONICSmat_heatmap.R** - Script to re-plot the normalized gene expression data from CONICSmat to the common format used in Figure 13c of the Results section.

**CNV_03_WGS_CN.R** - Script to visualize the results from shallow WGS for D2R2 day 168, 245 and 315 in Figure 12b.

**CNV_03_Method_Stats.R** - Script to compare the output from CONICSmat, HoneyBADGER and inferCNV to each other and assess the normalization quality. The results of this comparison are visualized in Figure 11c and 11d.

**CNV_04_Pseudobulk_Correlation.R** - Script to assess the performance of each method in comparison to the copy number estimations from shallow whole genome sequencing (sWGS). Here, pseudo-bulk modified expression is compared to logR values per gene. The output of this script is depicted in Figure 12c as well as Supplementary Figure 2.

**CNV_04_honeybadger.R** - Script to run HoneyBADGER on the data from D2R2 day 168 data.

**CNV_05_CellFrequency.R** - Script to calculate the frequency of altered cell per chromosome arm from CONICSmat, inferCNV and scDNA-seq data. The output of this script is represented in Figure 13d.

**Seq8_ECB_heatmaps.R** - Script to re-plot the inferCNV output in a format which is overlapping with scDNA and CONICSmat. The representation is shown in Figure 13b.

**inferCNV_analysis.R** - R script with the commands to run inferCNV.

**inferCNV_prepare_Seq8_input.R** - R script that takes the count matrix from Seurat, after the removal of bad quality cells and wrangles it into shape for the analysis script. Here, data from Sequencing 8, corresponding to D2R2 day 168 is specifically used.

**sbatch_Seq8_inferCNV.sh** - Bash script to run inferCNV on the data from D2R2 day 168 on a slurm cluster.

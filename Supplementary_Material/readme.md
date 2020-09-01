# Supplementary Material providing additional information about the thesis.

**LSI_01_Seurat_ZhangEtal.R** - Script to re-analyze the publicly available dataset from Zhang et al., 2019, Cell Reports using Seurat (<https://www.cell.com/cell-reports/pdf/S2211-1247(19)30525-X.pdf>).

**LSI_02_COMBAT.R** - Script for the removal of the organoid-specific gene signature by using reference data from Zhang et al.

**LSI_03_Clustering_UMAP_Zhang_v2.R** - Script for clustering of the scRNA-seq data from Zhang et al., generated with *LSI_01_Seurat_ZhangEtal.R*, which are then subsequently visualized using UMAP. The UMAP representation is then saved for the projection in *LSI_04_Project_Zhang_w_EL_Organoids.R*. This script generates Supplementary Figure 12a and Supplementary Figure 13a.

**LSI_04_Project_Zhang_w_EL_Organoids.R** - Here, the LSI projection itself is implemented. This script projects each individual organoid clone at the *Early* and *Late* time point into the gastric dataset form *LSI_03_Clustering_UMAP_Zhang_v2.R*. Subsequently, the nearest neighbors are quantified and used to quantify the NN cell type frequency. This script generates Supplementary Figure 12b and Supplementary Figure 12c, as well as Supplementary Figure 13b-d.

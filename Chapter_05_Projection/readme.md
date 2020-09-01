# Quantifying evolution towards malignancy in Early and Late organoids

**LSI_01_Seurat_SatheEtal.R** - Script to re-analyze the publicly available dataset from Sathe et al., 2020, Clinical Cancer Research using Seurat (<https://clincancerres.aacrjournals.org/content/clincanres/early/2020/04/03/1078-0432.CCR-19-3231.full.pdf>).

**LSI_02_COMBAT.R** - Script for the removal of the organoid-specific gene signature by using reference data from Zhang et al.

**LSI_03_Clustering_UMAP_Sathe_v2.R** - Script for clustering of the scRNA-seq data from Sathe et al., generated with *LSI_01_Seurat_SatheEtal.R*, which are then subsequently visualized using UMAP. The UMAP representation is then saved for the projection in *LSI_04_Project_Sathe_w_EL_Organoids.R*. This script generates Figure 20a and Supplementary Figure 11a.

**LSI_04_Project_Sathe_w_EL_Organoids.R** - Here, the LSI projection itself is implemented. This script projects each individual organoid clone at the *Early* and *Late* time point into the gastric dataset form *LSI_03_Clustering_UMAP_Sathe_v2.R*. Subsequently, the nearest neighbors are quantified and used to quantify the NN cell type frequency. This script generates Figure 20b,c and Supplementary Figure 11b-d.

The scripts represented here were adapted from Jeffrey Granja's [MPAL publication](https://github.com/GreenleafLab/MPAL-Single-Cell-2019).

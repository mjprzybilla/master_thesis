# LSI projection into gastric tissue to assess early malignant development of gastric organoids

**LSI_01_COMBAT.R** - Script for the removal of the organoid-specific gene signature by using reference data from Zhang et al.

**LSI_02_SatheEtal.R** - Script for clustering of the scRNA-seq data from Sathe et al., which are then subsequently visualized using UMAP. The UMAP representation is then saved for the projection in *LSI_03_Projection.R*. This script generates Figure 20b.

**LSI_03_Projection.R** - Here, the LSI projection is implemented. This script projects each individual organoid clone and the *Early* and *Late* time point into the gastric dataset form *LSI_02*. Subsequently, the nearest neighbors are quantified and used to quantify the NN cell type frequency. This script generates Figure 20c and 20d.

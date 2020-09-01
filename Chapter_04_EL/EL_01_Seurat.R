#############################################################################################################################
##                                                                                                                      
##  Analysis of single-cell RNA-seq data from Early & Late organoids using Seurat
##                                                                                                                      
##  Date: 23 July 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
#############################################################################################################################

## SEURAT
#' Seurat is an analysis tool developed by Rahul Satija's Lab at the New York Genome Center in Lower Mannhattan. This script is
#' based on their featured Seurat introductions which can be found [here](https://satijalab.org/seurat/). There are various vignettes
#' for distinct analytical workflows and procedures provided. The following code is a breakdown of what you can find there. 
#' In case of questions and unclear points in this Markdown document, it might be best to look into their original vignettes for further descriptions. 
#' The `Seurat` version used in this analysis is `Seurat v3.0` upwards. 
#' Also, it might be beneficial to check out their publication [Stuart, Butler et al., 2019, Cell](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8).
#' The first part of this markdown file will focus on single sample analysis. 
#' The fourth and last part will then explore differential gene expression analysis and further visualization of the data. 

#############################################################################################################################

# clear workspace
rm(list = ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "biomaRt", "scater", "DoubletFinder", "viridis", "dittoSeq",
                      "ggsignif", "rstatix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# set up functions which are used 
`%notin%` <- Negate(`%in%`)

# add these lines for the RTAM normalization
# this assumes that you have cloned the following github repo <https://github.com/TiedemannLab/sciCNV>
# path to the cloned github repo
path.code <- "~/sciCNV/sciCNV-Analysis"
source(file.path(path.code, "Mito_umi_gn.R"))
source(file.path(path.code, "RTAM_normalization.R"))
source(file.path(path.code, "sciCNV.R"))
source(file.path(path.code, "Scaling_CNV.R"))
source(file.path(path.code, "CNV_score.R"))
source(file.path(path.code, "sciCNV.R"))

#############################################################################################################################

## READ IN YOUR DATA OF INTEREST
#' Here, we read in all the 10x Genomics matrices, that belong to the set of samples we want to merge. 
#' Merging data is only recommended if you expect lower batch effects than in for instance the analysis of datasets from distinct
#' laboratories and technologies. Here, we assume that we have minimal batch effects in our data from the same set of experiments
#' and technologies. For the dataset from Zhang et al., 2019, Cell we use the raw count matrices from pre-processed Seurat objects. 

## EARLY LATE ORGANOIDS
c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric"

# define output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze"

# define a sample id for the merged object
sample.tmp <- "EL_Seq13_19_27"

# create sample-specific output directory with a "plots" subfolder
dir.create(paste0(o.dir, "/" , sample.tmp))
dir.create(paste0(o.dir, "/" , sample.tmp, "/plots"))

## EARLY LATE ORGANOIDS
sample.list <- list.files(c.dir, pattern = "filtered_feature_bc_matrix", recursive = T, full.names = T, include.dirs = T)
sample.list <- sample.list[grep("Sequencing13|Sequencing19|Sequencing27", sample.list)]

# get sample ids for EL ORGANOIDS
sample.ids <- str_split_fixed(sample.list, "/", 9)[,6]

# get the metadata for the EARLY LATE ORGANOID DATA
seurat.metadata.list <- list.files("/labs/ccurtis2/mjprzy/hashECB_data_updated/", pattern = "Hash_correlation.txt", recursive = T, full.names = T)
seurat.metadata.list <- seurat.metadata.list[grep("Sequencing13|Sequencing19|Sequencing27", seurat.metadata.list)]

# read in all metadata
seurat.metadata <- lapply(seurat.metadata.list, read.table, header = T)

# check each metadata file
for (i in 1:length(seurat.metadata.list)){
  
  # get metadata and change the col names
  metadata.tmp <- seurat.metadata[[i]]
  colnames(metadata.tmp)[1] <- "Cell_Barcode"
  colnames(metadata.tmp)[ncol(metadata.tmp)] <- "HashTag"
  
  # convert to characters
  metadata.tmp$HashTag <- as.character(metadata.tmp$HashTag)
  metadata.tmp$Cell_Barcode <- as.character(metadata.tmp$Cell_Barcode)
  
  # check if complete
  metadata.tmp <- metadata.tmp[complete.cases(metadata.tmp),]
  
  # remove multiplets
  metadata.tmp <- metadata.tmp[-grep("Multiplet", metadata.tmp$HashTag),]
  
  if (i == 1){ # 6077 EL
    
    # remove the samples which we dont want
    metadata.tmp <- metadata.tmp[-grep("BASE", metadata.tmp$HashTag),]
    
    # remove the samples which we dont want
    metadata.tmp <- metadata.tmp[-grep("6077_C1_LATE", metadata.tmp$HashTag),]
    metadata.tmp <- metadata.tmp[-grep("6077_D3_CDKN2A", metadata.tmp$HashTag),]
    
  } else if (i == 3){ # 0891 EL
    # remove the samples which we dont want
    metadata.tmp <- metadata.tmp[-grep("0891_P41_W15", metadata.tmp$HashTag),]
    metadata.tmp <- metadata.tmp[-grep("0891_PB2_W14", metadata.tmp$HashTag),]
    metadata.tmp <- metadata.tmp[-grep("0891_P31_W15", metadata.tmp$HashTag),]
    
    # rename ids
    metadata.tmp[metadata.tmp == "0891_WT_P2"] <- "0891_WT"
    metadata.tmp[metadata.tmp == "0891_P_B2_LATE_102219"] <- "0891_PB2_LATE"
    metadata.tmp[metadata.tmp == "0891_P_31_EARLY"] <- "0891_P31_EARLY"
    metadata.tmp[metadata.tmp == "0891_P_41_EARLY"] <- "0891_P41_EARLY"
    metadata.tmp[metadata.tmp == "0891_P_31_LATE"] <- "0891_P31_LATE"
    metadata.tmp[metadata.tmp == "0891_P_41_LATE"] <- "0891_P41_LATE"
    metadata.tmp[metadata.tmp == "0891_P_B2_EARLY_102219"] <- "0891_PB2_EARLY"
    
  } else { # 4230 EL
    
    # remove the samples which we dont want
    metadata.tmp <- metadata.tmp[-grep("BASE", metadata.tmp$HashTag),]
    
    # rename ids 
    metadata.tmp[metadata.tmp == "LATE_4230_APA6"] <- "4230_APA6_LATE"
    metadata.tmp[metadata.tmp == "EARLY_4230_PB6"] <- "4230_PB6_EARLY"
    metadata.tmp[metadata.tmp == "LATE_4230_APA4"] <- "4230_APA4_LATE"
    metadata.tmp[metadata.tmp == "EARLY_4230_APA6"] <- "4230_APA6_EARLY"
    metadata.tmp[metadata.tmp == "EARLY_4230_APA4"] <- "4230_APA4_EARLY"
    metadata.tmp[metadata.tmp == "LATE_4230_PB6"] <- "4230_PB6_LATE"
    metadata.tmp[metadata.tmp == "4230_WT_Passage_midi"] <- "4230_WT"
    
  }
  
  # remove duplicated barcodes
  metadata.tmp <- metadata.tmp[!duplicated(metadata.tmp$Cell_Barcode),]
  
  # replace the old version
  seurat.metadata[[i]] <- metadata.tmp
}

## EL ORGANOIDS - 10X MATRICES
sample.list <- lapply(sample.list, Read10X)

# list to store seurat.obj.bigects in
seurat.obj.big.list <- list()

############################################################################
##          Read in the data for the EL ORGANOID DATASET
############################################################################ 

# remove cells from the barcode matrices which are not present in the metadata
for(i in seq(length(sample.list))){
  
  # remove the -1 tag
  colnames(sample.list[[i]]) <- str_split_fixed(colnames(sample.list[[i]]), "-1",2)[,1]
  
  # subset matrices to cells which are in metadata
  sample.list[[i]] <- sample.list[[i]][,colnames(sample.list[[i]]) %in% seurat.metadata[[i]]$Cell_Barcode]
  sample.list[[i]] <- sample.list[[i]][,unique(colnames(sample.list[[i]]))]
  
}

# set sample ids for the seurat objects
sample.ids <- c("6077_EL", "4230_EL", "0891_EL")

# iterate over seurat objects and extract raw count matrices
for (i in 1:length(sample.list)){
  
  # extract count matrix
  c.mtx <- sample.list[[i]]
  
  # check dims
  print(dim(c.mtx))
  
  # create new seurat object
  seurat.obj.big <- CreateSeuratObject(counts = c.mtx, 
                                   project = sample.ids[i], 
                                   min.cells = 3, 
                                   min.features = 200)
  
  # add Cell_barcode column
  seurat.obj.big$Cell_Barcode <- rownames(seurat.obj.big@meta.data)
  
  # merge with the Hash metadata
  seurat.obj.big@meta.data <- merge(seurat.obj.big@meta.data, seurat.metadata[[i]], by = "Cell_Barcode")
  
  # store in list
  seurat.obj.big.list[[i]] <- seurat.obj.big
  
}

# combined the seurat objects together
seurat.obj.big <- merge(seurat.obj.big.list[[1]], y = seurat.obj.big.list[c(2:length(seurat.obj.big.list))], add.cell.ids = sample.ids, project = "EL_Seq13_19_27", merge.data = T)

# save identity
seurat.obj.big$sample_ident <- seurat.obj.big$orig.ident

# replace the orig.ident column with HashTag
seurat.obj.big$orig.ident <- seurat.obj.big$HashTag

# check colnames
head(colnames(seurat.obj.big))
tail(colnames(seurat.obj.big))

# check amount of cells
table(seurat.obj.big$orig.ident)

############################################################################
##                     Check out the raw matrix first
############################################################################

c.mtx <- seurat.obj.big@assays$RNA@counts

# calculate counts per cell
counts.per.cell <- Matrix::colSums(c.mtx)

# calculate nubmer of genes per cell
counts.per.gene <- Matrix::rowSums(c.mtx)

# caluclate genes per cell
genes.per.cell <- Matrix::colSums(c.mtx>0) # count gene only if it has non-zero reads mapped.

# calculate cells per genes
cells.per.gene <- Matrix::rowSums(c.mtx>0) # only count cells where the gene is expressed

# plot histograms looking at the calculated metrics
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_preSeurat_QCplots.pdf"), width = 20, height = 20)
hist(log10(counts.per.cell+1),main='counts per cell',col='wheat')
hist(log10(genes.per.cell+1), main='genes per cell', col='wheat')
plot(counts.per.cell, genes.per.cell, log='xy', col='wheat')
title('counts vs genes per cell')
hist(log10(counts.per.gene+1), main='counts per gene', col='wheat')
dev.off()

# plot each cell ranked by their number of genes detected per cell
# represents distribution of library complexity in the sequencing run
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_library_complexity.pdf"), width = 20, height = 20)
plot(sort(genes.per.cell), xlab='cell', log='y', main='genes per cell (ordered)')
dev.off()

#############################################################################################################################

## DATA PRE-PROCESSING
#' Now, we want to create a `Seurat object`, which we can further use to perform quality control and biological meaningful analysis.
#' Upon creation of the `Seurat object`, the number of genes and UMIs (nGene and nUMI) is calculated automatically for every object.
#' For non-UMI data, nUMI represents the sum of the non-normalized values within a cell. In the next steps, we further calculate the
#' percentage of `mitochondrial and ribosomal genes`, as well as the number of `housekeeping genes` per cell. This is adapted from an
#' older Seurat version, where this was done manually. In particular, these HK genes reflect commomn processes active in a cell and
#' hence provide a good global quality measure. The list of HK genes is adapted from [Tirosh et al., 2016, Nature](https://pubmed.ncbi.nlm.nih.gov/27806376/).

message("Start pre-processing..")

# Calculate the percentage of mitochondrial genes
seurat.obj.big[["percent.mt"]] <- PercentageFeatureSet(seurat.obj.big, pattern = "^MT-")

# Calculate percent ribosomal genes
seurat.obj.big[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj.big, pattern = "^RP")

# the quality metrices in the seurat object are stored in the metadata
head(seurat.obj.big@meta.data, 5)

# make a FeatureScatter plot to visualize feature-feature relationships (can be used generally for all columns in metadata)
plot1 <- FeatureScatter(seurat.obj.big, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.obj.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# save plot 1 and plot2 together
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_FeatureScatterplots.pdf"), width = 20, height = 20)
print(plot1 + plot2)
dev.off()

# Load the the list of house keeping genes
hkgenes <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/housekeeping_genes_tirosh.txt")
hkgenes <- as.vector(hkgenes$V1)

# remove hkgenes that were not found
hkgenes.found <- which(toupper(rownames(seurat.obj.big@assays$RNA@counts)) %in% hkgenes)

# calculate the number of housekeeping genes per cell and add it as metadata
n.expressed.hkgenes <- Matrix::colSums(seurat.obj.big@assays$RNA@counts[hkgenes.found, ] > 0)
seurat.obj.big <- AddMetaData(object = seurat.obj.big, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

# plot additional QC
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_QC_plots.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "n.exp.hkgenes"), ncol = 5))
dev.off()

#############################################################################################################################

## REMOVE CELLS NOT MATCHING THE QC THRESHOLDS
#' Once this step is implemented, we should visually examine the quality of the dataset. Following this, we will
#' remove cells based on the calculated QC metrics. However, it is recommended to use very conservative thresholds,
#' as we do not want to remove valuable biological insights. 

# subset the seurat object to the QC-passing cells
seurat.obj.big <- subset(seurat.obj.big, subset = nFeature_RNA > 500 & percent.mt < 20 & n.exp.hkgenes > 55 & percent.ribo < 40)

# plot post QC
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_postQC_plots.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "n.exp.hkgenes"), ncol = 5))
dev.off()

# check amount of cells after QC
table(seurat.obj.big$orig.ident)

#############################################################################################################################

## NORMALIZATION WITH SCTRANSFORM
#' Now that we are having our quality controlled cells in the dataset, we can implement the next step - normalization of the data.
#' The Seurat workflow previously featured a log-normalization, now recommending to use an alternative method, called `sctransform`.
#' In contrast to the previously implemented normalization method, `sctransform` uses a regularied negative binomial regression
#' for the normalization and variance stabilization. It does so by modeling the expression of each gene individually, then grouping
#' them together based on similarity in terms of the determined model. The regression model can incorporate various independent
#' variables, i.e. the sequencing depth, mitochondrial genes or the difference in cell cycle marker expression. For detailed information
#' have a look at the [Github repository](https://github.com/ChristophH/sctransform) or the
#' [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1).
#' 
#' In this step, we will again use the list of cell cycle genes from [Tirosh et al., 2015, Nature](), which is loaded automatically
#' within the Seurat package. We can segregate this list into markers of **G2/M phase** and markers of **S phase**. Based on this we
#' can determine the signals separating non-cycling cells and cycling cells. Using this knowledge, we can provide the cell cycle difference 
#' as an independent variable to `sctransform`, where we will regress out the differences in cell cycle phase amongst proliferating cells 
#' (which are often uninteresting).

# get cell cycle marker genes from Tirosh et al.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# normalize data with SCTransform() before applying the CellCycleScoring
seurat.obj.big <- SCTransform(seurat.obj.big, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))

# perform cell cycle analysis (make sure to specify the "assay" parameter)
seurat.obj.big <- CellCycleScoring(seurat.obj.big, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_RidgePlot.pdf"), width = 20, height = 20)
RidgePlot(seurat.obj.big, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

# calculate the difference between G2M and S pahse scores and regress it out
seurat.obj.big$CC.Difference <- seurat.obj.big$S.Score - seurat.obj.big$G2M.Score

# view cell cycle scores and phase assignments
head(seurat.obj.big[[]])

# normalise again but this time including also the cell cycle scores
seurat.obj.big <- SCTransform(seurat.obj.big, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference'))

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
seurat.obj.big <- RunPCA(seurat.obj.big, features = c(s.genes, g2m.genes))

# save with regressed out cell cycle difference
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_Dimplot.pdf"), width = 20, height = 20)
print(DimPlot(seurat.obj.big, reduction = "pca"))
dev.off()

# save intermediate object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_SCtransform_seurat_obj.rds"))

#############################################################################################################################

## NORMALIZATION WITH RTAM
#' In addition, we will use a previously released normalization method called RTAM to perform normalization. Check out the Supplementary material
#' of the bioRxiv preprint to get further information about the normalization method <https://www.biorxiv.org/content/10.1101/2020.02.10.942607v1>.
#' In brief, this method relies on the assumption that highly covered genes are more reliable for the normalization of all cells in contrast to also
#' incorporate lowly covered genes. 

# run RTAM
message("Run RTAM...")

# run normalization based on raw count RNA assay
RTAM.data <- RTAM_normalization(mat = seurat.obj.big@assays$RNA@counts, method = "RTAM2", Min_nGn  = 500, Optimizing = FALSE)
rownames(RTAM.data) <- rownames(seurat.obj.big@assays$RNA@counts)
colnames(RTAM.data) <- colnames(seurat.obj.big@assays$RNA@counts)

# remove non-expressed genes
RTAM.matrix <- as.matrix(RTAM.data)
colnames(RTAM.matrix) <- colnames(RTAM.data)
final.RTAM.matrix <- as.matrix(RTAM.matrix[which(rowSums(RTAM.matrix != 0) >= 1  ), ])

# add the RTAM normalized count matrix to the seurat object
seurat.obj.big[["RTAM"]] <- CreateAssayObject(data = final.RTAM.matrix)

# switch default assay to RTAM
DefaultAssay(object = seurat.obj.big) <- "RTAM"

# Find variable features for the RTAM assay
seurat.obj.big <- FindVariableFeatures(object = seurat.obj.big, assay = "RTAM", nfeatures = 3000)

# scale data for the RTAM assay
seurat.obj.big <- ScaleData(seurat.obj.big, assay = "RTAM", features = seurat.obj.big@assays$RTAM@var.features)

# switch default back to SCT
DefaultAssay(object = seurat.obj.big) <- "SCT"

message("Finished..")
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_RTAM_seurat_obj.rds"))

#############################################################################################################################

## DIMENSIONALITY REDUCTION
#' With both assays in place we can now run a pca to reduce the dimensionality of the data to our essential components. Each of the analysis
#' steps will be implemented on both assays, so that we can decide which one to use at the very end. Within the next part of this document,
#' we will evaluate how many prinicipal components of the dataset accurately explain most of the variance we observe in our data. To do so, we 
#' will use `Seurat's` visualization functions and manually investigate how many PCs are meaningful for downstream analysis. One of the functions,
#' `DimHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset cells and features are ordered according to there PCA scores.
#' Setting the argument 'cells' shows the 'extreme' cells on both ends of the spectrum. In the end, we will also base our decision on the examination
#' of an ElbowPlot, where PCs are ranked based on the percentage of variance they explain. The number of PCs will then be used for the classification
#' of clusters and nearest neighbors downstream.

seurat.obj.big <- readRDS(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_RTAM_seurat_obj.rds"))

## SCT ASSAY
# run principal component analysis
seurat.obj.big <- RunPCA(seurat.obj.big, npcs = 100, features = seurat.obj.big@assays$SCT@var.features)

# perform Independent component analysis as well
seurat.obj.big <- RunICA(seurat.obj.big, features = seurat.obj.big@assays$SCT@var.features)

# examine and visualize the PCA results in different ways
print(seurat.obj.big[["pca"]], dims = 1:5, nfeatures = 5)

# show PCA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_PCA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj.big, dims = 1:2, reduction = "pca"))
print(DimPlot(seurat.obj.big, reduction = "pca"))
dev.off()

# make dim heatmaps
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_PCA_heatmaps.pdf"), width = 20, height = 20)
print(DimHeatmap(seurat.obj.big, dims = 1, cells = 500, balanced = T))
print(DimHeatmap(seurat.obj.big, dims = 1:20, cells = 500, balanced = T))
dev.off()

# show ICA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_ICA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj.big, dims = 1:3, reduction = "ica"))
dev.off()

# Heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj.big, ndims = 100))
dev.off()

## RTAM
# switch default assay to RTAM
DefaultAssay(object = seurat.obj.big) <- "RTAM"

# run principal component analysis
seurat.obj.big <- RunPCA(seurat.obj.big, npcs = 100, assay = "RTAM", reduction.name = "pca_RTAM", reduction.key = "pca_RTAM_", features = seurat.obj.big@assays$RTAM@var.features)

# perform Independent component analysis as well
seurat.obj.big <- RunICA(seurat.obj.big, assay = "RTAM", reduction.name = "ica_RTAM", reduction.key = "ica_RTAM_", features = seurat.obj.big@assays$RTAM@var.features)

# examine and visualize the PCA results in different ways
print(seurat.obj.big[["pca_RTAM"]], dims = 1:5, nfeatures = 5)

# show PCA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_PCA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj.big, dims = 1:2, reduction = "pca_RTAM"))
print(DimPlot(seurat.obj.big, reduction = "pca_RTAM"))
dev.off()

# make dim heatmaps
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_PCA_heatmaps.pdf"), width = 20, height = 20)
print(DimHeatmap(seurat.obj.big, dims = 1, cells = 500, balanced = T, reduction = "pca_RTAM"))
print(DimHeatmap(seurat.obj.big, dims = 1:20, cells = 500, balanced = T, reduction = "pca_RTAM"))
dev.off()

# show ICA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_ICA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj.big, dims = 1:3, reduction = "ica_RTAM"))
dev.off()

# Heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj.big, ndims = 100, reduction = "pca_RTAM"))
dev.off()

# switch default back to SCT
DefaultAssay(object = seurat.obj.big) <- "SCT"

#############################################################################################################################

## CELL CLUSTERING AND EMBEDDINGS
#' With the dimensioality reduction in place, we can apply graph-based clustering method where we construct a KNN graph based
#' on the euclidean distance in PCA space. We further use modularity optimization techniques such as the louvain algorithm. In brief,
#' the higher the resolution we set in this clustering, the greater the number of clusters generally is. A resolution of 0.4 - 1.2
#' typically give good results for 3k cells. Generally, the goal of cell clustering and embedding is to learn the underlying manifold
#' of the data in order to place similar cells together in a low-dimensional space. 

## SCT
# find NNs
seurat.obj.big<- FindNeighbors(seurat.obj.big, dims = 1:50, assay = "SCT", reduction = "pca")

# find clusters
seurat.obj.big <- FindClusters(seurat.obj.big, resolution = 0.9)

# run umap
seurat.obj.big <- RunUMAP(seurat.obj.big, dims = 1:50, assay = "SCT", reduction = "pca")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_UMAPplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", label = T))
dev.off()

# run tsne
seurat.obj.big <- RunTSNE(seurat.obj.big, reduction = "pca", dims = 1:50)

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_TSNEplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "tsne", label = T))
dev.off()

## RTAM
# switch default assay to RTAM
DefaultAssay(object = seurat.obj.big) <- "RTAM"

# find NNs
seurat.obj.big<- FindNeighbors(seurat.obj.big, dims = 1:50, assay = "RTAM", reduction = "pca_RTAM")

# find clusters
seurat.obj.big <- FindClusters(seurat.obj.big, resolution = 0.9)

# run umap
seurat.obj.big <- RunUMAP(seurat.obj.big, dims = 1:50, assay = "RTAM", reduction = "pca_RTAM", reduction.name = "umap_RTAM", reduction.key = "umapRTAM_")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_UMAPplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap_RTAM", label = T))
dev.off()

# run tsne
seurat.obj.big <- RunTSNE(seurat.obj.big, dims = 1:50, assay = "RTAM", reduction = "pca_RTAM", reduction.name = "tsne_RTAM", reduction.key = "tsneRTAM_")

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_TSNEplot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "tsne_RTAM", label = T))
dev.off()

# switch default back to SCT
DefaultAssay(object = seurat.obj.big) <- "SCT"

# save intermediate object
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_DimRed_seurat_obj.rds"))

#############################################################################################################################

## DOUBLET DETECTION WITH DOUBLETFINDER
#' Uses a fully-processed seurat object (after tsne has been run) needs PCs (range 1:10 for example).
#' pN - number of artifical generated doublets (default 25%)
#' pK - definition PC neighborhodd size
#' nExp - defines threshold for doublet/singlet predictions
#' pK Identification with sctransform used (no ground-truth - can be run with ground-truth as well)
#' Further information: https://github.com/chris-mcginnis-ucsf/DoubletFinder

# follow the DoubletFinder vignette for detecting and removing doublets
sweep.res.list.seurat.obj.big <- paramSweep_v3(seurat.obj.big, PCs = 1:75, sct = T)
sweep.stats_seurat.obj.big <- summarizeSweep(sweep.res.list.seurat.obj.big, GT = FALSE)

# plot Mean-variance normalized bimodality coefficient (bcmvn) 
# ground-truth-agnostic metric that coincides with pK
bcmvn.seurat.obj.big <- find.pK(sweep.stats_seurat.obj.big)

# plot pk value
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_pKvalue_plot.pdf"))
pK=as.numeric(as.character(bcmvn.seurat.obj.big$pK))
BCmetric=bcmvn.seurat.obj.big$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
dev.off()

# get annotatios for clustering
annotations <- seurat.obj.big@meta.data$seurat_clusters

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.075*length(seurat.obj.big@meta.data$barcode))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seurat.obj.big <- doubletFinder_v3(seurat.obj.big, PCs = 1:75, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
seurat.obj.big <- doubletFinder_v3(seurat.obj.big, PCs = 1:75, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK_choose, "_", nExp_poi), sct = TRUE)

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Doublet_plot.pdf"))
print(DimPlot(seurat.obj.big, reduction = "umap", group.by = paste0("DF.classifications_0.25_",pK_choose, "_", nExp_poi.adj)))
dev.off()

# assign a new column name to the classification column
colnames(seurat.obj.big@meta.data)[ncol(seurat.obj.big@meta.data)] <- "Doublet_Score"

# remove doublets from seurat object
seurat.obj.big <- subset(seurat.obj.big, subset = Doublet_Score != "Doublet")

# save interemediate objcet
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_seurat_obj.rds"))

#############################################################################################################################

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
#' In this last part, we will create some plots and investigate distinct clusters for the expression of different genes, most
#' of which have been reported in the literature. We will use this manual investigation, to manually, but only roughly, 
#' assign possible identities to clusters. 

# # distinct Hash samples
# sample.tmp <- "EL_Seq13_19_27"

# # read seurat objects
seurat.obj.big <- readRDS(file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_Adapted_seurat_obj.rds"))
Idents(seurat.obj.big) <- "SCT_snn_res.1"

## CELL MARKERS FROM ZHANG ET AL.
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_featurePlot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "CD79A", "CD19", "OLFM4", "LGR5", "SOX2","CCKBR", "FABP1", "FABP2", "CA1", "VIL1", "MUC6", "TFF2", "CHGA", "TAC1", "TPH1", "CHGB"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CDK1", "MKI67", "CEACAM5", "CEACAM6", "PGC", "CXCL3", "IL8", "COL1A2", "LUM", "DCN", "PDPN", "FAP", "COL3A1", "COL6A1", "VWF", "ENG", "MCAM"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "CD14", "CD68", "CSF1R", "MYL2", "ACTA2", "CD14", "CD68", "CSF1R", "MYL2"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("ACTA2", "TPSAB1", "TPSB2", "PGA3", "PGA4", "LIPF", "CD2", "CD3D", "CD3E", "CD3G", "ATP4A", "ATP4B", "GAST", "GHRL", "SST"))
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("OLFM4", "PHLDA1", "LEFTY1")) # stem cells
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM6", "BAX", "CCND2")) # cancer cells
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM5", "FABP1", "CDH17")) # non-specific cancer cells - also in enterocytes
dev.off()

## EPITHELIAL AND NON-EPITHELIAL
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_Umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj.big, features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4"))
FeaturePlot(seurat.obj.big, features = c("CD4", "VIM", "ACTA2", "PTPRC"))
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_expressionPlot.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4")))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CD4", "VIM", "ACTA2", "PTPRC")))
dev.off()

# plot probability distributions across clusters for the top5 genes indicated in the Supplementary data from Zhang et al.
# However, markers shared by different cell types were removed from both.
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_exp_violin_plot.pdf"), width = 20, height = 20)
# check for PMCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GKN1", "GKN2", "MUC5AC", "TFF1", "DPCR1")) + labs(title = "PMCs"))
# check for MSCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("OLFM4", "RPL7", "CLDN4","TSPAN8", "REG1A")) + labs(title = "MSCs"))
# check for enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("FABP1", "FABP2", "RBP2", "ANPEP", "APOA4")) + labs(title = "Enterocytes"))
# check for GMC
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PRR4", "C6orf58","MUC6", "TFF2", "LTF")) + labs(title = "GMCs"))
# check for enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CHGA", "PCSK1N", "SCG5", "CHGB", "TPH1")) + labs(title = "Enteroendocrine")) 
# check for PCs
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("TOP2A", "MKI67", "UBE2C", "HMGB2", "PTTG1")) + labs(title = "PCs"))
# check for cancer cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("REG4","CLDN7","KRT18", "LGALS3", "CEACAM6")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CLDN3", "CST1", "MUC3A", "CLDN4", "PI3", "UBD")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CDH17", "PRAP1", "UBE2C", "CCL20", "LCN2", "SERPINB5")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("RRM2", "MYBL2", "MMP7", "TPX2", "MISP", "TMPRSS4")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("RMRP", "CLDN1", "GPRC5A", "CLRN3", "CXCL1", "MSLN")) + labs(title = "Cancer Cells"))
# check for neck like cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("MIA", "CXCL3", "CXCL2", "CXCL17", "CLU")) + labs(title = "Neck-like cells"))
# check for Goblet cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "ZG16")) + labs(title = "Goblet Cells"))
# check for Chief cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PGA3", "PGA4", "LIPF", "CHIA", "PDIA2")) + labs(title = "Chief Cells")) 
# check for parietal cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ATP4A", "ATP4B", "VEGFB")) + labs(title = "Parietal cells"))
# check for endocrine cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GAST", "GHRL", "SST")) + labs(title = "Endocrine cells"))
##GAO DATASET
# check for Chief Stem Cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("GIF", "LRIG1", "PROCR")) + labs(title = "Chief Stem cells"))
# check for Enteroendocrine
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("SST", "PYY", "PROX1", "TPH1", "REG4", "NEUROD1", "GIP")) + labs(title = "Enteroendocrine cells"))
# check for Enterocytes
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("CA2", "KRT20", "ALPI", "TREH", "LCT", "MME", "CDH1", "VIL1", "AQP8", "SI", "CDX2")) + labs(title = "Enterocytes"))
# check for paneth cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("NOTUM", "NUPR1", "PLA2G2A")) + labs(title = "Paneth Cells"))
# check for tuft cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("DCLKK1", "TRPM5", "PTGS1", "RGS13")) + labs(title = "Tuft Cells"))
# check for Stem cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PROM1", "MEX3A", "TERT", "LGR5", "SOX9", "CD24", "ALCAM", "PROCR")) + labs(title = "Stem Cells"))
# check for PROCRhigh progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("PPIH", "SGK1", "GPT2")) + labs(title = "PROCR High Progenitor Cells"))
# check for Neck progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("TFF2", "ECE1", "PROM1")) + labs(title = "Neck-Progenitor Cells"))
# check for HES1high progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("LINC01578", "PLXDC1")) + labs(title = "HES1 High Progenitor Cells"))
# check for parietal progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ATP4A", "ATP4B", "CFH", "FGD5")) + labs(title = "Parietal Progenitor Cells"))
# check for Pit progenitor cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("LRIG1", "AXIN2", "CD44", "ACTC1")) + labs(title = "Pit Progenitor Cells"))
# check for Pit Cells
print(VlnPlot(object = seurat.obj.big, assay = "SCT", features = c("ID1", "MBD1", "GREM", "RSPO2")) + labs(title = "Pit Cells"))
dev.off()

###########################################################################
#               ANALYSE THE GENE EXPRESSION PER CLUSTER
###########################################################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
seurat.obj.big.markers <- seurat.obj.big.markers[order(seurat.obj.big.markers$cluster, seurat.obj.big.markers$avg_logFC, decreasing = c(F,T)),]
top20_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(top20_seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_top20_cluster_markers_pct.filtered.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_cluster_markers_pct.filtered.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# seurat.obj.big.markers <- read.table(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_cluster_markers.txt"), col.names = T)

## TOP MARKERS
# top 5 markers (or all markers if less than 10) for each cluster.
top3_seurat.obj.big.markers <- seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_ExpressionHeatmap", ".pdf", sep = ""), width = 13, height = 15)
print(DoHeatmap(seurat.obj.big, features = top3_seurat.obj.big.markers$gene, size = 5, draw.lines = T, hjust = 0.4, angle = 45, disp.min = -3, disp.max = 3) + 
        scale_fill_gradientn(colors = c("blue", "lightgrey", "red")) + 
        NoLegend() +  theme(text = element_text(size = 16, color = "black")))
dev.off()

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_dataSlot_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, assay = "SCT", slot = "data", features = top5_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue")))
dev.off()

###########################################################################
#             ANALYSE THE GENE EXPRESSION PER CLONE
###########################################################################

# set idents to clone id
Idents(seurat.obj.big) <- "orig.ident" # EL combined

# check contribution of each clone to the overall cell number
table(seurat.obj.big$orig.ident)

# EL
levels(seurat.obj.big) <- c("6077_WT", "6077_C1_EARLY", "6077_C5_EARLY", "6077_C5_LATE", "6077_D3_EARLY", "6077_D3_LATE", "4230_WT", "4230_PB6_EARLY", "4230_PB6_LATE",
                            "4230_APA4_EARLY", "4230_APA4_LATE", "4230_APA6_EARLY", "4230_APA6_LATE", "0891_WT", "0891_PB2_EARLY", "0891_PB2_LATE", "0891_P31_EARLY",
                            "0891_P31_LATE", "0891_P41_EARLY", "0891_P41_LATE")

levels(seurat.obj.big) <- c("6077_WT", "4230_WT", "0891_WT", 
                            "6077_C1_EARLY", "6077_C5_EARLY", "6077_D3_EARLY", "4230_PB6_EARLY", "4230_APA4_EARLY", "4230_APA6_EARLY", "0891_PB2_EARLY", "0891_P31_EARLY", "0891_P41_EARLY",
                            "6077_C5_LATE", "6077_D3_LATE", "4230_PB6_LATE", "4230_APA4_LATE", "4230_APA6_LATE", "0891_PB2_LATE", "0891_P31_LATE", "0891_P41_LATE")

clone.seurat.obj.big.markers <- FindAllMarkers(seurat.obj.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_clone_seurat.obj.big.markers <- clone.seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top5_clone_seurat.obj.big.markers <- clone.seurat.obj.big.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(clone.seurat.obj.big.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_clone_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
clone.seurat.obj.big.markers <- read.table(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_clone_markers.txt"), header = T)

## TOP MARKERS
# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_top5_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, features = top5_clone_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")) + NoLegend())
dev.off()

# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_top5_SCT_dataSlot_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj.big, assay = "SCT", slot = "data", features = top5_clone_seurat.obj.big.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue"))+ NoLegend())
dev.off()

# add clone id metadata
seurat.obj.big$clone_id <- paste0(str_split_fixed(seurat.obj.big$HashTag, "_", 3)[,1], "_",str_split_fixed(seurat.obj.big$HashTag, "_", 3)[,2])

# remove these two clones from the heatmap
plotting.obj <- subset(seurat.obj.big, clone_id %notin% c("6077_C1", "0891_PB2"))

levels(plotting.obj) <- c("0891_WT", "0891_P31_EARLY", "0891_P31_LATE", "0891_P41_EARLY", "0891_P41_LATE",
                          "6077_WT", "6077_C5_EARLY", "6077_C5_LATE", "6077_D3_EARLY", "6077_D3_LATE",
                          "4230_WT", "4230_PB6_EARLY", "4230_PB6_LATE", "4230_APA4_EARLY",  "4230_APA4_LATE", "4230_APA6_EARLY", "4230_APA6_LATE")

clone.plotting.obj <- FindAllMarkers(plotting.obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
top5_clone_seurat.obj.big.markers <- clone.plotting.obj %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


## TOP MARKERS
# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_top5_subset_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10, height = 12)
print(DoHeatmap(plotting.obj, features = top5_clone_seurat.obj.big.markers$gene, size = 5, draw.lines = T, angle = 45, disp.min = 0, disp.max = 3) + 
        scale_fill_gradientn(colors = c("white", "red")) + 
        theme(text = element_text(size = 16, color = "black")))
dev.off()

###########################################################################
#                       CELL TYPE MARKER HEATMAP
###########################################################################
# remove these two clones from the heatmap
plotting.obj <- subset(seurat.obj.big, clone_id %notin% c("6077_C1", "0891_PB2"))

genes <- c("FABP1", "VIL1", "TFF3", "MUC5AC", "MUC6", "OLFM4", "SOX2",
           "MKI67", "BIRC5", "CDK1", "CD44", "TCF4", "ID1", "ID2", "ID3", "BMP2", "REG4", "KRT18", "CD9", "CLDN7",
           "CEACAM6", "S100A14", "LGALS3", "TSPAN8", "S100A16", "SERPINB5", "GDA", "DMBT1", "PCK1", "KRT20", "HSD17B2", "CDH17", "LGALS4", "MUC13",
           "S100P", "FCGBP", "CEACAM5", "GDA", "LYZ", "KRT20", "ADH1C", "AKR1B10",
           "CDCA7", "SLC6A14", "AADAC", "HSD17B2", "GCNT3", "CDX1", "CDX2", "CLDN3", "WFDC2")

mtor <- c("AK4", "ALDOA", "ASNS", "ATP5MC1", "AURKA", "BCAT1", "BHLHE40", "BTG2", "CALR", "CD9", "CDKN1A", "CYB5B", "DAPP1", "DDIT3", "DDIT4", "DDX39A", 
          "DHCR24", "DHCR7", "DHFR", "EBP", "EGLN3", "ENO1", "ERO1A", "FADS1", "G6PD", "GAPDH", "GBE1", "GCLC", "GLRX", "GOT1", "GPI", "HK2", "HMGCR", "HMGCS1",  "HSP90B1",
          "HSPA5", "HSPA9", "IDH1", "IDI1", "IFRD1", "INSIG1", "LDHA", "LDLR", "MTHFD2", "NUPR1", "P4HA1", "PGK1", "PHGDH", "PLK1", 
          "PLOD2", "PPP1R15A", "PRDX1", "PSAT1", "SCD", "SDF2L1", "SKAP2", "SLC1A5", "SLC2A1", "SLC7A11", "SLC7A5", "SORD", "SQLE", "SQSTM1", "STARD4", "SYTL2",
          "TMEM97", "TPI1", "TRIB3", "TUBA4A", "TXNRD1", "XBP1")

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_GeneList_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(plotting.obj, features = genes, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend())
dev.off()

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_mTOR_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(plotting.obj, features = mtor, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")) + NoLegend())
dev.off()

genes <- c("FABP1", "VIL1", "TFF3", "MUC5AC", "MUC6", "OLFM4", "SOX2",
           "MKI67", "BIRC5", "CDK1", "CD44", "TCF4", "ID1", "ID2", "ID3", "BMP2", "REG4", "KRT18", "CD9", "CLDN7",
           "CEACAM6", "S100A14", "LGALS3", "TSPAN8", "S100A16", "SERPINB5", "GDA", "DMBT1", "PCK1", "KRT20", "HSD17B2", "CDH17", "LGALS4", "MUC13",
           "S100P", "FCGBP", "CEACAM5", "GDA", "LYZ", "KRT20", "ADH1C", "AKR1B10",
           "CDCA7", "SLC6A14", "AADAC", "HSD17B2", "GCNT3", "CDX1", "CDX2", "CLDN3", "WFDC2", "AK4", "ALDOA", "ASNS", "ATP5MC1", "AURKA", "BCAT1", "BHLHE40", "BTG2", "CALR", "CD9", "CDKN1A", "CYB5B", "DAPP1", "DDIT3", "DDIT4", "DDX39A", 
           "DHCR24", "DHCR7", "DHFR", "EBP", "EGLN3", "ENO1", "ERO1A", "FADS1", "G6PD", "GAPDH", "GBE1", "GCLC", "GLRX", "GOT1", "GPI", "HK2", "HMGCR", "HMGCS1",  "HSP90B1",
           "HSPA5", "HSPA9", "IDH1", "IDI1", "IFRD1", "INSIG1", "LDHA", "LDLR", "MTHFD2", "NUPR1", "P4HA1", "PGK1", "PHGDH", "PLK1", 
           "PLOD2", "PPP1R15A", "PRDX1", "PSAT1", "SCD", "SDF2L1", "SKAP2", "SLC1A5", "SLC2A1", "SLC7A11", "SLC7A5", "SORD", "SQLE", "SQSTM1", "STARD4", "SYTL2",
           "TMEM97", "TPI1", "TRIB3", "TUBA4A", "TXNRD1", "XBP1")

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_TotalGeneList_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 16.5)
print(DoHeatmap(plotting.obj, features = genes, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend())
dev.off()


###########################################################################
#            GENERATE PSEUDOBULKS PER CLONE FOR EACH GENE
###########################################################################

## EL ORGANOIDS
# add timepoint metadata
seurat.obj.big$timepoint <- str_split_fixed(seurat.obj.big$HashTag, "_", 3)[,3]
seurat.obj.big@meta.data[grep("WT", seurat.obj.big$HashTag), "timepoint"] <- "WT"

# add clone id metadata
seurat.obj.big$clone_id <- paste0(str_split_fixed(seurat.obj.big$HashTag, "_", 3)[,1], "_",str_split_fixed(seurat.obj.big$HashTag, "_", 3)[,2])

# calculate average gene expression per clone id for Aziz
clone.averages <- AverageExpression(seurat.obj.big)
SCT.pseudobulk <- clone.averages$SCT

write.table(SCT.pseudobulk, "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/EL_Seq13_19_27/EL_Seq13_19_27_pseudobulk.txt", col.names = T, row.names = T, sep = "\t", quote = F)

############################################################################
##        MAKE PLOTS FOR THE UMAP EMBEDDING FOR EL ORGANOIDS
############################################################################
# save final cell type annotated seurat_obj
seurat.obj.big <- readRDS(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_CellType_seurat_obj.rds"))

# # set it to seurat clusters again
# Idents(seurat.obj.big) <- "seurat_clusters"

# get colours for the plot
cols <- c(viridis(60), magma(60), inferno(60))

## CLONE ID
plot1 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "clone_id") +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("CLONE IDs", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_CloneID.pdf"), width = 6, height = 5)
print(plot1 + NoLegend()) 
dev.off()

plot2 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "timepoint", pt.size = 1.5) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=40, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=40, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  scale_color_manual(labels = c("Early", "Late", "WT"), values = c("#EA9F37", "#782867", "#4F9E4C")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 8)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Timepoint.pdf"), width = 6, height = 5)
print(plot2 + NoLegend()) 
dev.off()

plot3 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "SCT_snn_res.0.6", pt.size = 1.5, label = T, label.size = 6) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=20, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Clusters", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Clusters.pdf"), width = 6, height = 5)
print(plot3 + NoLegend()) 
dev.off()

plotlist <- list(plot1, plot2, plot3)

# plot the list
nCol <- 3
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Comparison_plot.pdf"), width = 60, height = 20)
cowplot::plot_grid(plotlist = plotlist, ncol = nCol)
dev.off()

## PATIENT ID
plot4 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "old.ident", pt.size = 1.5) +
  xlab("UMAP-1") + 
  ylab("UMAP-2") +
  theme_classic() +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=20, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("CLONE IDs", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_PatientID.pdf"), width = 12, height = 10)
print(plot4) 
dev.off()

# show seurat color palette
# show_col(hue_pal()(12))

############################################################################
##                PLOT INDIVIDUAL MARKER GENES HERE
############################################################################

# # plot feature expression on a tSNE or PCA plot
# pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_MucosalMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
# FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "MUC1"), pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#B79F00"))
# dev.off()

## PROLIFERATIVE MARKERS
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ProliferativeMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("TOP2A", "MKI67"), pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#B79F00"))
dev.off()

## ENTEROENDOCRINE MARKERS
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_EnteroendocrineMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("PCSK1N"), pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#C77CFF")) + scale_color_gradientn(colours = c("lightgrey", "#00BA38"),  limits = c(0, 3))
dev.off()

## MUCOSAL STEM CELLS
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_MucosalStemMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("OLFM4"), pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#F8766D"))
dev.off()

## ENTEROCYTES
# plot feature expression on a tSNE or PCA plot
p1 <- FeaturePlot(seurat.obj.big, reduction = "umap", features = c("FABP1", "KRT20"), combine = F, pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#F564E3"))
fix.sc <- scale_color_gradientn(colours = c("lightgrey", "#F564E3"),  limits = c(0, 6))
p2 <- lapply(p1, function (x) x + fix.sc)
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_EnterocytesMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
CombinePlots(p2)
dev.off()

## NECK-LIKE CELLS
p1 <- FeaturePlot(seurat.obj.big, reduction = "umap", features = c("PGC", "REG1A"), pt.size = 1, combine = F, ncol = 2, order = T, cols = c("lightgrey", "#C77CFF"), min.cutoff = 0, max.cutoff = 6)
fix.sc <- scale_color_gradientn(colours = c("lightgrey", "#C77CFF"),  limits = c(0, 6))
p2 <- lapply(p1, function (x) x + fix.sc)
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_NecklikeMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
CombinePlots(p2)
dev.off()

## GOBLET CELLS
# plot feature expression on a tSNE or PCA plot
p1 <- FeaturePlot(seurat.obj.big, reduction = "umap", features = c("WFDC2", "MUC5B"), combine = F, pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#00B4F0"), min.cutoff = 0)
fix.sc <- scale_color_gradientn(colours = c("lightgrey", "#00B4F0"),  limits = c(0, 6))
p2 <- lapply(p1, function (x) x + fix.sc)
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GobletMarkers_FeaturePlot.pdf"), width = 11, height = 4.5)
CombinePlots(p2)
dev.off()

############################################################################
##           PLOT INDIVIDUAL MARKER GENES (EXTENDED LIST) HERE
############################################################################

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangMucosalTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "TFF2"), pt.size = 1, ncol = 2, order = T)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("PSCA", "S100P", "MUC1", "CA2"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangMucosalStemTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("OLFM4"), pt.size = 1, ncol = 2, order = T)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("TSPAN8", "HSPD1", "NPM1", "EDN1"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangEnterocytesTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("FABP1", "RBP2", "FABP2", "ANPEP"), pt.size = 1, ncol = 2, order = T)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("PHGR1", "KRT20", "PRAP1", "ALDOB"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangGlandMucosalTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("PGC", "MUC6", "CLU", "REG1A"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangEnteroendocrineTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MDK", "PCSK1N", "SCG5", "CPE"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangNeckLikeTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CXCL3", "LYZ", "PGC", "CXCL2"), pt.size = 1, ncol = 2, order = T)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MIA", "MMP1", "CXCL17"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangCancerTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("REG4", "CLDN7", "TSPAN8", "KRT18"), pt.size = 1, ncol = 2, order = T)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("LGALS3", "CEACAM6", "CD9", "TM4SF20"), pt.size = 1, ncol = 2, order = T)
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_ZhangGobletTop10_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("SPINK4", "TFF3", "ITLN1", "ZG16B"), pt.size = 1, ncol = 2, order = T)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("FCGBP", "WFDC2", "KLK1"), pt.size = 1, ncol = 2, order = T)
dev.off()

############################################################################
##         PLOT INDIVIDUAL MARKER GENES (SOME GENES OF INTERES) HERE
############################################################################

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest1_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("MUC5AC", "TFF1", "CEACAM6", "PGC"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest2_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("WFDC2", "MUC5B", "REG4", "OLFM4"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest3_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("TSPAN8", "VIL1", "AKR1C1", "LYZ"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest4_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("DMBT1", "FABP1", "GSTA1", "PSCA"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_MalignantMarker_FeaturePlot.pdf"), width = 11, height = 4.5)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c("CEACAM6", "CEACAM5", "KLK10"), pt.size = 1, ncol = 2, order = T, cols = c("lightgrey", "#00C08B"))
dev.off()

## PLOT DIFFERENTIALLY EXPRESSED GENES IN THESE CLUSTERS
genes <- unique(top3_seurat.obj.big.markers$gene)
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Top3Genes_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[1:9]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[10:18]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[19:27]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[28:36]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[37:45]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[46:54]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[55:63]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[64:72]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[73:81]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[82:90]), pt.size = 0.25, ncol = 3)
FeaturePlot(seurat.obj.big, reduction = "umap", features = c(genes[81:96]), pt.size = 0.25, ncol = 3)
dev.off()

############################################################################
##                  PLOT SOME MARKERS IN A DIFFERENT WAY
############################################################################
# set the factors here
seurat.obj.big$old.ident <- factor(seurat.obj.big$old.ident, levels = c("0891_EL", "6077_EL", "4230_EL"))
seurat.obj.big$SCT_snn_res.1 <- factor(seurat.obj.big$SCT_snn_res.1, levels = c(0:23))

## PLOT SUBFIGURE 18E
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_Cluster_CompositionPlot.pdf"), width = 7.5, height = 6.5)
dittoBarPlot(seurat.obj.big, "old.ident", group.by = "SCT_snn_res.1",
             main = NULL, 
             xlab = "Clusters", # NULL = remove
             ylab = "Fraction of Cells",
             legend.title = "Clusters",
             var.labels.rename = c("P1", "P3", "P2"),
             x.reorder = c(1,2, 13, 18:24, 3:12, 14:17),
             x.labels.rotate = FALSE)
dev.off()

## PLOT SUBFIGURE 19C
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_MucosalRidgePlot_PatientID.pdf"), width = 8, height = 4)
RidgePlot(seurat.obj.big, features = c("MUC5AC", "MUC1"), group.by = "old.ident")
dev.off()

# get the normalized matrix
matrix <- seurat.obj.big@assays$SCT@data

# make a dataframe with the TFF genes
TFF.df <- as.data.frame(t(matrix[c("TFF1", "TFF2", "TFF3"),]))
TFF.df$patient_id <- str_split_fixed(rownames(TFF.df), "_", 3)[,1]
TFF.melt <- melt(TFF.df)

# rename the ids
TFF.melt[TFF.melt$patient_id == "0891", "patient_id"] <- "P1"
TFF.melt[TFF.melt$patient_id == "6077", "patient_id"] <- "P2"
TFF.melt[TFF.melt$patient_id == "4230", "patient_id"] <- "P3"

# order
TFF.melt$patient_id <- factor(TFF.melt$patient_id, levels = c("P1", "P2", "P3"))

# Visualize the expression profile
my_comparisons <- list(c("P1", "P2"), c("P1", "P3"), c("P2", "P3"))
ggviolin(TFF.melt, x = "patient_id", y = "value", fill = "patient_id", 
         add = "jitter", legend = "none", facet.by = "variable", short.panel.labs = TRUE, add.params = list(size = 0.1, jitter = 0.2)) +
  rotate_x_text(angle = 45) +       # Add global annova p-value
  stat_compare_means(label = "p.signif", p.adjust.method = "bonferroni", method='t.test', 
                     comparisons = my_comparisons) +
  theme_classic() + 
  labs(x = "Patients", y = "Normalized Expression") +
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 20, face = "bold" ),
        plot.title = element_text(colour = "black", size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "None")
ggsave(paste0(o.dir, "/", sample.tmp, "/TFF_gene_VlnPlot.pdf"), width = 10, height = 4.5)

# make a dataframe with the "cancer" genes
malignant.df <- as.data.frame(t(matrix[c("CEACAM6", "REG4", "CD9"),]))
malignant.df$patient_id <- str_split_fixed(rownames(malignant.df), "_", 3)[,1]
malignant.melt <- melt(malignant.df)

# rename the ids
malignant.melt[malignant.melt$patient_id == "0891", "patient_id"] <- "P1"
malignant.melt[malignant.melt$patient_id == "6077", "patient_id"] <- "P2"
malignant.melt[malignant.melt$patient_id == "4230", "patient_id"] <- "P3"

# order
malignant.melt$patient_id <- factor(malignant.melt$patient_id, levels = c("P1", "P2", "P3"))

# 
malignant.gene.labels <- c("CEAMCAM6", "REG4", "CD9")

ggviolin(malignant.melt, x = "patient_id", y = "value", fill = "patient_id", 
          add = "jitter", legend = "none", facet.by = "variable", short.panel.labs = TRUE, add.params = list(size = 0.1, jitter = 0.2)) +
  rotate_x_text(angle = 45) +       # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     comparisons = my_comparisons) +
  theme_classic() + 
  labs(x = "Patients", y = "Normalized Expression") +
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 20, face = "bold" ),
        plot.title = element_text(colour = "black", size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "None")
ggsave(paste0(o.dir, "/", sample.tmp, "/malignant_gene_VlnPlot.pdf"), width = 10, height = 4.5)

############################################################################
##                  PLOT SOME MARKERS IN A DIFFERENT WAY
############################################################################

# gene set for scoring
genes.for.scoring <- list(c("AREG", "LGALS4", "CLDN4", "CEACAM6", "CEACAM5", "REG4", "CD9", "KLK6", "ZFAS1", "CST1", "AKR1C1",
                            "SERPINA1", "FTL", "LCN2", "ID3", "MUC13", "CDH17", "UCA1", "MSLN", "KLK8", "LY6D", "LY6E"))

# calcuate malignancy score
seurat.obj.big <- AddModuleScore(seurat.obj.big, features = genes.for.scoring, assay = "SCT", name = "Malignancy_Score")

# make a dataframe with the score and the identities from hashtags
plotting.df <- seurat.obj.big@meta.data[, c("Cell_Barcode", "orig.ident", "Malignancy_Score1")]

# rename the ids
plotting.df[plotting.df$orig.ident == "0891_WT", "orig.ident"] <- "P1WT"
plotting.df[plotting.df$orig.ident == "0891_PB2_EARLY", "orig.ident"] <- "P1C1_Early"
plotting.df[plotting.df$orig.ident == "0891_PB2_LATE", "orig.ident"] <- "P1C1_Late"
plotting.df[plotting.df$orig.ident == "0891_P31_EARLY", "orig.ident"] <- "P1C2_Early"
plotting.df[plotting.df$orig.ident == "0891_P31_LATE", "orig.ident"] <- "P1C2_Late"
plotting.df[plotting.df$orig.ident == "0891_P41_EARLY", "orig.ident"] <- "P1C3_Early"
plotting.df[plotting.df$orig.ident == "0891_P41_LATE", "orig.ident"] <- "P1C3_Late"

# rename the ids
plotting.df[plotting.df$orig.ident == "6077_WT", "orig.ident"] <- "P2WT"
plotting.df[plotting.df$orig.ident == "6077_C1_EARLY", "orig.ident"] <- "P2C1_Early"
plotting.df[plotting.df$orig.ident == "6077_C5_EARLY", "orig.ident"] <- "P2C2_Early"
plotting.df[plotting.df$orig.ident == "6077_C5_LATE", "orig.ident"] <- "P2C2_Late"
plotting.df[plotting.df$orig.ident == "6077_D3_EARLY", "orig.ident"] <- "P2C3_Early"
plotting.df[plotting.df$orig.ident == "6077_D3_LATE", "orig.ident"] <- "P2C3_Late"

# rename the ids
plotting.df[plotting.df$orig.ident == "4230_WT", "orig.ident"] <- "P3WT"
plotting.df[plotting.df$orig.ident == "4230_PB6_EARLY", "orig.ident"] <- "P3C1_Early"
plotting.df[plotting.df$orig.ident == "4230_PB6_LATE", "orig.ident"] <- "P3C1_Late"
plotting.df[plotting.df$orig.ident == "4230_APA4_EARLY", "orig.ident"] <- "P3C2_Early"
plotting.df[plotting.df$orig.ident == "4230_APA4_LATE", "orig.ident"] <- "P3C2_Late"
plotting.df[plotting.df$orig.ident == "4230_APA6_EARLY", "orig.ident"] <- "P3C3_Early"
plotting.df[plotting.df$orig.ident == "4230_APA6_LATE", "orig.ident"] <- "P3C3_Late"


# order
plotting.df$orig.ident <- factor(plotting.df$orig.ident, levels = c("P1WT", "P1C1_Early" , "P1C1_Late", "P1C2_Early", "P1C2_Late", "P1C3_Early", "P1C3_Late",
                                                                    "P2WT", "P2C1_Early", "P2C2_Early", "P2C2_Late", "P2C3_Early", "P2C3_Late",
                                                                    "P3WT", "P3C1_Early", "P3C1_Late", "P3C2_Early",  "P3C2_Late", "P3C3_Early", "P3C3_Late"))

# Visualize the expression profile
# my_comparisons <- list(c("P1", "P2"), c("P1", "P3"), c("P2", "P3"))
ggviolin(plotting.df, x = "orig.ident", y = "Malignancy_Score1", fill = "orig.ident", 
         add = "jitter", legend = "none", add.params = list(size = 0.1, jitter = 0.2)) +
  rotate_x_text(angle = 45) +       # Add global annova p-value
  #stat_compare_means(label = "p.signif", p.adjust.method = "bonferroni", method='t.test', 
   #                  comparisons = my_comparisons) +
  theme_classic() + 
  labs(x = "Clone IDs", y = "Malignancy Score") +
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold", angle = 45, hjust = 1),
        axis.title = element_text(colour = "black", size = 20, face = "bold"),
        plot.title = element_text(colour = "black", size = 20, face = "bold"),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "None")
ggsave(paste0(o.dir, "/", sample.tmp, "/MalignantScore_VlnPlot.pdf"), width = 10, height = 4.5)

############################################################################
##            CHECK THE DIFFERENCE BETWEEN EARLY AND LATE
############################################################################

# make a dataframe with the TFF genes
malignant.df <- as.data.frame(t(matrix[c("AREG", "LGALS4", "CLDN4", "CEACAM6", "CEACAM5", "REG4", "CD9", "KLK6", "ZFAS1", "CST1", "AKR1C1",
                                         "SERPINA1", "FTL", "LCN2", "ID3", "MUC13", "CDH17", "UCA1", "MSLN", "KLK8", "LY6D", "LY6E"),]))
malignant.df$barcode_id <- rownames(malignant.df)

# seurat metadata
metadata <- seurat.obj.big@meta.data[,c("orig.ident", "Cell_Barcode")]
metadata$barcode_id <- rownames(metadata)

# merge both information
malignant.meta.df <- merge(malignant.df, metadata, by = "barcode_id")

# make a dataframe with the "cancer" genes
malignant.meta.df$timepoint <- str_split_fixed(malignant.meta.df$orig.ident, "_", 3)[,3]
malignant.meta.df[grep("WT", malignant.meta.df$orig.ident), "timepoint"] <- "WT"

# rename the ids
malignant.meta.df[malignant.meta.df$orig.ident == "0891_WT", "orig.ident"] <- "P1WT"
malignant.meta.df[malignant.meta.df$orig.ident == "0891_PB2_EARLY", "orig.ident"] <- "P1C1_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "0891_PB2_LATE", "orig.ident"] <- "P1C1_Late"
malignant.meta.df[malignant.meta.df$orig.ident == "0891_P31_EARLY", "orig.ident"] <- "P1C2_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "0891_P31_LATE", "orig.ident"] <- "P1C2_Late"
malignant.meta.df[malignant.meta.df$orig.ident == "0891_P41_EARLY", "orig.ident"] <- "P1C3_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "0891_P41_LATE", "orig.ident"] <- "P1C3_Late"

# rename the ids
malignant.meta.df[malignant.meta.df$orig.ident == "6077_WT", "orig.ident"] <- "P2WT"
malignant.meta.df[malignant.meta.df$orig.ident == "6077_C1_EARLY", "orig.ident"] <- "P2C1_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "6077_C5_EARLY", "orig.ident"] <- "P2C2_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "6077_C5_LATE", "orig.ident"] <- "P2C2_Late"
malignant.meta.df[malignant.meta.df$orig.ident == "6077_D3_EARLY", "orig.ident"] <- "P2C3_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "6077_D3_LATE", "orig.ident"] <- "P2C3_Late"

# rename the ids
malignant.meta.df[malignant.meta.df$orig.ident == "4230_WT", "orig.ident"] <- "P3WT"
malignant.meta.df[malignant.meta.df$orig.ident == "4230_PB6_EARLY", "orig.ident"] <- "P3C1_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "4230_PB6_LATE", "orig.ident"] <- "P3C1_Late"
malignant.meta.df[malignant.meta.df$orig.ident == "4230_APA4_EARLY", "orig.ident"] <- "P3C2_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "4230_APA4_LATE", "orig.ident"] <- "P3C2_Late"
malignant.meta.df[malignant.meta.df$orig.ident == "4230_APA6_EARLY", "orig.ident"] <- "P3C3_Early"
malignant.meta.df[malignant.meta.df$orig.ident == "4230_APA6_LATE", "orig.ident"] <- "P3C3_Late"

# make patient variable
malignant.meta.df$patient_id <- str_split_fixed(malignant.meta.df$barcode_id, "_", 3)[,1]
malignant.meta.df[,c("barcode_id", "Cell_Barcode")] <- NULL
malignant.melt <- melt(malignant.meta.df, id.vars = c("orig.ident", "timepoint", "patient_id"))

# rename the ids
malignant.melt[malignant.melt$patient_id == "0891", "patient_id"] <- "P1"
malignant.melt[malignant.melt$patient_id == "6077", "patient_id"] <- "P2"
malignant.melt[malignant.melt$patient_id == "4230", "patient_id"] <- "P3"

# order
malignant.meta.df$orig.ident <- factor(malignant.meta.df$orig.ident, levels = c("P1WT", "P1C1_Early" , "P1C1_Late", "P1C2_Early", "P1C2_Late", "P1C3_Early", "P1C3_Late",
                                                                    "P2WT", "P2C1_Early", "P2C2_Early", "P2C2_Late", "P2C3_Early", "P2C3_Late",
                                                                    "P3WT", "P3C1_Early", "P3C1_Late", "P3C2_Early",  "P3C2_Late", "P3C3_Early", "P3C3_Late"))

# order
malignant.melt$patient_id <- factor(malignant.melt$patient_id, levels = c("P1", "P2", "P3"))

ggviolin(malignant.melt, x = "orig.ident", y = "value", fill = "orig.ident", 
         add = "jitter", legend = "none", facet.by = "variable", short.panel.labs = TRUE, add.params = list(size = 0.1, jitter = 0.2)) +
  rotate_x_text(angle = 45) +       # Add global annova p-value
  #stat_compare_means(label = "p.signif", method = "t.test",
   #                  comparisons = my_comparisons) +
  theme_classic() + 
  labs(x = "Patients", y = "Normalized Expression") +
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 20, face = "bold" ),
        plot.title = element_text(colour = "black", size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "None")
ggsave(paste0(o.dir, "/", sample.tmp, "/malignant_gene_VlnPlot.pdf"), width = 10, height = 4.5)

############################################################################
##            CHECK THE CONSERVED MARKERS ACROSS CLUSTERS
############################################################################

seurat.obj.big$pseudo_id <- "EL"
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat.obj.big,
                       ident.1 = cluster,
                       grouping.var = "pseudo_id",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:18), get_conserved)
write.table(conserved_markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_ConservedCluster_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

############################################################################
##    ANNOTATE THE EL ORGANOID DATA USING THE LITERATURE MARKER LIST
############################################################################

# use canonical markers to match the unbiased clustering to known cell types for Sathe dataset
new.cluster.ids <- c("Mucosal Cells", "Mucosal Cells", "Malignant Cells", "Mucosal Cells", "Goblet Cells", "MSCs",
                     "Mucosal Cells", "Stem Cells", "Mucosal Cells", "Mucosal Cells", "Malignant Cells", "Goblet Cells",
                     "Malignant Cells", "Pit Cells", "Neck-like Cells", "Enterocytes", "Malignant Cells", "Pit Cells", "Goblet Cells", 
                     "Proliferating Cells", "Proliferating Cells", "Malignant Cells")

# change cluster ids for Cell names
names(new.cluster.ids) <- levels(seurat.obj.big)
seurat.obj.big <- RenameIdents(seurat.obj.big, new.cluster.ids)

# store as cell type information
seurat.obj.big$CellType <- Idents(seurat.obj.big)

plot4 <- DimPlot(seurat.obj.big, reduction = "umap", pt.size = 1.5, label = F, label.size = 3) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=20, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Cell Types", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_CellTypes.pdf"), width = 6, height = 5)
print(plot4) 
dev.off()

# save final cell type annotated seurat_obj
saveRDS(seurat.obj.big, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_CellType_seurat_obj.rds"))

############################################################################
##                 PlOT UMAP REPRESENTATIONS HERE
###########################################################################

# get colours for the plot
cols <- viridis(option = "D", 20)

## CLONE ID
plot1 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "clone_id", pt.size = 1.5) +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("CLONE IDs", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_CloneID.pdf"), width = 6, height = 5)
print(plot1 + NoLegend()) 
dev.off()

plot2 <- DimPlot(seurat.obj.big, reduction = "umap", group.by = "timepoint", pt.size = 1.5) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=40, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=40, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  scale_color_manual(labels = c("Early", "Late", "WT"), values = c("#EA9F37", "#782867", "#4F9E4C")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 8)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Timepoint.pdf"), width = 6, height = 5)
print(plot2 + NoLegend()) 
dev.off()

plot3<- DimPlot(seurat.obj.big, reduction = "umap", pt.size = 1.5, label = T, label.size = 4) +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=20, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Cell Types", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Clusters.pdf"), width = 6, height = 5)
print(plot3 + NoLegend()) 
dev.off()

plot3<- DimPlot(seurat.obj.big, reduction = "umap", split.by = "old.ident", pt.size = 1.5, label = T, label.size = 4) +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=16, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=16, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=16)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Cell Types", override.aes = list(size = 6)))

# save single plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Split.pdf"), width = 12, height = 6.5)
print(plot3 + NoLegend())
dev.off()

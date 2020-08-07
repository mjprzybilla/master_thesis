#############################################################################################################################
##                                                                                                                      
##  Analysis of single-cell RNA-seq data from Hash-seq using Seurat
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
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "biomaRt", "scater", "DoubletFinder", "viridis")
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

## READ IN DATA OF INTEREST
#' The following code is adapted from the seurat vignette for the analysis of PBMCs provided [here](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html). 
#' The `Seurat` website is always the fountain of all knowledge, so if you have some issues with this code, having a look at the website might
#' help to figure out what is going on. First, we read in the 10x Genomics data matrix and look at some raw matrix QC metrics, as for instance,
#' counts per cell or number of genes per cell. 

# define the cellranger input directory of the sample of interest
# c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric/Sequencing5_cr3/Markdown/outs/"
# c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric/Sequencing13_6077_EL_cr3/Markdown/outs/"
# c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric/Sequencing19_4230_EL_Base/Markdown/outs/"
c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric/Sequencing27_0891EL/Markdown/outs/"

# define output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze"

# create o.dir
dir.create(o.dir)
setwd(o.dir)

# get sample.ids
sample.tmp <- unique(str_split_fixed(list.files(c.dir, full.names = T), "/", 8)[,6])

# which sample? 
message(sample.tmp)

# create sample-specific output directory with a "plots" subfolder
dir.create(paste0(o.dir, "/" , sample.tmp))
dir.create(paste0(o.dir, "/" , sample.tmp, "/plots"))

# read in genes x cells matrix from 10x genomics
c.mtx <- Read10X(paste0(c.dir, "/filtered_feature_bc_matrix"))

# read in metadata if present
# cell.metadata <- read_delim("/labs/ccurtis2/mjprzy/hashECB_data_updated/Sequencing5_4230_thaw/Hash/hash_outs/intermediate_output/Hash_correlation.txt", delim = "\t", col_names = T)
# cell.metadata <- read_delim("/labs/ccurtis2/mjprzy/hashECB_data_updated/Sequencing13_6077_all/Hash/hash_outs/intermediate_output/Hash_correlation.txt", delim = "\t", col_names = T)
# cell.metadata <- read_delim("/labs/ccurtis2/mjprzy/hashECB_data_updated/Sequencing19_4230_EL_Base/Hash/hash_outs/intermediate_output/Hash_correlation.txt", delim = "\t", col_names = T)
cell.metadata <- read_delim("/labs/ccurtis2/mjprzy/hashECB_data_updated/Sequencing27_0891EL_mnp05/Hash/hash_outs/intermediate_output/Hash_correlation.txt", delim = "\t", col_names = T)

# rename the first column
colnames(cell.metadata)[1] <- "Cell_Barcode"
cell.metadata$Cell_Barcode <- paste0(cell.metadata$Cell_Barcode, "-1")

# check dimensions
dim(c.mtx)

############################################################################
##                     CHECK OUT THE RAW MATRIX FIRST
############################################################################

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

# create seurat object with GO = gastric organoid
seurat.obj <- CreateSeuratObject(counts = c.mtx, project = paste0(sample.tmp, "_GO"), min.cells = 3, min.features = 200)

# add barcodes to metadta
seurat.obj@meta.data$Cell_Barcode <- rownames(seurat.obj@meta.data)

# get seurat metadata
seurat.metadata <- seurat.obj@meta.data

# adapt metadata
new.seurat.metadata <- merge(seurat.metadata, cell.metadata, by = "Cell_Barcode", all.x = T)
rownames(new.seurat.metadata) <- new.seurat.metadata$Cell_Barcode
seurat.obj@meta.data <- new.seurat.metadata

# Hashs
seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "Multiplet")

# for Seq 5
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "4230_APA4_080817")

# # for Seq13
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "6077_C1_LATE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "6077_D3_CDKN2A")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "6077_C1_LATE_BASE")

# # for Seq19
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "4230_WT_Passage_midi_BASE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "EARLY_4230_APA4_BASE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "LATE_4230_APA4_BASE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "EARLY_4230_PB6_BASE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "LATE_4230_PB6_BASE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "LATE_4230_APA6_BASE")
# seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "EARLY_4230_APA6_BASE")

# for Seq27
seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "0891_P41_W15")
seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "0891_PB2_W14")
seurat.obj <- subset(seurat.obj, subset = Sample_Origin != "0891_P31_W15")

# Calculate the percentage of mitochondrial genes
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")

# Calculate percent ribosomal genes
seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(seurat.obj, pattern = "^RP")

# the quality metrices in the seurat object are stored in the metadata
head(seurat.obj@meta.data, 5)

# make a FeatureScatter plot to visualize feature-feature relationships (can be used generally for all columns in metadata)
plot1 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# save plot 1 and plot2 together
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_FeatureScatterplots.pdf"), width = 20, height = 20)
print(plot1 + plot2)
dev.off()

# Load the the list of house keeping genes
hkgenes <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/housekeeping_genes_tirosh.txt")
hkgenes <- as.vector(hkgenes$V1)

# remove hkgenes that were not found
hkgenes.found <- which(toupper(rownames(seurat.obj@assays$RNA@counts)) %in% hkgenes)

# calculate the number of housekeeping genes per cell and add it as metadata
n.expressed.hkgenes <- Matrix::colSums(seurat.obj@assays$RNA@counts[hkgenes.found, ] > 0)
seurat.obj <- AddMetaData(object = seurat.obj, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

# plot additional QC
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_QC_plots.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "n.exp.hkgenes"), ncol = 5))
dev.off()

#############################################################################################################################

## REMOVE CELLS NOT MATCHING THE QC THRESHOLDS
#' Once this step is implemented, we should visually examine the quality of the dataset. Following this, we will
#' remove cells based on the calculated QC metrics. However, it is recommended to use very conservative thresholds,
#' as we do not want to remove valuable biological insights. 

# subset the seurat object to the QC-passing cells
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 500 & percent.mt < 20 & n.exp.hkgenes > 55 & percent.ribo < 40)

# plot post QC
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_postQC_plots.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "n.exp.hkgenes"), ncol = 5))
dev.off()

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
seurat.obj <- SCTransform(seurat.obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))

# perform cell cycle analysis (make sure to specify the "assay" parameter)
seurat.obj <- CellCycleScoring(seurat.obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_RidgePlot.pdf"), width = 20, height = 20)
RidgePlot(seurat.obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

# calculate the difference between G2M and S pahse scores and regress it out
seurat.obj$CC.Difference <- seurat.obj$S.Score - seurat.obj$G2M.Score

# view cell cycle scores and phase assignments
head(seurat.obj[[]])

# normalise again but this time including also the cell cycle scores
seurat.obj <- SCTransform(seurat.obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'CC.Difference'))

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
seurat.obj <- RunPCA(seurat.obj, features = c(s.genes, g2m.genes))

# save with regressed out cell cycle difference
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_CellCycle_Dimplot.pdf"), width = 20, height = 20)
print(DimPlot(seurat.obj, reduction = "pca"))
dev.off()

# save intermediate object
saveRDS(seurat.obj, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_SCtransform_seurat_obj.rds"))

#############################################################################################################################
# read in from sctransform step
# seurat.obj <- readRDS(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_SCtransform_seurat_obj.rds"))

## NORMALIZATION WITH RTAM
#' In addition, we will use a previously released normalization method called RTAM to perform normalization. Check out the Supplementary material
#' of the bioRxiv preprint to get further information about the normalization method <https://www.biorxiv.org/content/10.1101/2020.02.10.942607v1>.
#' In brief, this method relies on the assumption that highly covered genes are more reliable for the normalization of all cells in contrast to also
#' incorporate lowly covered genes. 

# run RTAM
message("Run RTAM...")

# run normalization based on raw count RNA assay
RTAM.data <- RTAM_normalization(mat = seurat.obj@assays$RNA@counts, method = "RTAM2", Min_nGn  = 500, Optimizing = FALSE)
rownames(RTAM.data) <- rownames(seurat.obj@assays$RNA@counts)
colnames(RTAM.data) <- colnames(seurat.obj@assays$RNA@counts)

# remove non-expressed genes
RTAM.matrix <- as.matrix(RTAM.data)
colnames(RTAM.matrix) <- colnames(RTAM.data)
final.RTAM.matrix <- as.matrix(RTAM.matrix[which(rowSums(RTAM.matrix != 0) >= 1  ), ])

# add the RTAM normalized count matrix to the seurat object
seurat.obj[["RTAM"]] <- CreateAssayObject(data = final.RTAM.matrix)

# switch default assay to RTAM
DefaultAssay(object = seurat.obj) <- "RTAM"

# Find variable features for the RTAM assay
seurat.obj <- FindVariableFeatures(object = seurat.obj, assay = "RTAM", nfeatures = 3000)

# scale data for the RTAM assay
seurat.obj <- ScaleData(seurat.obj, assay = "RTAM", features = seurat.obj@assays$RTAM@var.features)

# switch default back to SCT
DefaultAssay(object = seurat.obj) <- "SCT"

message("Finished..")
saveRDS(seurat.obj, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_RTAM_seurat_obj.rds"))

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

## SCT ASSAY
# run principal component analysis
seurat.obj <- RunPCA(seurat.obj, npcs = 100, features = seurat.obj@assays$SCT@var.features)

# perform Independent component analysis as well
seurat.obj <- RunICA(seurat.obj, features = seurat.obj@assays$SCT@var.features)

# examine and visualize the PCA results in different ways
print(seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)

# show PCA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_PCA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj, dims = 1:2, reduction = "pca"))
print(DimPlot(seurat.obj, reduction = "pca"))
dev.off()

# make dim heatmaps
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_PCA_heatmaps.pdf"), width = 20, height = 20)
print(DimHeatmap(seurat.obj, dims = 1, cells = 500, balanced = T))
print(DimHeatmap(seurat.obj, dims = 1:20, cells = 500, balanced = T))
dev.off()

# show ICA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_ICA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj, dims = 1:3, reduction = "ica"))
dev.off()

# Heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj, ndims = 100))
dev.off()

## RTAM
# switch default assay to RTAM
DefaultAssay(object = seurat.obj) <- "RTAM"

# run principal component analysis
seurat.obj <- RunPCA(seurat.obj, npcs = 100, assay = "RTAM", reduction.name = "pca_RTAM", reduction.key = "pca_RTAM_", features = seurat.obj@assays$RTAM@var.features)

# perform Independent component analysis as well
seurat.obj <- RunICA(seurat.obj, assay = "RTAM", reduction.name = "ica_RTAM", reduction.key = "ica_RTAM_", features = seurat.obj@assays$RTAM@var.features)

# examine and visualize the PCA results in different ways
print(seurat.obj[["pca_RTAM"]], dims = 1:5, nfeatures = 5)

# show PCA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_PCA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj, dims = 1:2, reduction = "pca_RTAM"))
print(DimPlot(seurat.obj, reduction = "pca_RTAM"))
dev.off()

# make dim heatmaps
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_PCA_heatmaps.pdf"), width = 20, height = 20)
print(DimHeatmap(seurat.obj, dims = 1, cells = 500, balanced = T, reduction = "pca_RTAM"))
print(DimHeatmap(seurat.obj, dims = 1:20, cells = 500, balanced = T, reduction = "pca_RTAM"))
dev.off()

# show ICA components with variable features
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_ICA_VizualizationPlots.pdf"), width = 20, height = 20)
print(VizDimLoadings(seurat.obj, dims = 1:3, reduction = "ica_RTAM"))
dev.off()

# Heuristic method - Elbow plot - ranking of pcs on the percentage of variance explained
# looking for an elbow, which should describe the number of pcs that describe the majority of true signal
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_ElbowPlot.pdf"), width = 20, height = 20)
print(ElbowPlot(seurat.obj, ndims = 100, reduction = "pca_RTAM"))
dev.off()

# switch default back to SCT
DefaultAssay(object = seurat.obj) <- "SCT"

#############################################################################################################################

## CELL CLUSTERING AND EMBEDDINGS
#' With the dimensioality reduction in place, we can apply graph-based clustering method where we construct a KNN graph based
#' on the euclidean distance in PCA space. We further use modularity optimization techniques such as the louvain algorithm. In brief,
#' the higher the resolution we set in this clustering, the greater the number of clusters generally is. A resolution of 0.4 - 1.2
#' typically give good results for 3k cells. Generally, the goal of cell clustering and embedding is to learn the underlying manifold
#' of the data in order to place similar cells together in a low-dimensional space. 

## SCT
# find NNs
seurat.obj<- FindNeighbors(seurat.obj, dims = 1:75, assay = "SCT", reduction = "pca")

# find clusters
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)

# run umap
seurat.obj <- RunUMAP(seurat.obj, dims = 1:75, assay = "SCT", reduction = "pca")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_UMAPplot.pdf"))
print(DimPlot(seurat.obj, reduction = "umap", label = T))
dev.off()

# run tsne
seurat.obj <- RunTSNE(seurat.obj, reduction = "pca", dims = 1:75)

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_TSNEplot.pdf"))
print(DimPlot(seurat.obj, reduction = "tsne", label = T))
dev.off()

## RTAM
# switch default assay to RTAM
DefaultAssay(object = seurat.obj) <- "RTAM"

# find NNs
seurat.obj<- FindNeighbors(seurat.obj, dims = 1:25, assay = "RTAM", reduction = "pca_RTAM")

# find clusters
seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)

# run umap
seurat.obj <- RunUMAP(seurat.obj, dims = 1:25, assay = "RTAM", reduction = "pca_RTAM", reduction.name = "umap_RTAM", reduction.key = "umapRTAM_")

# plot individual clusters
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_UMAPplot.pdf"))
print(DimPlot(seurat.obj, reduction = "umap_RTAM", label = T))
dev.off()

# run tsne
seurat.obj <- RunTSNE(seurat.obj, dims = 1:25, assay = "RTAM", reduction = "pca_RTAM", reduction.name = "tsne_RTAM", reduction.key = "tsneRTAM_")

# plot tsne
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_RTAM_TSNEplot.pdf"))
print(DimPlot(seurat.obj, reduction = "tsne_RTAM", label = T))
dev.off()

# switch default back to SCT
DefaultAssay(object = seurat.obj) <- "SCT"

# save intermediate object
saveRDS(seurat.obj, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_DimRed_seurat_obj.rds"))

#############################################################################################################################
# read seurat object
# seurat.obj <- readRDS(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_DimRed_seurat_obj.rds"))

## DOUBLET DETECTION WITH DOUBLETFINDER
#' Uses a fully-processed seurat object (after tsne has been run) needs PCs (range 1:10 for example).
#' pN - number of artifical generated doublets (default 25%)
#' pK - definition PC neighborhodd size
#' nExp - defines threshold for doublet/singlet predictions
#' pK Identification with sctransform used (no ground-truth - can be run with ground-truth as well)
#' Further information: https://github.com/chris-mcginnis-ucsf/DoubletFinder

# follow the DoubletFinder vignette for detecting and removing doublets
sweep.res.list.seurat.obj <- paramSweep_v3(seurat.obj, PCs = 1:75, sct = T)
sweep.stats_seurat.obj <- summarizeSweep(sweep.res.list.seurat.obj, GT = FALSE)

# plot Mean-variance normalized bimodality coefficient (bcmvn) 
# ground-truth-agnostic metric that coincides with pK
bcmvn.seurat.obj <- find.pK(sweep.stats_seurat.obj)

# plot pk value
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_pKvalue_plot.pdf"))
pK=as.numeric(as.character(bcmvn.seurat.obj$pK))
BCmetric=bcmvn.seurat.obj$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
dev.off()

# get annotatios for clustering
annotations <- seurat.obj@meta.data$seurat_clusters

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.075*length(seurat.obj@meta.data$barcode))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:75, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:75, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",pK_choose, "_", nExp_poi), sct = TRUE)

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_Doublet_plot.pdf"))
print(DimPlot(seurat.obj, reduction = "umap", group.by = paste0("DF.classifications_0.25_",pK_choose, "_", nExp_poi.adj)))
dev.off()

# assign a new column name to the classification column
colnames(seurat.obj@meta.data)[ncol(seurat.obj@meta.data)] <- "Doublet_Score"

# remove doublets from seurat object
seurat.obj <- subset(seurat.obj, subset = Doublet_Score != "Doublet")

# save interemediate objcet
saveRDS(seurat.obj, file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_seurat_obj.rds"))

#############################################################################################################################

## DIFFERENTIAL GENE EXPRESSION ANALYSIS
#' In this last part, we will create some plots and investigate distinct clusters for the expression of different genes, most
#' of which have been reported in the literature. We will use this manual investigation, to manually, but only roughly, 
#' assign possible identities to clusters. 

# # distinct Hash samples
# sample.tmp <- "Sequencing5_cr3"
# # sample.tmp <- "Sequencing13_6077_EL_cr3"
# # sample.tmp <- "Sequencing19_4230_EL_Base"
# # sample.tmp <- "Sequencing27_0891EL"
# 
# # # read seurat objects
# # seurat.obj <- readRDS(file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_seurat_obj.rds"))

# change the cluster identities from RTAM to SCT
Idents(seurat.obj) <- "SCT_snn_res.0.5"

## CELL MARKERS FROM ZHANG ET AL.
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_featurePlot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj, reduction = "umap", features = c("MUC5AC", "TFF1", "CD79A", "CD19", "OLFM4", "LGR5", "SOX2","CCKBR", "FABP1", "FABP2", "CA1", "VIL1", "MUC6", "TFF2", "CHGA", "TAC1", "TPH1", "CHGB"))
FeaturePlot(seurat.obj, reduction = "umap", features = c("CDK1", "MKI67", "CEACAM5", "CEACAM6", "PGC", "CXCL3", "IL8", "COL1A2", "LUM", "DCN", "PDPN", "FAP", "COL3A1", "COL6A1", "VWF", "ENG", "MCAM"))
FeaturePlot(seurat.obj, reduction = "umap", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "CD14", "CD68", "CSF1R", "MYL2", "ACTA2", "CD14", "CD68", "CSF1R", "MYL2"))
FeaturePlot(seurat.obj, reduction = "umap", features = c("ACTA2", "TPSAB1", "TPSB2", "PGA3", "PGA4", "LIPF", "CD2", "CD3D", "CD3E", "CD3G", "ATP4A", "ATP4B", "GAST", "GHRL", "SST"))
FeaturePlot(seurat.obj, reduction = "umap", features = c("OLFM4", "PHLDA1", "LEFTY1")) # stem cells
FeaturePlot(seurat.obj, reduction = "umap", features = c("CEACAM6", "BAX", "CCND2")) # cancer cells
FeaturePlot(seurat.obj, reduction = "umap", features = c("CEACAM5", "FABP1", "CDH17")) # non-specific cancer cells - also in enterocytes
dev.off()

## EPITHELIAL AND NON-EPITHELIAL
# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_Umap_Featureplot.pdf"), width = 20, height = 20)
FeaturePlot(seurat.obj, features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4"))
FeaturePlot(seurat.obj, features = c("CD4", "VIM", "ACTA2", "PTPRC"))
dev.off()

pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_EpithelialMarkers_expressionPlot.pdf"), width = 20, height = 20)
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("EPCAM", "KRT18", "MUC1", "KRT19", "CDH1", "CLDN4")))
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("CD4", "VIM", "ACTA2", "PTPRC")))
dev.off()

# plot probability distributions across clusters for the top5 genes indicated in the Supplementary data from Zhang et al.
# However, markers shared by different cell types were removed from both.
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_exp_violin_plot.pdf"), width = 20, height = 20)
# check for PMCs
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("GKN1", "GKN2", "MUC5AC", "TFF1", "DPCR1")) + labs(title = "PMCs"))
# check for MSCs
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("OLFM4", "RPL7", "CLDN4","TSPAN8", "REG1A")) + labs(title = "MSCs"))
# check for enterocytes
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("FABP1", "FABP2", "RBP2", "ANPEP", "APOA4")) + labs(title = "Enterocytes"))
# check for GMC
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("PRR4", "C6orf58","MUC6", "TFF2", "LTF")) + labs(title = "GMCs"))
# check for enteroendocrine
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("CHGA", "PCSK1N", "SCG5", "CHGB", "TPH1")) + labs(title = "Enteroendocrine")) 
# check for PCs
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("TOP2A", "MKI67", "UBE2C", "HMGB2", "PTTG1")) + labs(title = "PCs"))
# check for cancer cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("REG4","CLDN7","KRT18", "LGALS3", "CEACAM6")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("CLDN3", "CST1", "MUC3A", "CLDN4", "PI3", "UBD")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("CDH17", "PRAP1", "UBE2C", "CCL20", "LCN2", "SERPINB5")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("RRM2", "MYBL2", "MMP7", "TPX2", "MISP", "TMPRSS4")) + labs(title = "Cancer Cells"))
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("RMRP", "CLDN1", "GPRC5A", "CLRN3", "CXCL1", "MSLN")) + labs(title = "Cancer Cells"))
# check for neck like cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("MIA", "CXCL3", "CXCL2", "CXCL17", "CLU")) + labs(title = "Neck-like cells"))
# check for Goblet cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("SPINK4", "TFF3", "MUC2", "ITLN1", "ZG16")) + labs(title = "Goblet Cells"))
# check for Chief cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("PGA3", "PGA4", "LIPF", "CHIA", "PDIA2")) + labs(title = "Chief Cells")) 
# check for parietal cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("ATP4A", "ATP4B", "VEGFB")) + labs(title = "Parietal cells"))
# check for endocrine cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("GAST", "GHRL", "SST")) + labs(title = "Endocrine cells"))
##GAO DATASET
# check for Chief Stem Cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("GIF", "LRIG1", "PROCR")) + labs(title = "Chief Stem cells"))
# check for Enteroendocrine
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("SST", "PYY", "PROX1", "TPH1", "REG4", "NEUROD1", "GIP")) + labs(title = "Enteroendocrine cells"))
# check for Enterocytes
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("CA2", "KRT20", "ALPI", "TREH", "LCT", "MME", "CDH1", "VIL1", "AQP8", "SI", "CDX2")) + labs(title = "Enterocytes"))
# check for paneth cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("NOTUM", "NUPR1", "PLA2G2A")) + labs(title = "Paneth Cells"))
# check for tuft cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("DCLKK1", "TRPM5", "PTGS1", "RGS13")) + labs(title = "Tuft Cells"))
# check for Stem cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("PROM1", "MEX3A", "TERT", "LGR5", "SOX9", "CD24", "ALCAM", "PROCR")) + labs(title = "Stem Cells"))
# check for PROCRhigh progenitor cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("PPIH", "SGK1", "GPT2")) + labs(title = "PROCR High Progenitor Cells"))
# check for Neck progenitor cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("TFF2", "ECE1", "PROM1")) + labs(title = "Neck-Progenitor Cells"))
# check for HES1high progenitor cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("LINC01578", "PLXDC1")) + labs(title = "HES1 High Progenitor Cells"))
# check for parietal progenitor cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("ATP4A", "ATP4B", "CFH", "FGD5")) + labs(title = "Parietal Progenitor Cells"))
# check for Pit progenitor cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("LRIG1", "AXIN2", "CD44", "ACTC1")) + labs(title = "Pit Progenitor Cells"))
# check for Pit Cells
print(VlnPlot(object = seurat.obj, assay = "SCT", features = c("ID1", "MBD1", "GREM", "RSPO2")) + labs(title = "Pit Cells"))
dev.off()

###########################################################################
#               ANALYSE THE GENE EXPRESSION PER CLUSTER
###########################################################################

# find markers for every cluster compared to all remaining cells, report only the positive ones
seurat.obj.markers <- FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_seurat.obj.markers <- seurat.obj.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top10_seurat.obj.markers <- seurat.obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(seurat.obj.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_cluster_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# seurat.obj.markers <- read.table(paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_cluster_markers.txt"), header = T)

## TOP MARKERS
# top 5 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj, features = top5_seurat.obj.markers$gene, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")))
dev.off()

pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_dataSlot_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj, assay = "SCT", slot = "data", features = top5_seurat.obj.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue")))
dev.off()

###########################################################################
#             ANALYSE THE GENE EXPRESSION PER CLONE
###########################################################################
# # Seq5
# # add timepoint metadata
# seurat.obj@meta.data[seurat.obj@meta.data == "4230_APA4_052617"] <- "TP1"
# seurat.obj@meta.data[seurat.obj@meta.data == "4230_APA4_103117"] <- "TP2"
# seurat.obj@meta.data[seurat.obj@meta.data == "4230_APA4_012318"] <- "TP3"
# seurat.obj@meta.data[seurat.obj@meta.data == "4230_APA4_041718"] <- "TP4"

# set idents to clone id
Idents(seurat.obj) <- "Sample_Origin" # EL single

# Seq5
# levels(seurat.obj) <- c("4230_WT", "TP1", "TP2", "TP3", "TP4")

# Seq13
# levels(seurat.obj) <- c("6077_WT", "6077_C1_EARLY", "6077_C5_EARLY", "6077_C5_LATE", "6077_D3_EARLY", "6077_D3_LATE")

# Seq19
# levels(seurat.obj) <- c("4230_WT_Passage_midi", "EARLY_4230_PB6", "LATE_4230_PB6", "EARLY_4230_APA4", "LATE_4230_APA4", "EARLY_4230_APA6", "LATE_4230_APA6")

# Seq27
levels(seurat.obj) <- c("0891_WT_P2", "0891_P_B2_EARLY_102219", "0891_P_B2_LATE_102219", "0891_P_31_EARLY", "0891_P_31_LATE", "0891_P_41_EARLY", "0891_P_41_LATE")

clone.seurat.obj.markers <- FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_clone_seurat.obj.markers <- clone.seurat.obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# write table with the DGEs per cluster
write.table(clone.seurat.obj.markers, paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_clone_markers.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

## TOP MARKERS
# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj, features = top10_clone_seurat.obj.markers$gene, size = 3, angle = 45, disp.min = -3, disp.max = 3) + scale_fill_gradientn(colors = c("blue", "lightgrey", "red")))
dev.off()

# top 3 markers (or all markers if less than 10) for each cluster.
pdf(paste(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_SCT_dataSlot_cloneID_ExpressionHeatmap", ".pdf", sep = ""), width = 10.5, height = 9)
print(DoHeatmap(seurat.obj, assay = "SCT", slot = "data", features = top10_clone_seurat.obj.markers$gene, size = 3, angle = 45, disp.max = 3) + scale_fill_gradientn(colors = c("white", "blue")))
dev.off()

###########################################################################
#            GENERATE PSEUDOBULKS PER CLONE FOR EACH GENE
###########################################################################

# ## SEQ5 APA4 TIMECOURSE
# # calculate average gene expression per timepoint for Aziz
# clone.averages <- AverageExpression(seurat.obj)
# SCT.pseudobulk <- clone.averages$SCT
# 
# write.table(SCT.pseudobulk, "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/Sequencing5_cr3/Sequencing5_cr3_pseudobulk.txt", col.names = T, row.names = T, sep = "\t", quote = F)

############################################################################
##        MAKE PLOTS FOR THE UMAP EMBEDDING FOR EL ORGANOIDS
############################################################################

# add timepoint metadata
## SEQ13
# seurat.obj$timepoint <- str_split_fixed(seurat.obj$Sample_Origin, "_", 3)[,3]
# seurat.obj@meta.data[grep("WT", seurat.obj$Sample_Origin), "timepoint"] <- "WT"

## SEQ19
# seurat.obj$timepoint <- str_split_fixed(seurat.obj$Sample_Origin, "_", 3)[,1]
# seurat.obj@meta.data[grep("4230", seurat.obj$Sample_Origin), "timepoint"] <- "WT"

## SEQ27
seurat.obj$timepoint <- str_split_fixed(seurat.obj$Sample_Origin, "_", 5)[,4]
seurat.obj@meta.data[grep("WT", seurat.obj$Sample_Origin), "timepoint"] <- "WT"

# add clone id metadata
## SEQ13
# seurat.obj$clone_id <- paste0(str_split_fixed(seurat.obj$Sample_Origin, "_", 3)[,1], "_",str_split_fixed(seurat.obj$Sample_Origin, "_", 3)[,2])

##SEQ19 & SEQ27
seurat.obj$clone_id <- paste0(str_split_fixed(seurat.obj$Sample_Origin, "_", 4)[,1], "_",str_split_fixed(seurat.obj$Sample_Origin, "_", 4)[,2],  "_",str_split_fixed(seurat.obj$Sample_Origin, "_", 4)[,3])

## CLONE ID
plot1 <- DimPlot(seurat.obj, reduction = "umap", group.by = "clone_id", pt.size = 2) +
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

plot2 <- DimPlot(seurat.obj, reduction = "umap", group.by = "timepoint", pt.size = 2) +
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

plot3 <- DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", pt.size = 2, label = T, label.size = 6) +
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

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest1_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj, reduction = "umap", features = c("MUC5AC", "TFF1", "CEACAM6", "PGC"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest2_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj, reduction = "umap", features = c("WFDC2", "MUC5B", "REG4", "OLFM4"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest3_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj, reduction = "umap", features = c("TSPAN8", "VIL1", "AKR1C1", "LYZ"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GenesOfInterest4_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj, reduction = "umap", features = c("DMBT1", "FABP1", "GSTA1", "PSCA"), pt.size = 1, ncol = 2)
dev.off()

# plot feature expression on a tSNE or PCA plot
pdf(paste0(o.dir, "/" , sample.tmp, "/plots/", sample.tmp, "_UMAP_GeneHeterogeneity_FeaturePlot.pdf"), width = 11, height = 9)
FeaturePlot(seurat.obj, reduction = "umap", features = c("MUC5AC", "TFF1", "CEACAM6", "WFDC2", "REG4", "DMBT1", "FABP1", "GSTA1", "LYZ"), pt.size = 1, ncol = 3)
dev.off()


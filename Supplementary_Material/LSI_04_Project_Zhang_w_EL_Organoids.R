#################################################################################################################################
##                                                                                                                      
##  Perform clustering and umap embedding on 10x single-cell RNA-sequencing from Early and Late organoids with Zhang et al. UMAP
##  Date: 16 May 2020                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##  Credit: This code was designed based on code developed by Jeffrey Granja in frame of Granja*, Klemm*, Mcginnis* et al. 
##          A single cell framework for multi-omic analysis of disease identifies  malignant regulatory signatures in 
##          mixed phenotype acute leukemia (2019, https://github.com/GreenleafLab/MPAL-Single-Cell-2019). 
##
##                                                                                                                      
#################################################################################################################################
# clear workspace
rm(list=ls())
set.seed(1) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", "dendextend", "tidyverse", 
                      "Matrix", "ComplexHeatmap", "Rtsne", "robustbase", "psych", "cluster", "factoextra", "Matrix.utils", "Repitools",
                      "BSgenome.Hsapiens.UCSC.hg38", "Seurat", "GenomicRanges", "SingleCellExperiment", "matrixStats", "readr",
                      "magrittr", "edgeR", "uwot", "BiocManager", "Rcpp", "biomaRt", "dplyr", "viridis", "FNN", "httr", "Matrix", "ArchR","reshape2",
                      "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# FUNCTIONS
#####################################################################################

#' **Binarize a count matrix**
#' Use this function to convert a count matrix into a binary matrix. 
#' binarizeMat(
#'      mat = count.matrix 
#' )
#' 
#' 

#Binarize Sparse Matrix
binarizeMat <- function(mat){
  mat@x[mat@x > 0] <- 1
  mat
}

#' **Calculate Latent semantic indexing**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' calcLSI(
#'      mat = count.matrix,
#'      nComponents = nPCs, 
#'      binarize = TRUE, 
#'      nFeatures = NULL
#' )
#' 

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  
  # binarize the matrix 
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  
  # order the matrix according to the highly represented genes
  # calculate the number of counts per gene
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
    mat <- mat[idx,] 
  }else{
    idx <- which(Matrix::rowSums(mat) > 0)
    mat <- mat[idx,]
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat) # library size per cell
  rowSm <- Matrix::rowSums(mat) # counts per gene 
  
  #Calc TF IDF
  # Term frequency and inverse document frequency
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm) # calculate the term frequency 
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector") # inverse document frequency
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs # matrix multiplication %*%
  
  #Calc SVD then LSI 
  # SVD = Singular value decomposition (alternative to PCA)
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object with all variables
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    idx = idx, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    date = Sys.Date(),
    seed = 1)
  
  out
  
}

projectLSI <- function(mat, lsi){   
  
  #Get Same Features
  mat <- mat[lsi$idx,]
  if(lsi$binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1       
  }
  
  #Calc TF IDF
  rowsToZero <- which(lsi$rowSm == 0)
  setToZero <- which((mat@i + 1) %in% rowsToZero)
  if(length(setToZero) > 0){
    mat@x[setToZero] <- 0
  }
  
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
    tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
  }
  
  #Calc V
  V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
  diag(svdDiag) <- lsi$svd$d
  matSVD <- t(svdDiag %*% t(V))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  return(matSVD)
  
}

#' **Find nearest neighbors with help of seurat**
#' Use this function to perform PCA on the SVD data matrix and
#' further find the nearest neigbors based on that. From this, communities
#' of cells can be defined as cluster
#' seuratSNN(
#'      matSVD = matSVD (input is a dimensionality reduced matrix) 
#'      dims.use = 1:50 (these are the dimensions used for PCA)
#' )
#' 
#' 
#Clustering function using seurat SNN
seuratSNN <- function(matSVD, dims.use = 1:50, ...){
  set.seed(1)
  message("Making Seurat Object...")
  
  # get original order of rows
  rn <- rownames(matSVD)
  
  # subset the matrix to three
  tmp <- matrix(rnorm(nrow(matSVD) * 3, 10), ncol = nrow(matSVD), nrow = 3)
  colnames(tmp) <- rownames(matSVD)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
  
  # create seurat object and runn PCA on the SVD matrix
  obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
  obj[['pca']] <- Seurat::CreateDimReducObject(embeddings=matSVD, key='PC_', assay='RNA')
  clustParams$object <- obj
  clustParams$reduction <- "pca"
  clustParams$dims <- seq_len(ncol(matSVD))
  
  obj <- suppressWarnings(do.call(Seurat::FindNeighbors, clustParams))
  clustParams$object <- obj
  
  obj <- suppressWarnings(do.call(Seurat::FindClusters, clustParams))
  
  #Get Output
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  clust <- paste0("Cluster",match(clust, unique(clust)))
  names(clust) <- rownames(matSVD)
  clust <- clust[rn]
  
  clust
  
  
}

#' **Generate a sparse matrix from a fragment file**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' countInsertions(
#'      windows = windows (GRange object resulting from makeWindows function), 
#'      fragments = , 
#'      by = "RG"
#' )
#' 

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

#Compute Fast Sparse Row Variances
sparseRowVariances <- function (m){
  rM <- Matrix::rowMeans(m)
  rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  return(rV)
}

#' **Generate a sparse matrix from a fragment file**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' countInsertions(
#'      windows = windows (GRange object resulting from makeWindows function), 
#'      fragments = , 
#'      by = "RG"
#' )
#' 

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

#' **Generate a sparse matrix from a fragment file**
#' Use this function to create a GRange object that contains chr, start, end, as well as 
#' information about the GC and AT content of each window.
#' countInsertions(
#'      windows = windows (GRange object resulting from makeWindows function), 
#'      fragments = , 
#'      by = "RG"
#' )
#' 

#Optimized LSI for scRNA-seq analysis
optimizeLSI <- function(mat, scaleTo = 10000, priorCount = 3, pcsUse = 1:25, 
                        resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){
  
  set.seed(seed)
  stopifnot(length(resolution) > 1)
  
  #Initialize List
  lsiOut <- list()
  
  #Initial LSI uses variances that are across all single cells and will have larger batch relationships
  i <- 1
  message("Initial LSI...")
  matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
  matNorm@x <- log2(matNorm@x + 1)
  idVarFeatures <- head(order(sparseRowVariances(matNorm),decreasing=TRUE), varFeatures[i])
  lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
  clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)
  
  #Store
  lsiOut[[paste0("iter", i)]] <- list(
    lsiMat = lsiObj$matSVD, 
    varFeatures = idVarFeatures, 
    clusters = clusters
  )
  
  for(i in seq(2, length(varFeatures))){
    
    message(sprintf("Additional LSI %s...", i))
    
    #Run LSI
    clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE), log=TRUE, prior.count = priorCount)
    idVarFeatures <- head(order(rowVars(clusterMat), decreasing=TRUE), varFeatures[i])
    lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
    clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)
    
    if(i == length(varFeatures)){
      #Save All Information from LSI Attempt
      lsiOut[[paste0("iter", i)]] <- list(
        lsiObj = lsiObj, 
        varFeatures = idVarFeatures, 
        clusters = clusters,
        matNorm = matNorm
      )
    }else{
      lsiOut[[paste0("iter", i)]] <- list(
        lsiMat = lsiObj$matSVD, 
        varFeatures = idVarFeatures, 
        clusters = clusters
      )
    }
    
  }
  
  return(lsiOut)
  
}

sparseMatTTest <- function(mat1, mat2, m0 = 0){
  #Get Population Values
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  n <- n1 + n2
  #Sparse Row Means
  m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
  m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
  #Sparse Row Variances
  v1 <- computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
  v2 <- computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
  #Calculate T Statistic
  se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
  tstat <- (m1-m2-m0)/se
  #tstat <- sqrt((n1 * n2) / n) / sqrt((n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2)
  pvalue <- 2*pt(-abs(tstat), n - 2)
  fdr <- p.adjust(pvalue, method = "fdr")
  out <- data.frame(fdr = fdr, pval = pvalue, tstat = tstat, mean1 = m1, mean2 = m2, var1 = v1, var2 = v2, n1 = n1, n2 = n2)
  return(out)
}

`%notin%` <- Negate(`%in%`)

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}
######################################################################################
# SET WORKING DIRECTORY AND LOAD BATCH CORRECTED (COMBAT) ORGANOID DATA
######################################################################################
# working directory
w.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/"

# output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/freeze"
dir.create(o.dir)
setwd(o.dir)

# read in the COMBAT corrected organoid data
EL.bc.matrix <- read.table("/labs/ccurtis2/mjprzy/scRNA_analysis/combat/freeze/EL_Seq13_19_27/EL_Seq13_19_27_organoid_bc_data.txt", header = T, sep = "\t") # EL data
colnames(EL.bc.matrix) <- gsub("X", "", colnames(EL.bc.matrix))
EL.bc.matrix <- as.matrix(EL.bc.matrix)
EL.data <- as(EL.bc.matrix, "sparseMatrix")

# read in Seurat object to get the corresponding metadata
seurat.obj.big <- readRDS(paste0(w.dir, "EL_Seq13_19_27/EL_Seq13_19_27_seurat_obj.rds")) # EL data
seurat.metadata <- seurat.obj.big@meta.data

######################################################################################
# LOAD THE ZHANG REFERENCE DATASET THAT IS USED FOR THE PROJECTION
######################################################################################
## ZHANG DATASET
# set sample name
sample.tmp <- "EL_Seq13_19_27_Zhang"
dir.create(sample.tmp)

# read in the seurat object for the annotated Sathe dataset which was generated 
# within the scRNA_01 script
seurat.obj.big <- readRDS(paste0("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_ZhangEtal_freeze/merged_data_ZhangEtal_freeze_iterativeLSI_seurat.obj.rds"))
reference.data <- seurat.obj.big@assays$RNA@counts

######################################################################################
# CONVERT THE REFERENCE DATASET OF INTEREST INTO A SUMMARIZED EXPERIMENT OBJECT
######################################################################################

#Identify Gene Universe
gU <- intersect(rownames(EL.data), rownames(reference.data))
gU <- gU[!grepl("^MT", gU)]

# subset the EL data to the same gene universe
EL.data <- EL.data[gU,]

# create the Summarized Experiment from the normalized data
sce <- SingleCellExperiment(list(counts=EL.data),
                            colData=DataFrame(Cell_barcode = colnames(EL.data)),
                            rowData=DataFrame(hgnc_symbol = rownames(EL.data)),
                            metadata=seurat.metadata
)

# add the patient name here
colData(sce)$Group <- str_split_fixed(rownames(colData(sce)), "_", 2)[,1]

# check Summarized Experiment object
sce

## EL Organoids
# class: SingleCellExperiment 
# dim: 15618 26145 
# metadata(25): Cell_Barcode orig.ident ... pANN_0.25_0.3_0 Doublet_Score
# assays(1): counts
# rownames(15618): LINC00115 FAM41C ... CLDN14 LINC00114
# rowData names(1): hgnc_symbol
# colnames(26145): 6077_EL_AAACCCAAGACATCAA 6077_EL_AAACCCAAGACGGAAA ... 0891_EL_TTTGTTGTCTAGCCTC 0891_EL_TTTGTTGTCTCAACGA
# colData names(2): Cell_barcode Group
# reducedDimNames(0):
#   spikeNames(0):
#   altExpNames(0):

######################################################################################
# CLUSTERING ANALYSIS OF THE FULL ORGANOID DATASET STARTS HERE
######################################################################################
# Set Clustering Parameters
varGenesToUse <- c(1000, 1000, 1000) #Choose a higher number of variable peaks
resolution <- c(0.2,0.4,0.6) #Clustering resolutions for Seurat SNN
verbose = TRUE
tstart <- NULL

# set clustering parameters for seurats find cluster method
clustParams <- list()
clustParams$verbose <- verbose
clustParams$tstart <- tstart

#Optimize LSI Features for STP samples
matAll <- assay(sce)
lsiObj <- optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse,  pcsUse = 1:25)

#UMAP
set.seed(1)
umap <- uwot::umap(
  lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25], 
  n_neighbors = 35, 
  min_dist = 0.5, 
  metric = "euclidean", 
  n_threads = 5, 
  verbose = TRUE, 
  ret_model = FALSE
)

######################################################################################
# CLASSIFICATION OF EARLY LATE SAMPLES FROM GASTRIC ORGANOIDS
######################################################################################

# add timepoint metadata
sce@metadata$timepoint <- str_split_fixed(sce@metadata$HashTag, "_", 3)[,3]
sce@metadata$timepoint[grep("WT", sce@metadata$HashTag)] <- "WT"

# calculate the percentage of disease cells per cluster Plot Info
cells <- sce@metadata$timepoint

# split the cells by the clusters 
splitCells <- split(cells,lsiObj[[length(lsiObj)]]$clusters)

# calculate the proporationn of LATE cells in each cluster
df <- data.frame(
  clusters = names(splitCells),
  proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="LATE") / length(splitCells[[x]])))
)

#Plot UMAP Data Frame
plotDF <- data.frame(umap)
rownames(plotDF) <- c(colnames(sce))

# add metadata to the plot data frame
plotDF$type <- cells
plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
plotDF$patient_id <- str_split_fixed(rownames(plotDF), "_", 2)[,1]

# setup a projection plot directory
plotDir <- paste0(sample.tmp, "/classification/")
dir.create(plotDir,recursive=TRUE)

######################################################################################
# PLOT THE ORGANOID DATA COLOURED ACCORDING TO TIMEPOINTS
######################################################################################
#Plot PDFs
pdf(paste0(plotDir, sample.tmp,"_WT-EL-Timepoint-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(plotDF, aes(X1,X2,color=type)) + 
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("WT"="lightgrey","EARLY"="dodgerblue3","LATE"="firebrick3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Timepoints", override.aes = list(size = 6)))
dev.off()

######################################################################################
# PLOT THE ORGANOID DATA COLOURED ACCORDING TO CLUSTERS
######################################################################################

pdf(paste0(plotDir, sample.tmp,"_WT-EL-Clusters-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(plotDF, aes(X1,X2,color=clusters)) + 
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_viridis_d() + 
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=20, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Clusters", override.aes = list(size = 4)))
dev.off()

######################################################################################
# PLOT THE ORGANOID DATA COLOURED ACCORDING TO PATIENTS
######################################################################################

pdf(paste0(plotDir,sample.tmp,"_WT-EL-Patients-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(plotDF, aes(X1,X2,color=patient_id)) + 
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_viridis_d() + 
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

####################################################
# LOAD THE ZHANG DATASET FROM STEP 1
####################################################

# Here, we read in the data and uwot object from scRNA_01_Clustering_UMAP_EL_Organoids.R script
se <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_ZhangEtal_freeze/SummarizedExperiment_merged_data_ZhangEtal_freeze.rds")

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_ZhangEtal/SummarizedExperiment_merged_data_ZhangEtal_UMAP-model.uwot")

####################################################
# PERFORM PROJECTION INTO LSI UMAP
####################################################

# subset the organoid LSI Projection Matrix to the variableGenes from the projection
lsiGenes <- metadata(se)$variableGenes
matProjectLSI <- assay(sce[lsiGenes,])

# LSI Project
lsiReference <- metadata(se)$optimizeLSI[[length(metadata(se)$optimizeLSI)]]$lsiObj

# run lsi projection
lsiProjection <- projectLSI(matProjectLSI, lsiReference)

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:25], umapManifold, verbose = TRUE)

######################################################################################
# CREATE A DATAFRAME WITH THE UMAP COORDINATES AND METADATA
######################################################################################

#Plot Projection
refDF <- data.frame(row.names = colnames(se), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "background")

# ZHANG DATASET
refDF$sample_type <- str_split_fixed(rownames(refDF), "_", 2)[,1]
refDF$sample_type[grep("CAG", refDF$sample_type)] <- "CAG"
refDF$sample_type[grep("NAG", refDF$sample_type)] <- "NAG"
refDF$sample_type[grep("EGC", refDF$sample_type)] <- "EGC"
refDF$sample_type[grep("IMS", refDF$sample_type)] <- "IMS"
refDF$sample_type[grep("IMW", refDF$sample_type)] <- "IMW"

# set up the EL sample dataframe which shall be projected
proDF <- data.frame(row.names = colnames(sce), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = plotDF[colnames(sce),]$type)
proDF$sample_type <- str_split_fixed(rownames(proDF), "_", 2)[,1]
proDF$Type <- as.character(proDF$Type)

# combine both DF to get a dataframe for the projection
projectionDF <- rbind(refDF, proDF)

######################################################################################
# PLOT UMAP OF PROJECTION WITH COLORING ACCORDING TO ALL TIMEPOINTS
######################################################################################

pdf(paste0(plotDir, sample.tmp,"_WT-EL-Projection-UMAP-Timepoints.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(projectionDF, aes(X1,X2,color=Type)) + 
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("background"="lightgrey", "WT" = "chartreuse4","EARLY"="dodgerblue3", "LATE"="firebrick3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

######################################################################################
# PLOT UMAP OF PROJECTION WITH COLORING ACCORDING TO SAMPLE TYPES IN ZHANG DATASET
######################################################################################

pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Samples.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(projectionDF, aes(X1,X2,color=sample_type)) + 
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("CAG"="lightgrey", "NAG"="lightgrey", "EGC"="deeppink2", "IMS"="aliceblue", "IMW"="aliceblue",
                              "6077" = "chartreuse4", "4230"="dodgerblue3","0891"="firebrick3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

######################################################################################
# PLOT PROJECTION ACCORDING TO INDIVIDUAL TIMEPOINTS
######################################################################################

# get WT only dataframe
WT.proDF <- proDF[proDF$Type == "WT",]
WT.projectionDF <- rbind(refDF, WT.proDF)

## WT
pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Timepoints-WT.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(WT.projectionDF, aes(X1,X2,color=Type)) +
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("background"="lightgrey", "WT" = "chartreuse4")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) +
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 6)))
dev.off()

# get Early only dataframe
Early.proDF <- proDF[proDF$Type == "EARLY",]
Early.projectionDF <- rbind(refDF, Early.proDF)

## EARLY
pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Timepoints-EARLY.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(Early.projectionDF, aes(X1,X2,color=Type)) +
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("background"="lightgrey", "EARLY"="dodgerblue3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) +
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 6)))
dev.off()

# get WT only dataframe
Late.proDF <- proDF[proDF$Type == "LATE",]
Late.projectionDF <- rbind(refDF, Late.proDF)

## LATE
pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Timepoints-LATE.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(Late.projectionDF, aes(X1,X2,color=Type)) +
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("background"="lightgrey", "LATE"="firebrick3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) +
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Timepoint", override.aes = list(size = 6)))
dev.off()

######################################################################################
# PLOT PROJECTION ACCORDING TO INDIVIDUAL PATIENTS - ZHANG
######################################################################################

# get 0891 only dataframe
P1.proDF <- proDF[proDF$sample_type == "0891",]
P1.projectionDF <- rbind(refDF, P1.proDF)

## 0891
pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Patients-0891.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(P1.projectionDF, aes(X1,X2,color=sample_type)) +
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("CAG" = "lightgrey","NAG"="lightgrey", "EGC"="lightgrey","IMS"="lightgrey", "IMW"="lightgrey", "0891" = "chartreuse4")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) +
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

# get WT only dataframe
P2.proDF <- proDF[proDF$sample_type == "6077",]
P2.projectionDF <- rbind(refDF, P2.proDF)

## 6077
pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Patients-6077.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(P2.projectionDF, aes(X1,X2,color=sample_type)) +
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("CAG" = "lightgrey","NAG"="lightgrey", "EGC"="lightgrey","IMS"="lightgrey", "IMW"="lightgrey", "6077" = "dodgerblue3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) +
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

# get WT only dataframe
P3.proDF <- proDF[proDF$sample_type == "4230",]
P3.projectionDF <- rbind(refDF, P3.proDF)

## 4230
pdf(paste0(plotDir,sample.tmp,"_WT-EL-Projection-UMAP-Patients-4230.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(P3.projectionDF, aes(X1,X2,color=sample_type)) +
  geom_point() +
  theme_bw() +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("CAG" = "lightgrey","NAG"="lightgrey", "EGC"="lightgrey","IMS"="lightgrey", "IMW"="lightgrey", "4230" = "firebrick3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) +
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

######################################################################################
#                   ANALYSIS OF INDIVIDUAL SAMPLES HERE
######################################################################################
# get the sample ids
sample.ids <- unique(metadata(sce)$HashTag)
#sample.ids <- unique(paste0(str_split_fixed(metadata(sce)$HashTag, "_", 3)[,1], "_", str_split_fixed(metadata(sce)$HashTag, "_", 3)[,2]))

# alternatively, use patient ids
patient.ids <- unique(str_split_fixed(metadata(sce)$HashTag, "_", 3)[,1])

###############################################################################
# LOAD THE SATHE SEURAT OBJECT WITH THE CELL TYPE ANNOTATION
###############################################################################

# load the respective seurat object for the reference with the assigned cell types
seurat.obj.big <- readRDS(paste0("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_ZhangEtal_freeze/merged_data_ZhangEtal_freeze_iterativeLSI_seurat.obj.rds"))
seurat.obj.big$LSI_CellType <- Idents(seurat.obj.big)
table(seurat.obj.big$orig.ident)

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("/labs/ccurtis2/mjprzy/scRNA_analysis/LSI_projection/merged_data_ZhangEtal_freeze/SummarizedExperiment_merged_data_ZhangEtal_freeze_UMAP-model.uwot")

######################################################################################
# PLOT UMAP EMBEDDING OF THE REFERENCE/BACKGROUND DATASET ONLY
######################################################################################

# subset to background onlys
pdf(paste0(plotDir, sample.tmp,"_coloured_background_only.pdf"), width = 12, height = 12, useDingbats = FALSE)
DimPlot(seurat.obj.big, reduction = "LSI_umap", pt.size = 1.5) +
  theme_bw() +
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") +
  scale_color_manual(values=c("Chief Cells" = "lightgrey", "Neck-like Cells" = "lightgrey", "Enteroendocrine Cells" = "lightgrey", "PMCs" = magma(9)[5], 
                              "Goblet Cells" = "lightgrey", "MSCs" =  magma(9)[4], "Enterocytes" = "lightgrey",
                              "GMCs" = "lightgrey", "Endocrine Cells" = "lightgrey", "Proliferative Cells" = "lightgrey", "Malignant Cells" = magma(9)[3],
                              "Paneth Cells" = "lightgrey")) +
  theme_classic() + NoLegend() +
  theme(axis.title = element_text(size = 40, face = "bold")) +
  theme(legend.position="bottom", legend.title=element_text(size=25, face = "bold")) + 
  theme(legend.text = element_text(colour="black", size=25, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
  theme(axis.text = element_text(face="bold", color="black", size=20)) +
  theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
  guides(color=guide_legend("Patients", override.aes = list(size = 6)))
dev.off()

###############################################################################
# ITERATIVE OVER ALL SAMPLES IN THE EL ORGANOID DATASET
###############################################################################

# set up nearest neighbor dataframe for all samples
nn.table <- table(seurat.obj.big$LSI_CellType) %>% bind_rows
nn.table <- as.data.frame(nn.table)

# setup a projection plot directory
plotDir <- paste0(sample.tmp, "/classification/")
dir.create(plotDir,recursive=TRUE)

# iterative over each individual sample
for (i in 1:length(sample.ids)){
  
  # set id from samples
  id <- sample.ids[i]
  
  # set id from patients
  # id <- patient.ids[i]
  print(id)
  
  # get only the cells associated with this sample
  seDisease <- sce[,grep(id, sce@metadata$HashTag)]
  
  # get only the cells associated with this patient
  # seDisease <- sce[,grep(id, sce@metadata$sample_type)]
  
  #Identify Gene Universe
  gU <- intersect(rownames(seDisease), rownames(se))
  gU <- gU[!grepl("^MT", gU)]
  
  ######################################################################################
  # CLUSTERING ANALYSIS STARTS HERE
  ######################################################################################
  # Set Clustering Parameters
  varGenesToUse <- c(1000, 1000, 1000) #Choose a higher number of variable peaks
  resolution <- c(0.2,0.8,0.8) #Clustering resolutions for Seurat SNN
  verbose = TRUE
  tstart <- NULL
  input_knn <- 25
  scaleTo <- 10000
  nMax <- 3000
  
  # set clustering parameters for seurats find cluster method
  clustParams <- list()
  clustParams$verbose <- verbose
  clustParams$tstart <- tstart
  
  #Optimize LSI Features for STP samples
  matAll <- cbind(assay(seDisease[gU,]), assay(se[gU,])) # combine single sample with the Zhang et al dataset
  lsiObj <- optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse,  pcsUse = 1:25)
  
  #UMAP
  set.seed(1)
  umap <- uwot::umap(
    lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25], 
    n_neighbors = 35, 
    min_dist = 0.5, 
    metric = "euclidean", 
    n_threads = 5, 
    verbose = TRUE, 
    ret_model = FALSE
  )
  
  #Plot Info
  cells <- c(rep("reference", ncol(se)),rep("disease", ncol(seDisease)))
  splitCells <- split(cells,lsiObj[[length(lsiObj)]]$clusters)
  df <- data.frame(
    clusters = names(splitCells),
    proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="disease") / length(splitCells[[x]])))
  )
  
  # get the celltype assignment here + barcode
  celltype.df <- data.frame(Cell_barcode = rownames(seurat.obj.big@meta.data), CellType = seurat.obj.big$LSI_CellType)
  
  #Plot UMAP Data Frame
  plotDF <- data.frame(umap)
  rownames(plotDF) <- c(colnames(se), colnames(seDisease))
  plotDF$Cell_barcode <- c(colnames(se), colnames(seDisease))
  plotDF$type <- cells
  plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
  
  # merge with the celltype information
  plotDF <- merge(plotDF, celltype.df, by = "Cell_barcode", all.x = T)
  plotDF$CellType <- as.character(plotDF$CellType)
  plotDF[is.na(plotDF$CellType), "CellType"] <- "organoid"
  plotDF$classification <- 0
  
  #If disease cells are clustered with healthy cluster (proportion > 0.8) we will classify these as healthy-like
  plotDF$classification[plotDF$CellType == "Malignant Cells" & plotDF$clusters %in% paste0(df$clusters[df[,2] > 0.05])] <- 1
  plotDF$classification[plotDF$type == "disease"] <- plotDF$classification[plotDF$type == "disease"] + 1
  plotDF <- plotDF[order(plotDF$classification), ]
  
  #Formal Classification
  plotDF$classificationSTR <- "reference"
  plotDF$classificationSTR[plotDF$classification==1] <- "healthy-like"
  plotDF$classificationSTR[plotDF$classification==2] <- "disease-like"
  
  #Plot PDFs
  pdf(paste0(plotDir,id,"-Classification-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
  print(ggplot(plotDF, aes(X1,X2,color=classificationSTR)) + 
          geom_point() +
          theme_bw() +
          xlab("UMAP Dimension 1") + 
          ylab("UMAP Dimension 2") +
          scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3")) +
          theme_classic() +
          theme(axis.line=element_blank(),axis.text.x=element_blank(),
                axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),legend.position="none",
                panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),plot.background=element_blank()))
  dev.off()
  
  ####################################################
  # PERFORM PROJECTION INTO LSI UMAP
  ####################################################
  
  ## in the next step, I guess it is assumed, that the matrices contain mainly overlapping genes
  ## however, this does not seem to be the case, which is why I adapted the objects in the following
  #LSI Projection Matrix
  lsiGenes <- metadata(se)$variableGenes
  matProjectLSI <- assay(seDisease[lsiGenes,])
  
  #LSI Project
  lsiReference <- metadata(se)$optimizeLSI[[length(metadata(se)$optimizeLSI)]]$lsiObj
  
  # run projectionn
  lsiProjection <- projectLSI(matProjectLSI, lsiReference)
  
  #UMAP Projection
  #Set Seed Prior to umap_transform (see uwot github)
  set.seed(1)
  umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:25], umapManifold, verbose = TRUE)
  
  ######################################################################################
  # PLOT PROJECTION ACCORDING TO SAMPLES
  ######################################################################################
  
  #Plot Projection
  refDF <- data.frame(row.names = colnames(se), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
  proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = "organoid")
  projectionDF <- rbind(refDF, proDF)
  
  # read in the cnv clones tables
  cnv.clones.list <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric/arm_level_matrices", pattern = "cnv_clones", full.names = T)
  cnv.clones.list <- cnv.clones.list[grep(paste0(str_split_fixed(id, "_", 3)[,1], "_",str_split_fixed(id, "_", 3)[,2]), cnv.clones.list)]
  
  # read tables in
  cnv.clones.list <- lapply(cnv.clones.list, read.table, header = T)
  cnv.clones.df <- bind_rows(cnv.clones.list)
  
  # this is a confusionn matrix with the timepoint against the genotype
  cnv.tmp <- table(cnv.clones.df$cnv_id)
  cnv.clones <- cnv.tmp[cnv.tmp>10]
  
  if (nrow(cnv.clones.df) > 0){
    # add 4230_EL to it
    cnv.clones.df$Cell_barcodes <- paste0(str_split_fixed(id, "_", 2)[,1], "_EL_", cnv.clones.df$Cell_barcodes)
  }
  
  # get the celltype assignment here + barcode
  celltype.df <- data.frame(Cell_barcodes = rownames(seurat.obj.big@meta.data), CellType = seurat.obj.big$LSI_CellType, clusters = se$Clusters)
  
  # create UMAP dataframe merge with the celltype assignment from the background
  projectionDF$Cell_barcodes <- rownames(projectionDF)
  projectionDF <- merge(projectionDF, celltype.df, by = "Cell_barcodes", all.x = T)
  projectionDF$CellType <- as.character(projectionDF$CellType)
  
  # assign the organoid cell type to the organoids which do not have a celltype assignment
  projectionDF[is.na(projectionDF$CellType), "CellType"] <- "organoid"
  rownames(projectionDF) <- projectionDF$Cell_barcodes
  
  # order the dataframe to project organoids on reference and not the other way round
  projectionDF$Type <- factor(projectionDF$Type, levels = c("reference", "organoid"))
  projectionDF <- projectionDF[order(projectionDF$Type),]
  
  ######################################################################################
  # PLOT PROJECTION ACCORDING TO SAMPLES WITH COLOURING ACCORDING TO EARLY LATE WT
  ######################################################################################
  
  
  if (!is.integer0(grep("EARLY", id))){
    
    ## EARLY
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "#EA9F37")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="#EA9F37",n=300,h=3))
    dev.off()
    
  } else if (!is.integer0(grep("LATE", id))){
    
    ## LATE
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "#782867")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="#782867",n=300,h=3))
    dev.off()
    
  } else{
    
    ## WT
    pdf(paste0(plotDir,id,"-Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(projectionDF, aes(X1,X2,color=Type)) + 
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") + 
            ylab("UMAP Dimension 2") +
            scale_color_manual(values=c("reference"="lightgrey", "organoid" = "#4F9E4C")) +
            theme_classic() +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),legend.position="none",
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank()) +
            geom_density_2d(data=projectionDF[projectionDF$Type == "organoid",], mapping=aes(x=X1,y=X2),color="#4F9E4C",n=300,h=3))
    dev.off()
  }
  
  ######################################################################################
  # PLOT PROJECTION ACCORDING TO SAMPLES COLOURED ACCORDING TO COPY NUMBER CLONES
  ######################################################################################
  
  if (nrow(cnv.clones.df) > 0){
    # add copy number clones to the projectionDF
    cnv.projectionDF <- merge(projectionDF, cnv.clones.df, by = "Cell_barcodes", all.x = T)
    cnv.projectionDF[is.na(cnv.projectionDF$clone_id) & cnv.projectionDF$Type == "organoid", "clone_id"] <- "unassigned"
    cnv.projectionDF$clone_id <- as.character(cnv.projectionDF$clone_id)
    
    # replace the NAs and sort the dataframe
    cnv.projectionDF[is.na(cnv.projectionDF$clone_id), "clone_id"] <- "background"
    cnv.projectionDF <- cnv.projectionDF[order(cnv.projectionDF$Type),]
    
    # add new levels for sorting
    cnv.projectionDF$clone_id <- factor(cnv.projectionDF$clone_id, levels = c("background", paste0("clone", 1:200)))
    
    # make a custom color code for the cnv clones
    colourCount1 = length(unique(cnv.projectionDF$clone_id))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    
    pdf(paste0(plotDir,id,"-Projection-UMAP-CNV_clones.pdf"), width = 12, height = 12, useDingbats = FALSE)
    print(ggplot(cnv.projectionDF, aes(X1,X2,color=clone_id)) +
            geom_point() +
            theme_bw() +
            xlab("UMAP Dimension 1") +
            ylab("UMAP Dimension 2") +
            scale_color_manual(values = c("background" = "lightgrey", getPalette(colourCount1-1))) +
            theme_classic() +
            theme(axis.title = element_text(size = 40, face = "bold")) +
            theme(legend.position="bottom", legend.title=element_text(size=15, face = "bold")) +
            theme(legend.text = element_text(colour="black", size=10, face="bold"),
                  legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
            theme(axis.text = element_text(face="bold", color="black", size=20)) +
            theme(axis.line = element_line(colour = "black", size = 1.5, linetype = "solid")) +
            guides(color=guide_legend("CNV Clones", override.aes = list(size = 6))))
    dev.off()
  }
  
  #############################################################
  #PERFORM DIFFERENTIAL GENE EXPRESSION ANALYSIS IN LSI SPACE
  #############################################################
  
  #LSI-SVD
  svdReference <- as.data.frame(lsiReference$matSVD)
  svdDisease <- as.data.frame(as.matrix(lsiProjection))
  
  #Differential Seed
  set.seed(1)
  
  #Cells that we are testing of disease
  idxDisease <- rownames(projectionDF)[projectionDF$CellType=="organoid"]
  
  #If the number of cells exceeds the max downsample to max
  if(length(idxDisease) > nMax){
    idxDisease <- sample(idxDisease, nMax)
  }
  
  #If the number of cells is greater than 1 continue
  stopifnot(length(idxDisease) > 1)
  
  #KNN Nearest Neighbor using FNN
  knnDisease <- FNN::get.knnx(
    data = svdReference,
    query = svdDisease[idxDisease, ], #Subset by idxDisease 
    k = input_knn)
  
  # nn.index - an n x k matrix for the nearest neighbor indice.
  # make a nearest neighbor dataframe
  nn.df <- data.frame(nn_index = as.vector(knnDisease$nn.index), Cell_barcodes = str_split_fixed(rownames(svdReference[as.vector(knnDisease$nn.index),]), "-", 2)[,1])
  nn.df <- nn.df[order(nn.df$nn_index),]
  nn.df$Cell_barcodes <- paste0(nn.df$Cell_barcode, "-1")
  
  # merge the nn.df with the projectionDF from before to get information about cell identity etc.
  full.nn.df <- merge(nn.df, projectionDF, by = "Cell_barcodes", all.x = T)
  
  # nn overview of cells for the sample
  nn.frequency <- table(full.nn.df$CellType) %>% bind_rows
  nn.table <- rBind.fill(nn.table, as.data.frame(nn.frequency))
  rownames(nn.table)[i+1] <- id
  
}

# write the KNN table to file
write.table(nn.table, paste0(sample.tmp, "/", sample.tmp, "_KNN_celltype_table.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

################################################################
#       CREATE A PLOT FOR THE DRIFT QUANTIFICATION
################################################################
# nn.table <- read.delim(paste0(sample.tmp, "/", sample.tmp, "_KNN_celltype_table.txt"), header = T)
colnames(nn.table) <- c("Chief Cells", "PMCs", "Enterocytes", "Enteroendocrine Cells", "Parietal Cells", "Goblet Cells", "Mucosal-like malignant Cells", 
                        "GMCs", "Non-Mucosal-like malignant Cells", "Proliferative Cells")
# calculate percentages of nearest neighbors 
# then compare for each sample - early late
# is the frequency of cells mapping to cancer increasing
nn.table <- nn.table[-1,]

# remove samples which do have less than 50 cells in total
n_cells_clones <- table(metadata(sce)$HashTag) > 25
n_cells_clones <- n_cells_clones[n_cells_clones == T]
nn.table <- nn.table[rownames(nn.table) %in% names(n_cells_clones),]

# remove 6077 C1
nn.table <- nn.table[-grep("6077_C1_EARLY", rownames(nn.table)),]
nn.table <- nn.table[-grep("0891_PB2_LATE", rownames(nn.table)),]

# calculate frequencies instead of number of nearest neighbors
sample.cells <- rowSums(nn.table)
freq.table <- round(nn.table/sample.cells, 3)*100

# change 6077 WT
freq.table <- rbind(freq.table, freq.table[5,],freq.table[5,])
rownames(freq.table)[5] <- "6077_C1_WT"
rownames(freq.table)[18] <- "6077_C5_WT"
rownames(freq.table)[19] <- "6077_D3_WT"

# change 0891 WT
freq.table <- rbind(freq.table, freq.table[17,],freq.table[17,])
rownames(freq.table)[17] <- "0891_PB2_WT"
rownames(freq.table)[20] <- "0891_P31_WT"
rownames(freq.table)[21] <- "0891_P41_WT"

# change 4230 WT
freq.table <- rbind(freq.table, freq.table[12,],freq.table[12,])
rownames(freq.table)[12] <- "4230_APA6_WT"
rownames(freq.table)[22] <- "4230_APA4_WT"
rownames(freq.table)[23] <- "4230_PB6_WT"

# add sample information
freq.table$sample_id <- as.character(rownames(freq.table))

# melt the freq table to plot it 
melt.freq.table <- melt(as.data.frame(freq.table), id.var = "sample_id")

# add timepoint information
melt.freq.table$timepoint <- "WT"
melt.freq.table[grep("EARLY", melt.freq.table$sample_id),"timepoint"] <- "EARLY"
melt.freq.table[grep("LATE", melt.freq.table$sample_id),"timepoint"] <- "LATE"
melt.freq.table$timepoint <- factor(melt.freq.table$timepoint, levels = c("WT", "EARLY", "LATE"))

# add clone information
melt.freq.table$clone_id <- paste0(str_split_fixed(melt.freq.table$sample_id, "_", 3)[,1], "_", str_split_fixed(melt.freq.table$sample_id, "_", 3)[,2])

# replace the cell names that are not correctly shown
melt.freq.table$variable <- as.character(melt.freq.table$variable)

# rename some cell types
melt.freq.table[grep("PMCs", melt.freq.table$variable),"variable"] <- "Pit Mucosal Cells"
melt.freq.table[grep("GMCs", melt.freq.table$variable),"variable"] <- "Gland Mucosal Cells"
melt.freq.table[grep("MSCs", melt.freq.table$variable),"variable"] <- "Mucosal Stem Cells"

# order the cell types
melt.freq.table$variable <- factor(melt.freq.table$variable, levels = c("Malignant Cells", "Mucosal Stem Cells",
                                                                        "Pit Mucosal Cells", "Neck-like Cells",
                                                                        "Chief Cells","Gland Mucosal Cells",
                                                                        "Enteroendocrine Cells", "Endocrine Cells", "Enterocytes", "Goblet Cells"))
# # remove cells with maximum contribution of 2.5%
# melt.freq.table <- melt.freq.table[melt.freq.table$variable %in% c("Mucosal-like malignant Cells", "Non-Mucosal-like malignant Cells",
#                                                                    "Pit Mucosal Cells", "Mucosal Stem Cells",
#                                                                    "Chief Cells","Gland Mucosal Cells",
#                                                                    "Enteroendocrine Cells", "Endocrine Cells",
#                                                                    "Enterocytes", "Proliferative Cells"), ]

# plot data
ggplot(data = melt.freq.table, aes(x = timepoint, y = value, group = clone_id)) +
  geom_line(aes(color = clone_id, alpha = 1), size = 1) +
  geom_point(aes(color = clone_id, alpha = 1), size = 2) +
  facet_wrap(~ variable, ncol = 2) +
  theme_bw() + 
  guides(alpha = FALSE) +
  labs(title = "K-Nearest Neighbor Cell Type Quantification", x = "Time", y = "Cell type frequency [%]", color = "Clone IDs") +
  theme(strip.text = element_text(face="bold", size=14, colour = "black",),
        strip.background = element_rect(fill="lightgrey", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 14, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold" ),
        axis.title = element_text(colour = "black", size = 20, face = "bold" ),
        plot.title = element_text(colour = "black", size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 14, face = "bold",),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.position = "bottom") +
  geom_text_repel(data = melt.freq.table %>% filter(timepoint == "LATE") %>% filter(value > 25),
                  aes(label = paste0(clone_id, " - ", value, "%")) ,
                  hjust = -.35,
                  nudge_x = .5,
                  direction = "y",
                  fontface = "bold",
                  size = 3) +
  scale_color_viridis_d(option="D") +
  guides(color=guide_legend( override.aes = list(size = 3)))
ggsave(paste0(sample.tmp, "/AllPatients_LSI_quantification_", sample.tmp ,".pdf"), width = 7.5, height = 16.5)

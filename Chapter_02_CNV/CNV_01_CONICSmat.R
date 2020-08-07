################################################################################################################################################
##                                                                                                                      
##  RUN CONICSmat FOR SEQUENCING 8 - P2C2R2T2                                                                                           
##                                                                                                                      
##  Date: 11 March 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script uses the Seurat object from Sequencing 8, corresponding to Sample P2C2R2 at time point 2.
##           
##                                                                                                                      
################################################################################################################################################
# clear workspace
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("beanplot", "mixtools", "pheatmap", "tidyverse", "zoo", "squash", "biomaRt",
                      "CONICSmat", "Rtsne", "scran", "stringr", "Matrix.utils", "Seurat", "reshape2", "readr", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                          DEFINE FUNCTIONS
############################################################################

# get count matrix with row and colnames
getCountMat <- function(cellranger_outs_folder){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  features.path <- file.path(matrix_dir, "features.tsv.gz")
  matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")c
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  return( mat )
}

# get the barcodes for the matrices
getBarcodeNames <- function(cellranger_outs_folder){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  return( as.character(barcode.names[,1]) )
}

# merge the matrices together
check_and_merge <- function(sparse_matrix_list){
  if(length(sparse_matrix_list) < 2){
    return(sparse_matrix_list[[1]])
  }
  combtab <- combn(x = seq(length(sparse_matrix_list)),m = 2)
  
  genes_in_common <- list()
  for(i in seq(ncol(combtab))){
    genes_in_common[[i]] <- intersect(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]]))
  }
  
  cg <- genes_in_common[[1]]
  if(length(genes_in_common) > 1 ){
    for(i in 2:length(genes_in_common)){
      cg <- intersect(cg,genes_in_common[[i]])
    }
  }
  
  filter_genes <- function(x, genes){
    x <- x[genes,]
    return(x)
  }
  
  sparse_matrix_list <- lapply(sparse_matrix_list, filter_genes, cg)
  
  check.cols <- c()
  for(i in seq(ncol(combtab))){
    check.cols <- c(check.cols, identical(rownames(sparse_matrix_list[[combtab[1,i]]]), rownames(sparse_matrix_list[[combtab[2,i]]])) )
  }
  if(all(check.cols)){
    outmat <- do.call(cbind, sparse_matrix_list)
    return(outmat)
  } else {
    return(NULL)
  }
}

`%notin%` <- Negate(`%in%`)

############################################################################
##                    SET DIRECTORIES AND IDS
############################################################################
# define working directory
w.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis"

# get samples
sample.list <- list.files(w.dir, pattern = "Sequencing8", recursive = T, full.names = T, include.dirs = T)

# get seurat objects
seurat.list <- sample.list[grep("cr3_seurat_obj.rds", sample.list)]

# get output directory
o.dir <- "/labs/ccurtis2/mjprzy/correlation_test"
setwd(o.dir)

# get sample.ids
sample.ids <- c("P2C2R2T2")

# set sample name
sample.tmp <- sample.ids[1]

# create an output directory
dir.create(sample.tmp)
setwd(sample.tmp)

message(paste0("Start time: ", Sys.time()))

############################################################################
##                    READ IN DATA OF INTEREST
############################################################################

## the input for CONICSmat is a genes X cells matrix of log2(CPM/10+1) values
# read in object for the sample of interest
seurat.obj <- readRDS(seurat.list[1])

# read in metadata
seurat.metadata <- seurat.obj@meta.data

# grep the positive and negative cells
pos.seurat.metadata <- seurat.metadata

# get seurat object for the observation matrix
mat.pos <- seurat.obj@assays$RNA@counts

# get the positive cells
names.use.pos <- colnames(mat.pos)[colnames(mat.pos) %in% pos.seurat.metadata$Cell_Barcode]

# get seurat obj for the reference sample
WT.seurat.obj <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/Sequencing13_6077_EL_cr3/Sequencing13_6077_EL_cr3_seurat_obj.rds")
neg.seurat.metadata <- WT.seurat.obj@meta.data
neg.seurat.metadata <- neg.seurat.metadata[which(neg.seurat.metadata$Sample_Origin == "6077_WT"),]

# for WT.seurat.obj
mat.neg <- WT.seurat.obj@assays$RNA@counts

# get the negative cells
names.use.neg <- colnames(mat.neg)[colnames(mat.neg) %in% neg.seurat.metadata$Cell_Barcode]

# make matrix lists
tumor.matrix <- mat.pos[, names.use.pos]
tumor.matrix <- normMat(as.matrix(tumor.matrix))
WT.matrix <- mat.neg[, names.use.neg]
WT.matrix <- normMat(as.matrix(WT.matrix))
colnames(WT.matrix) <- paste0(colnames(WT.matrix), "_WT")

# merge matrices
all.matrix <- check_and_merge(sparse_matrix_list = list(tumor.matrix, WT.matrix))

# check dims
dim(tumor.matrix)
dim(all.matrix)
tail(colnames(all.matrix))

# define a vector with distinct labels - here tumor and WT
cell.labels <- unlist(str_split_fixed(colnames(all.matrix),"_",2)[,2])
cell.labels[1:ncol(tumor.matrix)] <- "tumor"

# read in chromsome arm positions file
regions <- read.table("/home/mjprzy/correlation_test/chr_arm_position_grch38.txt", sep = "\t", row.names = 1, header = T)

############################################################################
##                      CONICSmat PROCESSING
############################################################################

# get chromosomal positions of genes in the matrix
gene_pos <- getGenePositions(rownames(all.matrix))

# filter genes which are present in less than 5 cells
filtered.all.matrix <- filterMatrix(all.matrix, gene_pos[,"hgnc_symbol"], minCells=5)

# calculate the normalization factor for each cell
# normalizatinon factor centers the gene expression in each gene around the mean
normFactor <- calcNormFactors(filtered.all.matrix)

# next step is to determine if the average gene expression any of the regions show a bimodal distribution across cells
# First, each cell is centered using the previously calculated normalization factor. 
# Then, the z-score of the centered gene expression across all cells is calculated. 
# Based on these z-scores, a Gaussian mixture model is calculated with the mixtools package
# only calulcated results for regions with more than 100 expressed genes
l <- plotAll(filtered.all.matrix, normFactor, regions, gene_pos, paste0(sample.tmp,"_GO_CNVs"), repetitions = 10, postProb = 0.75)

# just chr7
r=plotChrEnichment(filtered.all.matrix, 7, normFactor, gene_pos, repetitions = 200, postProb = 0.75)

# results are visualized as apdf file, with the a matrix including Bayesian information criterion (BIC)
# and the adjusted likelihood ratio test (LRT) is written
# for BIC, the model with the lowest BIC is preferred

# visualize heatmap of the posterior probabilites of cells
# posterior probabilites are z-scored
# visualization of component2, which is the component with the larger mean - 
# cells with higher expression at that locus will appear in red, lower expression in blue
hi <- plotHistogram(l, filtered.all.matrix, clusters = 1, zscoreThreshold=4, patients = cell.labels)
pdf(paste0(sample.tmp, "_prefilter_histogram.pdf"))
plotHistogram(l, filtered.all.matrix, clusters = 1, zscoreThreshold=4, patients = cell.labels)
dev.off()

# To obtain a final assignmant as malignant or non-malignant cells, we first filter uninformative, noisy regions 
# based on the results of the likelihood ratio test and the BIC for each region.
lrbic <- read.table(paste0(sample.tmp, "_GO_CNVs_BIC_LR.txt"), sep="\t", header=T, row.names=1, check.names=F)
colnames(lrbic)

# get candidate regions
candRegions <- rownames(lrbic)[which(lrbic[,"BIC difference"]>30 & lrbic[,"LRT adj. p-val"] < 0.01)]

# gernate a histogram for the most informative regions again
hi <- plotHistogram(l[,candRegions], filtered.all.matrix, clusters=1, zscoreThreshold=4, patients = cell.labels)
pdf(paste0(sample.tmp, "_postfilter_histogram.pdf"))
plotHistogram(l[,candRegions], filtered.all.matrix, clusters=1, zscoreThreshold=4, patients = cell.labels)
dev.off()

# we can now assign a label as malignant or non-malignant to each cell based on the clusters returned by the heatmap function.
WT.cell.index <- grep("WT", names(hi))
names(WT.cell.index) <- names(hi[WT.cell.index])
tumor.cell.index <- c(1:ncol(tumor.matrix))
names(tumor.cell.index) <- names(hi[tumor.cell.index])
hi[WT.cell.index] <- 1
hi[tumor.cell.index] <- 2
normal <- WT.cell.index
tumor <- tumor.cell.index

############################################################################
##              FILTER CHROMOSOME ARMS AND BINARIZE ASSIGNMENT
############################################################################
# now we can plot the posterior probabilities again, but with statistics for normal & tumor cells
redu <- plotAll(filtered.all.matrix, normFactor, regions[candRegions,], gene_pos, paste0(sample.tmp, "_GO_CNVs_with_info.pdf"), repetitions = 10, postProb = 0.75, normal=normal,tumor=tumor)

# just chr3
r <- plotChrEnichment(filtered.all.matrix, 3, normFactor, gene_pos, groups1 = normal, groups2 = tumor, repetitions = 200, postProb = 0.75)

# thresholding on the posterior probabilites, we can generate a binnary matrix 
# 1 indicates the presence of a CNV & 0 the absence
bin_mat <- binarizeMatrix(redu, normal, tumor, 0.8)
write_rds(bin_mat, "CONICSmat.binMatrix.rds")

# plot the whole binary matrix
pdf(paste0(sample.tmp, "_binary_matrix.pdf"))
plotBinaryMat(bin_mat, cell.labels, normal, tumor)
dev.off()

# subset matrix
bin_mat <- bin_mat[complete.cases(bin_mat),]
dim(bin_mat)

# summarize the resulting matrix to order cells by "mutations"
geno <- sort(table(apply(bin_mat, 1, paste, collapse=",")), decreasing = T)
geno
geno_sort=geno[sort(names(geno[which(geno>15)]))]
length(geno)

############################################################################
##                    BREAKPOINT VISUALIZATION
############################################################################
# visualize the breakpoints for each chrosomome
# just check at the most important regions, not all segments
# uses pooling of cells and somooting the expressin of genes in a sliding window of 101 windows
# prior to plotting, all ratios are median-centered, assuming that less than 50% of the genome are altered by CNVs
plotAllChromosomes (mat = filtered.all.matrix, normal = normal, tumor = tumor, windowsize = 101, gene_pos = gene_pos, fname = paste0(sample.tmp, "_chr_breakpoints"), patients = cell.labels, breakpoints = regions)

# generate a matrix of corrected p-values for each cell and each CNV by using a cell population without 
# CNVs for each region. Use parameter estimation (mean and sdv) to calculate a p-value (adjusted with Benjamini-Hochberg)
# the matrix r holds adjusted p-values for each cell and CNV candidate region
r <- generatePvalMat(filtered.all.matrix, regions[candRegions,], normFactor, normal, tumor, gene_pos, threshold=0.8)

# binr is a binarized version with a p-val < 0.1
binr <- ifelse(r > 0.1, 0,1)

# plot distribution of  adjusted p-values for CNV presence in tumor cells in selected regions
pdf(paste0(sample.tmp, "_boxplot.pdf"))
boxplot(r)
dev.off()

############################################################################
##                PLOT HEATMAP AND SAVE MATRIX TO FILES
############################################################################
# lastly, we can generate a visualization of chromosomal alteration in each single cell
# calculate the average expression in normal cells in a smoothing window (sorted by genomic position)
# center the expresion ratio to a normal by the mean
# windowsize determines the size of the smoothing window, expThresh determines the minimum average expression in control and tumor cells
final.matrix <- plotChromosomeHeatmap(filtered.all.matrix, normal = normal, 
                                      plotcells = 1:length(cell.labels), gene_pos = gene_pos, 
                                      windowsize = 151, chr=T, expThresh=0.2, thresh = 10, 
                                      colo=c(rep("blue",length(normal)),rep("red",length(tumor))), retMat = T, plotdendrogram = T)

# plot heatmap
png(paste(sample.tmp,"_ExpHeatmap.png",sep=""), width = 1000,height = 600)
plotChromosomeHeatmap(filtered.all.matrix[,c(normal,tumor)],normal = 1:length(normal),plotcells = 1:(length(normal)+length(tumor)),chr=T, gene_pos = gene_pos,windowsize = 151,thresh = 0.6, expThresh = 0.2,colo = c(rep("blue",length(normal)),rep("red",length(tumor))))
abline(h=length(normal),lty=16)
dev.off()

# transpose the output matrix to have genes x cells
t.final.matrix <- t(final.matrix)
t.final.matrix <- t.final.matrix[,-1]

# remove normal cells from the matrix
tumor.final.matrix <- t.final.matrix[,colnames(t.final.matrix) %in% colnames(tumor.matrix)]
colnames(tumor.final.matrix) <- str_split_fixed(colnames(tumor.final.matrix), "_", 2)[,1]

# save matrix
write_rds(tumor.final.matrix, "CONICSmat.matrix.rds")

# message the end time
message(paste0("End time: ", Sys.time()))

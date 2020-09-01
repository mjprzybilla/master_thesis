#############################################################################################################################
##                                                                                                                      
##  RUN HONEYBADGER ON SEQUENCING 8 CORRESPONDING TO P2C2R2 AT TIMEPOINT 2
##                                                                                                                      
##  Date: 06 August 2020                                                                                                                    
##  
##  Author: Moritz Przybilla
##                                                                                                                      
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "tidyverse", "Matrix", "Seurat", "GenomeInfoDb", "BSgenome", "devtools", "HoneyBADGER", "biomaRt")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# set working direcotry 
setwd("/labs/ccurtis2/mjprzy/correlation_test")
message(paste0("Start time: ", Sys.time()))

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
##                  SET UP THE DATA MATRICES OF INTEREST
############################################################################

# get the mart object
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl')

# define working directory
w.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis"

# get samples
sample.list <- list.files(w.dir, pattern = "Sequencing8", recursive = T, full.names = T, include.dirs = T)

# get seurat objects
seurat.list <- sample.list[grep("cr3_seurat_obj.rds", sample.list)]

# get seurat obj for the respective sample
seurat.obj <- readRDS(seurat.list[1])
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
mat.pos.list <- mat.pos[, names.use.pos]
mat.neg.list <- mat.neg[, names.use.neg]

# add labels for observatoin and ref
colnames(mat.pos.list) <- paste0(colnames(mat.pos.list), "_pos")
colnames(mat.neg.list) <- paste0(colnames(mat.neg.list), "_neg")

# combine both matrices together
all_mat <- check_and_merge(sparse_matrix_list = list(mat.pos.list, mat.neg.list))
normal_mat <- mat.neg.list

############################################################################
##                  RUN THE FIRST STAGE OF HONEYBADGER
############################################################################

hb_norm <- new('HoneyBADGER', name='hGasHB')
hb_norm$setGexpMats(all_mat, normal_mat, mart.obj, filter=TRUE, scale=TRUE, verbose=TRUE,
                    minMeanBoth =0.1, minMeanTest =0.1, minMeanRef =0.1) #5k genes

png("ECB_min10_ok.png",width = 24, height = 8, units = 'in', res = 600)
honeySmooth <- hb_norm$plotGexpProfile(zlim =c(-0.6,0.6),setOrder=FALSE,cellOrder=all_cells,returnPlot=TRUE)
dev.off()

message(paste0("End time: ", Sys.time()))

# write.table(honeySmooth[[11]],"HoneySmooth_chr13.txt",sep="\t",quote = FALSE)

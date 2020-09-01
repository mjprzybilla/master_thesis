################################################################################################################################################
##                                                                                                                      
##  Prepare input for inferCNV for Sequencing 19                                                                                        
##                                                                                                                      
##  Date: 11 March 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script shall produce the inferCNV inputs to use with inferCNV_analysis.R afterwards.   
##           
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("Matrix", "tidyverse", "biomaRt", "httr", "Seurat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
set_config(config(ssl_verifypeer = 0L))

############################################################################
##                              Define Functions
############################################################################

# get count matrix with row and colnames
getCountMat <- function(cellranger_outs_folder){
  matrix_dir = file.path(cellranger_outs_folder,"filtered_feature_bc_matrix")
  barcode.path <- file.path(matrix_dir, "barcodes.tsv.gz")
  features.path <- file.path(matrix_dir, "features.tsv.gz")
  matrix.path <- file.path(matrix_dir, "matrix.mtx.gz")
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
##                              Read in data
############################################################################

# get input directory for tumor samples
tumor.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze"

# list files
tumor.dir <- list.files(tumor.dir, full.names = T, pattern = "Sequencing19")

# get output directory
o.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric/freeze"

# create output directory
dir.create(o.dir)
setwd(o.dir)

# get sample names
sample.ids <- basename(tumor.dir)

# make list of working directories
w.dir_list <- file.path(o.dir, paste0(sample.ids, "_WT"))

# w.dir <- w.dir_list[1]
for(w.dir in w.dir_list){
  
  # create sample working/output directory
  dir.create(w.dir)
  message(w.dir)
  
  # setwd 
  setwd(w.dir)
  
  ############################################################################
  ##   Get count matrices from which tumor barcodes will be extracted
  ############################################################################
  
  ## THIS IS THE OBJECT OF INTEREST
  # get seurat obj for the respective sample
  seurat.obj <- readRDS(paste0(tumor.dir,"/", sample.ids, "_seurat_obj.rds"))
  seurat.metadata <- seurat.obj@meta.data
  
  # rename ids 
  seurat.metadata[seurat.metadata == "LATE_4230_APA6"] <- "4230_APA6_LATE"
  seurat.metadata[seurat.metadata == "EARLY_4230_PB6"] <- "4230_PB6_EARLY"
  seurat.metadata[seurat.metadata == "LATE_4230_APA4"] <- "4230_APA4_LATE"
  seurat.metadata[seurat.metadata == "EARLY_4230_APA6"] <- "4230_APA6_EARLY"
  seurat.metadata[seurat.metadata == "EARLY_4230_APA4"] <- "4230_APA4_EARLY"
  seurat.metadata[seurat.metadata == "LATE_4230_PB6"] <- "4230_PB6_LATE"
  seurat.metadata[seurat.metadata == "4230_WT_Passage_midi"] <- "4230_WT"
  
  # grep the positive cells
  pos.seurat.metadata <- seurat.metadata[which(seurat.metadata$Sample_Origin != "4230_WT"),]
  
  # get seurat object for the observation matrix
  mat.pos <- seurat.obj@assays$RNA@counts
  
  # get the positive cells
  names.use.pos <- colnames(mat.pos)[colnames(mat.pos) %in% pos.seurat.metadata$Cell_Barcode]
  
  # get seurat obj for the reference sample
  WT.seurat.obj <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/Sequencing5_cr3/Sequencing5_cr3_seurat_obj.rds")
  neg.seurat.metadata <- WT.seurat.obj@meta.data
  neg.seurat.metadata <- neg.seurat.metadata[which(neg.seurat.metadata$Sample_Origin == "4230_WT"),]
  
  # for WT.seurat.obj
  mat.neg <- WT.seurat.obj@assays$RNA@counts
  
  # get the negative cells
  names.use.neg <- colnames(mat.neg)[colnames(mat.neg) %in% neg.seurat.metadata$Cell_Barcode]
  
  # make matrix lists
  mat.pos.list <- list(mat.pos[, names.use.pos])
  mat.neg.list <- list(mat.neg[, names.use.neg])
  
  # add label for cells from the same sample
  for(i in seq(length(mat.pos.list))){
    colnames(mat.pos.list[[i]]) <- paste0(colnames(mat.pos.list[[i]]), "_pos")
  }
  
  for(i in seq(length(mat.neg.list))){
    colnames(mat.neg.list[[i]]) <- paste0(colnames(mat.neg.list[[i]]), "_neg")
  }
  
  ############################################################################
  ##
  ############################################################################
  
  # merge the matrix lists together
  mat.pos <- check_and_merge(sparse_matrix_list = mat.pos.list)
  mat.neg <- check_and_merge(sparse_matrix_list = mat.neg.list)
  
  # for seurat obj
  pos.seurat.metadata$Cell_Barcode <- paste0(pos.seurat.metadata$Cell_Barcode, "_pos")
  
  # 4230 EL
  original.mtx <- mat.pos
  
  # for 4230 EL
  RNA.sample.ids <- unique(pos.seurat.metadata$Sample_Origin)
  RNA.sample.ids <- unique(paste0(str_split_fixed(RNA.sample.ids, "_", 3)[,1], "_", str_split_fixed(RNA.sample.ids, "_", 3)[,2]))
  
  # RNA.sample.tmp <- RNA.sample.ids[1]
  for (RNA.sample.tmp in RNA.sample.ids){
    
    # which sample?
    message(RNA.sample.tmp)
    
    # get sample metadata
    sub.sample.seurat.metadata <- pos.seurat.metadata[grep(RNA.sample.tmp, pos.seurat.metadata$Sample_Origin),]
    
    # reorder cells according to timepoints
    EL.order <- c(paste0(RNA.sample.tmp,"_", c("EARLY", "LATE")))
    sub.sample.seurat.metadata$Sample_Origin <- factor(sub.sample.seurat.metadata$Sample_Origin, levels = EL.order)
    sub.sample.seurat.metadata <- sub.sample.seurat.metadata[with(sub.sample.seurat.metadata,order(Sample_Origin)),]
    
    # create output folder for sample
    dir.create(paste0(w.dir, "/", RNA.sample.tmp))
    
    # subset original.mtx to sub.sample.seurat.metadata
    names.use <- colnames(original.mtx)[(colnames(original.mtx) %in% sub.sample.seurat.metadata$Cell_Barcode)]
    mat.tmp <- original.mtx[, names.use]
    
    # subset metadata to only the cells present in the matrix
    sub.sample.seurat.metadata <- sub.sample.seurat.metadata[sub.sample.seurat.metadata$Cell_Barcode %in% colnames(mat.tmp),]
    
    # merge matrices
    mat.matrix.sparse <- check_and_merge(sparse_matrix_list = list(mat.tmp, mat.neg))
    
    # check on the dimensions of the matrix
    # they will be displayed as a whole string
    print(dim(mat.tmp))
    print(dim(mat.neg))
    print(dim(mat.matrix.sparse))
    
    # make a cell annotation table for 6077_EL with ordered cells
    cellAnnotations <- data.frame(Cell_Barcode = colnames(mat.matrix.sparse))
    cellAnnotations <- merge(cellAnnotations, sub.sample.seurat.metadata[,c("Cell_Barcode", "Sample_Origin")], by = "Cell_Barcode", all.x = T)
    
    # set timepoint order
    cellAnnotations$Sample_Origin <- as.character(cellAnnotations$Sample_Origin)
    cellAnnotations$Sample_Origin <- factor(cellAnnotations$Sample_Origin, levels = EL.order)
    cellAnnotations <- cellAnnotations[with(cellAnnotations,order(Sample_Origin)),]
    
    # add "malignant" and "normal"
    cellAnnotations$Sample_Origin <- paste0("malignant_", cellAnnotations$Sample_Origin)
    cellAnnotations[cellAnnotations == "malignant_NA"] <- "normal"
    cellAnnotations <- distinct(cellAnnotations)
    
    # write df
    write.table(cellAnnotations, file= paste0(w.dir, "/", RNA.sample.tmp, "/cellAnnotations_", RNA.sample.tmp, ".txt"), col.names = F, row.names = F,quote=F, sep="\t")
    
    # remove non-unique barcodes
    mat.matrix.sparse <- mat.matrix.sparse[,unique(colnames(mat.matrix.sparse))]
    mat.matrix.sparse <- mat.matrix.sparse[,colnames(mat.matrix.sparse) %in% cellAnnotations$Cell_Barcode]
    
    # save merged sample matrix
    save(mat.matrix.sparse, file = paste0(w.dir, "/", RNA.sample.tmp, "/sc.10x.counts_", RNA.sample.tmp, ".RData"), compress = T)
    
    # make a gene ordering file
    if(TRUE){
      
      mart = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
      ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
                   "hgnc_symbol", rownames(mat.matrix.sparse), mart)
      gene_ordering_file <- ann[,c(1:4)]
      
      # add chr and subset to essential chromosomes
      gene_ordering_file <- gene_ordering_file[gene_ordering_file$chromosome_name %in% c(1:22), ]
      
      # order numerical
      chr_order <- c(1:22)
      gene_ordering_file$chromosome_name <- factor(gene_ordering_file$chromosome_name, levels = chr_order)
      gene_ordering_file <- gene_ordering_file[with(gene_ordering_file, order(chromosome_name, start_position)),]
      
      # write gene ordering file
      gene_ordering_file <- gene_ordering_file[-which(duplicated(gene_ordering_file$hgnc_symbol)),]
      write.table(gene_ordering_file, file= paste0(w.dir, "/", RNA.sample.tmp, "/gene_ordering_file.txt"), col.names = F, row.names = F,quote=F, sep="\t")
    }
  }
}
  

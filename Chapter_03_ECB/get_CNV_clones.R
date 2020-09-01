#############################################################################################################################
##                                                                                                                      
##  Infer copy number clones from inferCNV output                                                                                            
##                                                                                                                      
##  Date: 13 May 2020                                                                                                                   
##  
##  Author: Moritz Przybilla
##
##
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", "dendextend", "tidyverse", 
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "fishplot", "Signac", "BiocManager", 
                      "biomaRt", "httr", "ComplexHeatmap", "rstatix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

`%notin%` <- Negate(`%in%`)
#####################################################################################
# READ IN DATA OF INTEREST
#####################################################################################
setwd("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/")
dir.create("arm_level")

# get all the BayesNet matrices subclusters
BayesNet.files <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze", pattern = "infercnv.19_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", full.names = T, recursive = T)
# BayesNet.files <- BayesNet.files[-grep("g8|g18|g27", BayesNet.files)]
BayesNet.files <- BayesNet.files[grep("g18|g27", BayesNet.files)]
# BayesNet.files <- BayesNet.files[-c(1,2)]

# get all the BayesNet matrices cells
# BayesNet.files <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric", pattern = "infercnv.19_HMM_predHMMi6.hmm_mode-cells.Pnorm_0.5.repr_intensities.observations.txt", full.names = T, recursive = T)

# and the sample ids
sample.ids <- str_split_fixed(BayesNet.files, "/", 9)[,8]
# sample.ids[3] <- c("Sequencing8_C5T2R2_cr3_WT")
# sample.ids[c(19)] <- c("Sequencing5_cr3_WT")

# iterate over all samples 
i <- 1
for (i in 1:length(sample.ids)){
  
  # read in the observations file with the copy number states filtered by the BayesNet output from infercnv
  BayesNet.matrix <- as.matrix(read.table(BayesNet.files[i], header = T, sep = " "))
  print(dim(BayesNet.matrix))
  
  # sample
  sample.tmp <- sample.ids[i]
  message(sample.tmp)
  
  #####################################################################################
  # ADAPT HMM MATRIX TO GET AMP AND DEL CALLS
  #####################################################################################
  # read in the chromosome arm positions
  chr.arm.pos <- read.table("/labs/ccurtis2/mjprzy/correlation_test/chromosome_arm_positions_GRCH38.txt", header = T)
  chr.arm.pos$Chrom <- paste0("chr", chr.arm.pos$Chrom)
  chr.arm.pos.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)
  
  # list of genes in matrix
  genes.of.interest <- data.frame(hgnc_symbol = rownames(BayesNet.matrix))
  
  # for peer certificate
  set_config(config(ssl_verifypeer = 0L))
  
  # use biomart to get hgnc_symbols
  mart = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
  ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
               "hgnc_symbol", as.character(genes.of.interest$hgnc_symbol), mart)
  
  ann$chromosome_name <- paste0("chr", ann$chromosome_name)
  
  # remove duplicated gene names
  ann <- distinct(ann, ensembl_gene_id, .keep_all = TRUE)
  
  # melt df/matrix with exp values
  melted.BayesNet.matrix <- reshape2::melt(BayesNet.matrix)
  colnames(melted.BayesNet.matrix) <- c("hgnc_symbol", "barcode", "copy_number_state")
  melted.BayesNet.matrix <- merge(melted.BayesNet.matrix, ann, by = "hgnc_symbol")
  colnames(melted.BayesNet.matrix) <- c("hgnc_symbol", "barcode", "copy_number_state", "chr", "start", "end", "ensembl_gene_id")
  
  # make GRange object for infercnv genes
  inferCNV.genes <- makeGRangesFromDataFrame(melted.BayesNet.matrix, keep.extra.columns = T)
  inferCNV.genes$gene_length <- width(inferCNV.genes)
  
  # then merge each dataframe information by overlap
  # find the overlapping regions of the WGS segments with the genes
  inferCNV.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, inferCNV.genes)
  
  # get the essential columns and make a new dataframe
  chr.arm.inferCNV <- data.frame("chr_arm" = inferCNV.GR.merge$Idf,
                                 "hgnc_symbol" = inferCNV.GR.merge$hgnc_symbol,
                                 "Cell_barcode" = inferCNV.GR.merge$barcode,
                                 "copy_number_state" = inferCNV.GR.merge$copy_number_state,
                                 "gene_length" = inferCNV.GR.merge$gene_length,
                                 "chr_arm_length" = inferCNV.GR.merge$Length,
                                 stringsAsFactors = F)
  
  # convert each row
  chr.arm.inferCNV$chr_arm <- as.character(chr.arm.inferCNV$chr_arm)
  chr.arm.inferCNV$Cell_barcode <- as.character(chr.arm.inferCNV$Cell_barcode)
  
  # calculate the weighting factor for each gene
  # however, since not the complete arm is coding and covered with genes, 
  # we will calculate a within arm length of contributiong genes
  arm.gene.length <- aggregate(chr.arm.inferCNV[,c("gene_length")], by = list(chr.arm.inferCNV$chr_arm, chr.arm.inferCNV$Cell_barcode), FUN = sum)
  covered.gene.length <- distinct(arm.gene.length[,c(1,3)])
  
  # merge with chromsome arm length
  arm.covered.length <- merge(covered.gene.length, chr.arm.pos[,c(1,5)], by.x = "Group.1", by.y = "Idf")
  arm.covered.length$percentage <- round(arm.covered.length$x/arm.covered.length$Length, 3)
  
  # add the covered.arm.length to the chr.arm.inferCNV dataframe for the adapted weighting factor
  chr.arm.inferCNV <- merge(chr.arm.inferCNV, arm.covered.length, by.x = "chr_arm", by.y = "Group.1")
  
  # remove some uninteresting columns and rename them 
  chr.arm.order <- unique(chr.arm.pos$Idf)
  chr.arm.inferCNV$chr_arm <- factor(chr.arm.inferCNV$chr_arm, levels = chr.arm.order)
  chr.arm.inferCNV <- chr.arm.inferCNV[with(chr.arm.inferCNV,order(chr_arm)),]
  chr.arm.inferCNV <- chr.arm.inferCNV[,c(1,3:7,9)]
  colnames(chr.arm.inferCNV) <- c("chr_arm", "Cell_barcode", "copy_number_state", "gene_length", "chr_arm_length", "gene_covered_length", "gene_arm_percentage")
  
  # calculate weighting factor 
  chr.arm.inferCNV$weight_factor <- chr.arm.inferCNV$gene_length/chr.arm.inferCNV$gene_covered_length
  
  # multiply prob with weight factor
  chr.arm.inferCNV$weighted_copy_number_state <- chr.arm.inferCNV$copy_number_state * chr.arm.inferCNV$weight_factor
  
  # and sum it up to get the copy number state per arm
  arm.level.CN.state <- aggregate(chr.arm.inferCNV[,ncol(chr.arm.inferCNV)], by = list(chr.arm.inferCNV$chr_arm, chr.arm.inferCNV$Cell_barcode), FUN = sum)
  colnames(arm.level.CN.state) <- c("chr_arm", "Cell_barcode", "weighted_copy_number_state")
  
  # calculate standard deviation
  w.cns.sd <- sd(arm.level.CN.state$weighted_copy_number_state)
  
  # assign states again based on deviation from sd
  arm.level.CN.state$weighted_copy_number_state[arm.level.CN.state$weighted_copy_number_state > 1-w.cns.sd & arm.level.CN.state$weighted_copy_number_state < 1+w.cns.sd] <- 1
  arm.level.CN.state$weighted_copy_number_state[arm.level.CN.state$weighted_copy_number_state < 1-w.cns.sd] <- 0
  arm.level.CN.state$weighted_copy_number_state[arm.level.CN.state$weighted_copy_number_state > 1+w.cns.sd & arm.level.CN.state$weighted_copy_number_state < 2+w.cns.sd] <- 2
  arm.level.CN.state$weighted_copy_number_state[arm.level.CN.state$weighted_copy_number_state > 2+w.cns.sd] <- 3
  
  # reshape the melted dataframe into an matrix 
  inferCNV.arm.matrix <- reshape2::acast(arm.level.CN.state, Cell_barcode ~ chr_arm, value.var = "weighted_copy_number_state", fun.aggregate = mean)
  
  # make more intuitive by adding 1, so that diploid = 2 and arm x cells
  inferCNV.arm.matrix <- t(inferCNV.arm.matrix + 1)
  colnames(inferCNV.arm.matrix) <- str_split_fixed(colnames(inferCNV.arm.matrix), "_pos", 2)[,1]
  
  # write matrix
  write.table(inferCNV.arm.matrix, paste0("arm_level/",sample.tmp, "_arm_level_matrix.txt"), quote = F, sep = "\t", col.names = T, row.names = T)
  
  # Order cells with hierarchical clustering
  dist.centered.matrix <- dist(t(inferCNV.arm.matrix), method = "euclidean")
  hc <- hclust(dist.centered.matrix, method = "ward.D2")
  cellOrder <- hc$order
  
  # set parameters
  pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
  zlim <- c(0, 4) # upper and lower limit
  
  # order and have a sneak peak
  inferCNV.arm.matrix.clustered <- t(inferCNV.arm.matrix)[cellOrder,]
  png(paste0("arm_level/", sample.tmp, "_armlevel_heatmap.png"), width = 13, height = 4, units = 'in', res = 600, type = "cairo")
  par(mfrow=c(1,1), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
  image(t(as.matrix(inferCNV.arm.matrix.clustered)), col=pcol, zlim=zlim, useRaster = T)
  dev.off()
  
  #####################################################################################
  # GENERATE CLONES HERE PER CELL ASSIGNMENT
  #####################################################################################
  
  # remove chromosomes 21 and 22
  inferCNV.arm.matrix.clustered <- as.data.frame(inferCNV.arm.matrix.clustered[, colnames(inferCNV.arm.matrix.clustered) %notin% c("21q", "22q")])
  
  # replace numbers by symbols
  inferCNV.arm.matrix.clustered[which(inferCNV.arm.matrix.clustered == 2, arr.ind = T)] <- "d"
  inferCNV.arm.matrix.clustered[which(inferCNV.arm.matrix.clustered == 3, arr.ind = T)] <- "+"
  inferCNV.arm.matrix.clustered[which(inferCNV.arm.matrix.clustered == 1, arr.ind = T)] <- "-"
  inferCNV.arm.matrix.clustered[which(inferCNV.arm.matrix.clustered == 0, arr.ind = T)] <- "2-"
  inferCNV.arm.matrix.clustered[which(inferCNV.arm.matrix.clustered == 4, arr.ind = T)] <- "2+"
  
  # add clone column
  inferCNV.arm.matrix.clustered$cnv_clone <- "NA"
  
  # iterate over each cell and get the cnv profile
  for (i in 1:nrow(inferCNV.arm.matrix.clustered)){
    
    # genotype dummy
    this.genotype <- paste(paste0(colnames(inferCNV.arm.matrix.clustered)[c(1:(ncol(inferCNV.arm.matrix.clustered)-1))],
                                  inferCNV.arm.matrix.clustered[i,c(1:(ncol(inferCNV.arm.matrix.clustered)-1))]), collapse= ":")
    inferCNV.arm.matrix.clustered[i,"cnv_clone"] <- this.genotype # changed + to :
  }
  
  # this is a confusionn matrix with the timepoint against the genotype
  scRNA.tmp <- table(inferCNV.arm.matrix.clustered$cnv_clone)

  # rename the different clones
  cnv.clones <- unique(inferCNV.arm.matrix.clustered$cnv_clone)
  names(cnv.clones) <- paste0("clone", c(1:length(cnv.clones)))
  
  # make a new dataframe with cell id + cnv_id + clone_id
  df <- data.frame(Cell_barcodes = rownames(inferCNV.arm.matrix.clustered), cnv_id = inferCNV.arm.matrix.clustered$cnv_clone, clone_id = "clone", stringsAsFactors = F)
  
  # iterate over each clone and replace the clone_id name
  for (i in 1:length(cnv.clones)){
    df[which(df$cnv_id == cnv.clones[i]), "clone_id"] <- names(cnv.clones[i])
  }
  
  # get the significant clones
  cnv.clones <- data.frame(clone_id = names(cnv.clones), cnv_id = as.vector(cnv.clones), stringsAsFactors = F)
  scRNA.clones <- scRNA.tmp[scRNA.tmp > 10]
  
  # subset the matrix to the essential clones
  cnv.clones <- cnv.clones[cnv.clones$cnv_id %in% names(scRNA.clones),]
  
  # subset the cnv_id to the essential alterations
  cnv.clones$cnv_id <- str_remove_all(cnv.clones$cnv_id, pattern = "(\\d+)qd:")
  cnv.clones$cnv_id <- str_remove_all(cnv.clones$cnv_id, pattern = "(\\d+)pd:")
  cnv.clones$cnv_id <- str_remove_all(cnv.clones$cnv_id, pattern = ":(\\d+)qd")
  cnv.clones$cnv_id <- str_remove_all(cnv.clones$cnv_id, pattern = ":(\\d+)pd")
  
  # same for the df
  df$cnv_id <- str_remove_all(df$cnv_id, pattern = "(\\d+)qd:")
  df$cnv_id <- str_remove_all(df$cnv_id, pattern = "(\\d+)pd:")
  df$cnv_id <- str_remove_all(df$cnv_id, pattern = ":(\\d+)qd")
  df$cnv_id <- str_remove_all(df$cnv_id, pattern = ":(\\d+)pd")
  
  # write clones to files
  write.table(cnv.clones, paste0("arm_level/", sample.tmp, "_translationFile_cnv_clones.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
  write.table(df, paste0("arm_level/", sample.tmp, "_cnv_clones.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
}

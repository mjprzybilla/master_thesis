################################################################################################################################################
##                                                                                                                      
##  Prepare input for inferCNV                                                                                            
##                                                                                                                      
##  Date: 11 March 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
##  Summary: This script shall produce the inferCNV inputs to use the inferCNV object afterwards.   
##           
##                                                                                                                      
################################################################################################################################################
# clear workspace
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("beanplot", "mixtools", "pheatmap", "tidyverse", "zoo", "squash", "biomaRt",
                      "CONICSmat", "Rtsne", "scran", "stringr", "Matrix.utils", "Seurat", "reshape2", "readr", "stringr",
                      "biomaRt", "Matrix", "stringr", "GenomicRanges", "scater")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

############################################################################
##                    Load data and set working directory
############################################################################
# define working directory
w.dir <- "/labs/ccurtis2/mjprzy/correlation_test"

# read in normalized CONICSmat matrix
norm.gExp.matrix <- as.matrix(read_rds(paste0(w.dir, "/P2C2R2T2/CONICSmat.matrix.rds")))
colnames(norm.gExp.matrix) <- str_split_fixed(colnames(norm.gExp.matrix), "-", 2)[,1]

# get sample.idss
sample.ids <- c("P2C2R2T2")

# read in annotation file
annotation.file <- read.table("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/cellAnnotations_Sequencing8_C5T2R2_cr3_WT.txt", header = F, sep = "\t", stringsAsFactors = F)
annotation.file$V1 <- str_split_fixed(annotation.file$V1, "-", 2)[,1]
annotation.file <- annotation.file[annotation.file$V2 != "normal",]

# list of genes in matrix
genes.of.interest <- data.frame(hgnc_symbol = rownames(norm.gExp.matrix))

# for peer certificate
set_config(config(ssl_verifypeer = 0L))

# use biomart to get hgnc_symbols
mart = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
             "hgnc_symbol", as.character(genes.of.interest$hgnc_symbol), mart)

ann$chromosome_name <- paste0("chr", ann$chromosome_name)

# remove duplicated gene names
ann <- distinct(ann, ensembl_gene_id, .keep_all = TRUE)

# wrangle gene ordering file into shape
gene.ordering.file <- ann
gene.ordering.file$coordinates <- paste(gene.ordering.file$chromosome_name, gene.ordering.file$start_position, gene.ordering.file$end_position, sep = ":")
gene.ordering.file <- gene.ordering.file[gene.ordering.file$chromosome_name %in% paste0("chr", c(1:22)),]
rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol
colnames(gene.ordering.file) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id", "coordinates")

# read in metadata
seurat.obj <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/Sequencing8_C5T2R2_cr3/Sequencing8_C5T2R2_cr3_seurat_obj.rds")
metadata <- seurat.obj@meta.data
metadata <- metadata[,c(1,5,6,7)]

#####################################################################################
#                   Generate input for the image function to plot
#####################################################################################
# center the infercnv output to 0
gexp.norm <- norm.gExp.matrix

# only take genes which are present in the matrix
genes <- makeGRangesFromDataFrame(gene.ordering.file)

# generate new object genes.interest with genes and respective coordinates
genes.interest <- as.data.frame(genes)
rownames(genes.interest) <- names(genes)
mat <- gexp.norm

# check for each chromosome if there are genes present
gExp.array <- tapply(1:nrow(genes.interest),as.factor(genes.interest$seqnames),function(ii) {
  na.omit(mat[rownames(genes.interest)[ii[order((genes.interest[ii,]$start+genes.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
})

# subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
chrs <- paste0("chr", 1:22)
gExp.array <- gExp.array[chrs]

# we want to scale the output to the respective amount of genes and their widths
widths <- sapply(gExp.array, nrow); widths <- widths/max(widths)*100

# generate the layout accordingly
plot.layout <- layout(matrix(seq(1,length(gExp.array)),1,length(gExp.array),byrow=TRUE), widths=widths)

# we also want to order or image hierarchically
setOrder <- T
if(setOrder) {
  avgd <- do.call(rbind, lapply(names(gExp.array),function(nam) {
    d <- gExp.array[[nam]]
    d <- colMeans(d)
    d
  }))
  hc <- hclust(dist(t(avgd)))
  cellOrder <- hc$order
}

# set parameters for the image
adapted.gExp.list <- gExp.array
pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
zlim <- c(-0.65, 0.65) # upper and lower limit

# limit the gExp to maximums according to zlim
limit.gExp.list <- lapply(names(adapted.gExp.list),function(nam) {
  d <- adapted.gExp.list[[nam]]
  d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
  return(d)
})

#####################################################################################
#                                 Make the plot
#####################################################################################

png("/labs/ccurtis2/mjprzy/correlation_test/P2C2R2T2/CONICSmat_P2C2R2T2_heatmap.png" , width = 13, height = 4, units = 'in', res = 600)
par(mfrow=c(1,22), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
## plot chromosomes
box()
for (i in 1:length(limit.gExp.list)){
  message(paste0("chr", i))
  d <- limit.gExp.list[[i]]
  d <- d[, cellOrder] 
  # image the respective chr 
  image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
  box()
}
dev.off()

png("/labs/ccurtis2/mjprzy/correlation_test/P2C2R2T2/CONICSmat_P2C2R2T2_heatmap_legend.png"  , width = 13, height = 4, units = 'in', res = 600, type = "cairo")
par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
dev.off()


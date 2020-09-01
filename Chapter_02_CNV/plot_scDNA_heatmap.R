#############################################################################################################################
##                                                                                                                      
##  Create Heatmap visualization for scCNV data output from construct_ploidyMatrix_scDNA.R                                                                                    
##                                                                                                                      
##  Date: 20 April 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
############################################################################################################################

# clear workspace
rm(list=ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if nnot available
list.of.packages <- c("tidyverse", "Seurat", "biomaRt", "RColorBrewer", "hdf5r", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# define the input directory containing the seurat objects
w.dir <- "/labs/ccurtis2/mjprzy/correlation_test/"

# get sample name
sample.tmp <- "6077_P_C5_ECB_R2_T2"
message(sample.tmp)

###########################################################################
#                   READ IN DATA OF INTEREST FROM inferCNV
###########################################################################
## Read in scCNV data for 6077 C5 R2 T2

# read in ploidy matrix from scDNA (construct_ploidyMatrix_scDNA.R)
scDNA.matrix <- read.table("/labs/ccurtis2/mjprzy/correlation_test/cnv.cells.mat.txt", header = T, stringsAsFactors = F)

# get dimensions
print(dim(scDNA.matrix))

# get the segments and make GRange from it 
seg.coordinates <- data.frame(chr = paste0("chr", str_split_fixed(scDNA.matrix$id.pos, ":",2)[,1]), 
                              start = str_split_fixed(str_split_fixed(scDNA.matrix$id.pos, ":",2)[,2], "-", 2)[,1], 
                              end = str_split_fixed(str_split_fixed(scDNA.matrix$id.pos, ":",2)[,2], "-", 2)[,2])
seg.coordinates.GRange <- makeGRangesFromDataFrame(seg.coordinates, keep.extra.columns = T)

# remove id column
rownames(scDNA.matrix) <- paste(seg.coordinates$chr, seg.coordinates$start, seg.coordinates$end, sep = "-")
scDNA.matrix <- scDNA.matrix[,-1]

#####################################################################################
# GENERATE A HEATMAP FOR THE LOG2FC
#####################################################################################

# only take features which are present in the matrix
# windowSummary is an grange object with windows names as rownames, the chromsome as seqnames and ranges
# convert back to dataframe again with genes as rownames seqnames start end width
windows.of.interest <- as.data.frame(seg.coordinates.GRange)
rownames(windows.of.interest) <- paste(windows.of.interest$seqnames, windows.of.interest$start, windows.of.interest$end, sep = "-")
mat <- as.matrix(scDNA.matrix)

# check for each chromosome if there are features present
# the matrix is splitted into a large array with features x cells 
tl <- tapply(1:nrow(windows.of.interest), as.factor(windows.of.interest$seqnames), function(ii) {
  na.omit(mat[rownames(windows.of.interest)[ii[order((windows.of.interest[ii,]$start+windows.of.interest[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
})

# subset to chromosome 3,4,9,13 and 20 (chr of interest in general)
chrs <- paste0("chr", c(1:22))
tl <- tl[chrs] # only these will be still there

# we want to scale the output to the respective amount of genes and their widths
# so here, the widths for each window are determined
setWidths <- T
chromSizes <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:24]
if(setWidths) {
  # widths <- sapply(tl, nrow); widths <- widths/max(widths)*100
  # Can also set by known chromosome size widths:
  widths <- end(chromSizes)[1:24]/1e7
} else {
  widths <- rep(1, length(tl))
}

# generate the layout accordingly
l <- layout(matrix(seq(1,length(tl)),1,length(tl),byrow=TRUE), widths=widths)

#####################################################################################
# PERFORM HIERARCHICAL CLUSTERING ON THE MATRIX TO ORDER THE CELLS
#####################################################################################
# calculate average distance
# adapted from Jean Fan fot the code https://github.com/JEFworks/HoneyBADGER
avgd <- do.call(rbind, lapply(names(tl),function(nam) {
  d <- tl[[nam]]
  d <- colMeans(d)
  d
}))

# Order cells with hierarchical clustering
dist.centered.matrix <- dist(t(avgd), method = "euclidean")
hc <- hclust(dist.centered.matrix, method = "ward.D2")
cellOrder <- hc$order

#####################################################################################
# GENERATE A HEATMAP FOR THE LOG2FC
#####################################################################################

# set parameters for the image
tlsub <- tl
pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
zlim <- c(0, 4) # upper and lower limit

# limit the gExp to maximums according to zlim
tlsmooth <- lapply(names(tlsub),function(nam) {
  d <- tlsub[[nam]]
  d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
  return(d)
})

png("/labs/ccurtis2/mjprzy/10x_scCNV_heatmap.png" , width = 13, height = 4, units = 'in', res = 600)
par(mfrow=c(1,22), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
for (i in 1:length(tlsmooth)){
  message(paste0("chr", i))
  d <- tlsmooth[[i]]
  d <- d[, cellOrder]
  # image the respective chr 
  image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=paste0("chr", i), useRaster = T)
  box()
}
dev.off()

png("/labs/ccurtis2/mjprzy/10x_scCNV_heatmap_legend.png", width = 13, height = 4, units = 'in', res = 600, type = "cairo")
par(mfrow=c(1,24), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.5)
image.plot(legend.only=TRUE, zlim=zlim, col = pcol, legend.shrink = 0.4, legend.width = 10, smallplot=  c(.65, .9, .5, .8)) # first argument of small plot = width of legend, third is making the length, the higher the smaller
dev.off()

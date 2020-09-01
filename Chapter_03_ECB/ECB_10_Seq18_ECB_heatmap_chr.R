################################################################################################################################################
##                                                                                                                      
##  Plot zoom-in on specific chromosomes from inferCNV for Sequecing 18                                                                 
##                                                                                                                      
##  Date: 11 April 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                    
##                                                                                                                      
##                                                                                                                      
################################################################################################################################################
# clear workspace beforehand
rm(list = ls())
set.seed(2020)

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("biomaRt", "Matrix", "stringr", "GenomicRanges", "scater")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
#                                  load data
#####################################################################################

# list all necessary information from infercnv
Seq18.norm.gExp.matrices <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing18_D3_APA4_Early_WT", pattern = "infercnv.median_filtered.observations.txt", recursive = T, full.names = T)
Seq18.annotation.files <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing18_D3_APA4_Early_WT", pattern = "cellAnnotations", recursive = T, full.names = T)
Seq18.gene.ordering.files <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing18_D3_APA4_Early_WT", pattern = "gene_ordering_file.txt", recursive = T, full.names = T)

# read in metadata
seurat.obj <- readRDS("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/Sequencing18_D3_APA4_Early/Sequencing18_D3_APA4_Early_seurat_obj.rds")
metadata <- seurat.obj@meta.data
metadata <- metadata[,c(1,5,6,7)]

i <- 2
for (i in 1:length(Seq18.norm.gExp.matrices)){
  
  # get sample id
  sample.tmp <- str_split_fixed(Seq18.norm.gExp.matrices[i], "/", 9)[,8]
  print(sample.tmp)
  
  # read in normalized matrix
  norm.gExp.matrix <- as.matrix(read.table(Seq18.norm.gExp.matrices[i], header = T, sep = " "))
  colnames(norm.gExp.matrix) <- str_split_fixed(colnames(norm.gExp.matrix), "\\.", 2)[,1]
  
  # read in annotation file
  annotation.file <- read.table(Seq18.annotation.files[i], header = F, sep = "\t", stringsAsFactors = F)
  annotation.file$V1 <- str_split_fixed(annotation.file$V1, "-", 2)[,1]
  annotation.file <- annotation.file[annotation.file$V2 == "malignant_0",]

  # read in gene_ordering file
  gene.ordering.file <- read.table(Seq18.gene.ordering.files[i], header = F, col.names = c("hgnc_symbol", "chr", "start", "end"), sep = "\t")
  gene.ordering.file$coordinates <- paste(paste0("chr", gene.ordering.file$chr), gene.ordering.file$start, gene.ordering.file$end, sep = ":")
  gene.ordering.file$chr <- paste0("chr", gene.ordering.file$chr)
  rownames(gene.ordering.file) <- gene.ordering.file$hgnc_symbol
  
  # subset metadata
  sample.metadata <- metadata[metadata$HashTag == sample.tmp,]
  
  # subset norm.gExp.matrix by annotation.file
  colnames(norm.gExp.matrix) <- str_split_fixed(colnames(norm.gExp.matrix), "\\.", 2)[,1]
  ECB.zero.norm.gExp.matrix <- norm.gExp.matrix[, colnames(norm.gExp.matrix) %in% annotation.file$V1]
  
  # make GRange object
  genes <- makeGRangesFromDataFrame(gene.ordering.file)
  
  #####################################################################################
  #                      Generate image plot for the matrix
  #####################################################################################
  # center the infercnv output to 0
  gexp.norm <- ECB.zero.norm.gExp.matrix -1
  
  # only take genes which are present in the matrix
  genes <- genes[rownames(gexp.norm)]
  
  # generate new object gos with genes and respective coordinates
  gos <- as.data.frame(genes)
  rownames(gos) <- names(genes)
  mat <- gexp.norm
  
  # check for each chromosome if there are genes present
  tl <- tapply(1:nrow(gos),as.factor(gos$seqnames),function(ii) {
    na.omit(mat[rownames(gos)[ii[order((gos[ii,]$start+gos[ii,]$end)/2,decreasing=F)]],,drop=FALSE])
  })
  
  # subset to chromosome 4 and 20 (chr of interest in general)
  chrs <- c("chr3", "chr9", "chr13", "chr16")
  tl <- tl[chrs]
  
  # we want to scale the output to the respective amount of genes and their widths
  setWidths <- T
  if(setWidths) {
    widths <- sapply(tl, nrow); widths <- widths/max(widths)*100
    ##Can also set by known chromosome size widths:
    ##https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
    ##widths <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 51304566, 48129895)/1e7
  } else {
    widths <- rep(1, length(tl))
  }
  
  # generate the layout accordingly
  l <- layout(matrix(seq(1,length(tl)),1,length(tl),byrow=TRUE), widths=widths)
  
  # we also want to order or image hierarchically
  setOrder <- T
  if(setOrder) {
    avgd <- do.call(rbind, lapply(names(tl),function(nam) {
      d <- tl[[nam]]
      d <- colMeans(d)
      d
    }))
    hc <- hclust(dist(t(avgd)))
    cellOrder <- hc$order
  }
  
  # set parameters for the image
  tlsub <- tl
  pcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256) # color panel
  window.size <- 31 # window size for smoothing
  zlim <- c(-0.2, 0.2) # upper and lower limit
  
  # limit the gExp to maximums according to zlim
  tlsmooth <- lapply(names(tlsub),function(nam) {
    d <- tlsub[[nam]]
    d[d< zlim[1]] <- zlim[1]; d[d>zlim[2]] <- zlim[2];
    return(d)
  })
  
  png(paste0("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing18_D3_APA4_Early_WT/", sample.tmp, "/", sample.tmp, "_w_ECB_chr3_9_13_16.png") , width = 10, height = 2.5, units = 'in', res = 600)
  par(mfrow=c(1,4), mar = c(0.5,0.2,3.0,0.2), mgp = c(2,0.65,0), cex = 0.8)
  ## plot chromosomes
  box()
  for (i in 1:length(tlsmooth)){
    message(chrs[i])
    d <- tlsmooth[[i]]
    d <- d[, cellOrder] 
    # image the respective chr 
    image(seq_len(nrow(d)), seq_len(ncol(d)), d, col=pcol, zlim=zlim, xlab="", ylab="", axes=F, main=chrs[i], useRaster = T)
    box()
  }
  dev.off()
}

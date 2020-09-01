#############################################################################################################################
##                                                                                                                      
##  CORRELATION ANALYSIS BETWEEN PSEUDOBULK MODIFIED EXPRESSION FROM INFERCNV AND CONICSMAT AND SHALLOW WGS 
##                                                                                                                      
##  Date: 18 March 2020                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                      
############################################################################################################################
# clear workspace
rm(list=ls())
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("tidyverse", "ggplot2", "GenomicRanges", "Matrix", "Signac", "Seurat",
                      "patchwork", "reshape2", "matrixStats", "biomaRt", "ggpubr", "RColorBrewer",
                      "viridis", "cowplot", "grid", "ggplotify")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# define output directory
o.dir <- "/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/"
dir.create(o.dir)

# set working diretory
setwd(o.dir)

#######################################################################################                                                                                                                    
##                      READ IN THE PREPROCESSED DATA                                                                                           
#######################################################################################

# read in the dataf.frames for each individual method 
inferCNV.data <- read.table("inferCNV_complete_df.txt", header = T)
CONICS.data <- read.table("CONICS_complete_df.txt", header = T)

# read in the chromosome arm positions
chr.arm.pos <- read.table("chr_arm_position_grch38.txt", header = T, stringsAsFactors = F)
chr.arm.pos$Chrom <- paste0("chr", chr.arm.pos$Chrom)
chr.arm.pos.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)

#######################################################################################                                                                                                                    
##           SPLIT DATA INTO COPY NUMBER AMP AND DEL - INFERCNV
#######################################################################################

# read in the dataframes for each individual method 
inferCNV.GRange <- makeGRangesFromDataFrame(inferCNV.data, keep.extra.columns = T)
CONICS.GRange <- makeGRangesFromDataFrame(CONICS.data, keep.extra.columns = T)

# then merge each dataframe information by overlap
# find the overlapping regions of the WGS segments with the genes
inferCNV.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, inferCNV.GRange)
CONICS.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, CONICS.GRange)

# subset it to altered chromosomes
subset.inferCNV.GR.merge <- inferCNV.GR.merge[inferCNV.GR.merge$Idf %in% c("3p", "3q", "4p", "4q", "9p", "9q",  "13p", "13q", "20q"),]
subset.CONICS.GR.merge <- CONICS.GR.merge[CONICS.GR.merge$Idf %in%c("3p", "3q", "4p", "4q", "9p", "9q", "13p", "13q", "20q"),]

## AMPLIFICATIONS
# copy number gain genes infercnv
amp.inferCNV.GR.subset <- subset.inferCNV.GR.merge[subset.inferCNV.GR.merge$Idf %in% c("3q", "9q", "20q"),]

# get the essential columns and make a new dataframe
amp.inferCNV <- data.frame("hgnc_symbol" = amp.inferCNV.GR.subset$hgnc_symbol,
                           "Cell_barcode" = amp.inferCNV.GR.subset$barcode,
                           "exp_value" = amp.inferCNV.GR.subset$exp_value,
                           "mean_logR" = amp.inferCNV.GR.subset$mean_logR,
                           "chr" = amp.inferCNV.GR.subset$Idf,
                           "gene_length" = amp.inferCNV.GR.subset$gene_length,
                           stringsAsFactors = F)

amp.inferCNV <- aggregate(amp.inferCNV[,c(3,4,6)], list(amp.inferCNV$hgnc_symbol, amp.inferCNV$chr), mean)
colnames(amp.inferCNV) <- c("hgnc_symbol", "chr", "exp_value", "mean_logR", "gene_length")
amp.inferCNV$chr <- factor(amp.inferCNV$chr, levels = c("3p", "3q", "4p", "4q", "9p", "9q",  "13p", "13q", "20q"))

## DELETIONS
# copy number loss genes infercnv
del.inferCNV.GR.subset <- subset.inferCNV.GR.merge[subset.inferCNV.GR.merge$Idf %in% c("3p", "4p", "4q", "9p", "13p", "13q") ,]

# get the essential columns and make a new dataframe
del.inferCNV <- data.frame("hgnc_symbol" = del.inferCNV.GR.subset$hgnc_symbol,
                           "Cell_barcode" = del.inferCNV.GR.subset$barcode,
                           "exp_value" = del.inferCNV.GR.subset$exp_value,
                           "mean_logR" = del.inferCNV.GR.subset$mean_logR,
                           "chr" = del.inferCNV.GR.subset$Idf,
                           "gene_length" = del.inferCNV.GR.subset$gene_length,
                           stringsAsFactors = F)

del.inferCNV <- aggregate(del.inferCNV[,c(3,4,6)], list(del.inferCNV$hgnc_symbol, del.inferCNV$chr), mean)
colnames(del.inferCNV) <- c("hgnc_symbol", "chr", "exp_value", "mean_logR", "gene_length")
del.inferCNV$chr <- factor(del.inferCNV$chr, levels = c("3p", "3q", "4p", "4q", "9p", "9q",  "13p", "13q", "20q"))

#######################################################################################                                                                                                                    
##           SPLIT DATA INTO COPY NUMBER AMP AND DEL - CONICSMAT
#######################################################################################
## AMPLIFICATIONS
# copy number gain genes conics
amp.CONICS.GR.subset <- subset.CONICS.GR.merge[subset.CONICS.GR.merge$Idf %in% c("3q", "9q", "20q"),]

# get the essential columns and make a new dataframe
amp.CONICS <- data.frame("hgnc_symbol" = amp.CONICS.GR.subset$hgnc_symbol,
                         "Cell_barcode" = amp.CONICS.GR.subset$barcode,
                         "exp_value" = amp.CONICS.GR.subset$exp_value,
                         "mean_logR" = amp.CONICS.GR.subset$mean_logR,
                         "chr" = amp.CONICS.GR.subset$Idf,  
                         "gene_length" = amp.CONICS.GR.subset$gene_length,
                         stringsAsFactors = F)

amp.CONICS <- aggregate(amp.CONICS[,c(3,4,6)], list(amp.CONICS$hgnc_symbol, amp.CONICS$chr), mean)
colnames(amp.CONICS) <- c("hgnc_symbol", "chr", "exp_value", "mean_logR", "gene_length")
amp.CONICS$chr <- factor(amp.CONICS$chr, levels = c("3p", "3q", "4p", "4q", "9p", "9q",  "13p", "13q", "20q"))

## DELETIONS
# copy number loss genes conics
del.CONICS.GR.subset <- subset.CONICS.GR.merge[subset.CONICS.GR.merge$Idf %in% c("3p", "4p", "4q", "9p",  "13p", "13q"),]

# get the essential columns and make a new dataframe
del.CONICS <- data.frame("hgnc_symbol" = del.CONICS.GR.subset$hgnc_symbol,
                         "Cell_barcode" = del.CONICS.GR.subset$barcode,
                         "exp_value" = del.CONICS.GR.subset$exp_value,
                         "mean_logR" = del.CONICS.GR.subset$mean_logR,
                         "chr" = del.CONICS.GR.subset$Idf,
                         "gene_length" = del.CONICS.GR.subset$gene_length,
                         stringsAsFactors = F)

del.CONICS <- aggregate(del.CONICS[,c(3,4,6)], list(del.CONICS$hgnc_symbol, del.CONICS$chr), mean)
colnames(del.CONICS) <- c("hgnc_symbol", "chr", "exp_value", "mean_logR", "gene_length")
del.CONICS$chr <- factor(del.CONICS$chr, levels = c("3p", "3q", "4p", "4q", "9p", "9q",  "13p", "13q", "20q"))

#######################################################################################                                                                                                                    
##                      Plot splitted data for each method                                                                                        
#######################################################################################
# inferCNV
ggscatter(amp.inferCNV, x = "exp_value", y = "mean_logR", size = "gene_length", color = "chr", alpha = 0.4,
          add = "reg.line", cor.method = "spearman", cor.coef = TRUE,
          conf.int = TRUE, add.params = list(color = "black", fill = "darkgray"),
          xlab = "Modified pseudobulk expression", ylab = "logR per Gene") +
  theme(legend.position="right") + labs(size = "Gene size [bp]", color = "Chromosome") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"))
ggsave("inferCNV_amp_genes_scatterplot.pdf", width = 9, height = 7.5)

ggscatter(del.inferCNV, x = "exp_value", y = "mean_logR", size = "gene_length", color = "chr", alpha = 0.4,
          add = "reg.line", cor.method = "spearman", cor.coef = TRUE,
          conf.int = TRUE, add.params = list(color = "black", fill = "darkgray"),
          xlab = "Modified pseudobulk expression", ylab = "logR per Gene") + theme(legend.position="right") +
  labs(size = "Gene size [bp]", color = "Chromosome") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"))
ggsave("inferCNV_del_genes_scatterplot.pdf", width = 9, height = 7.5)

# CONICS
ggscatter(amp.CONICS, x = "exp_value", y = "mean_logR", size = "gene_length", color = "chr", alpha = 0.4,
          add = "reg.line", cor.method = "spearman", cor.coef = TRUE,
          conf.int = TRUE, add.params = list(color = "black", fill = "darkgray"),
          xlab = "Modified pseudobulk expression", ylab = "logR per Gene") + theme(legend.position="right") +
  labs(size = "Gene size [bp]", color = "Chromosome") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"))
ggsave("CONICS_amp_genes_scatterplot.pdf", width = 9, height = 7.5)

ggscatter(del.CONICS, x = "exp_value", y = "mean_logR", size = "gene_length", color = "chr", alpha = 0.4,
          add = "reg.line", cor.method = "spearman", cor.coef = TRUE,
          conf.int = TRUE, add.params = list(color = "black", fill = "darkgray"),
          xlab = "Modified pseudobulk expression", ylab = "logR per Gene") + theme(legend.position="right") +
  labs(size = "Gene size [bp]", color = "Chromosome") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"))
ggsave("CONICS_del_genes_scatterplot.pdf", width = 9, height = 7.5)

#' -1 indicates a strong negative correlation : this means that every time x increases, y decreases (left panel figure)
#' 0 means that there is no association between the two variables (x and y) (middle panel figure)
#' 1 indicates a strong positive correlation : this means that y increases with x (right panel figure)


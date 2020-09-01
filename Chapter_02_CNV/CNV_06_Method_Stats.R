#############################################################################################################################
##                                                                                                                      
##  ASSESS STANDARD STATISTICS FROM THE DISTINCT CNV ESTIMATION METHODS HONEYBADGER, INFERCNV AND CONICSMAT                                                                               
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
                      "viridis")
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
##                    READ IN AND PREPROCESS SHALLOW WGS FILES                                                                                          
#######################################################################################
# read in CNV calls from WGS
wgs.cnv.table <- read.table("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/WGS_6077_P_C5_Rep2_T3_sorted_CNV_Smooth_Segment_per_bin.txt", header = T, stringsAsFactors = F)

# remove feature column and rename last one
wgs.cnv.table$feature <- NULL
colnames(wgs.cnv.table) <- c("chr", "start", "end", "LogR")

# add "chr" to chr column
wgs.cnv.table$chr <- paste0("chr", wgs.cnv.table$chr)

# make data frame a Grange object
CNV.GRange <- makeGRangesFromDataFrame(wgs.cnv.table, keep.extra.columns = T)
CNV.GRange$length <- width(CNV.GRange)

# read in chromosome arm segments
chr.arm.pos <- read.table("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/chr_arm_position_grch38.txt", header = T, stringsAsFactors = F)

# make GRange object
chr.arm.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)

#######################################################################################                                                                                                                    
##                     READ IN AND PREPROCESS HONEYBADGER                                                                                   
#######################################################################################
# read in the hb gene expression matrix
hb_exp.matrix <- readRDS("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/honeybadger_NormMatrix.rds")

# get dimensions
print(dim(hb_exp.matrix))

# remove wt cells
names.use <- grep("_pos", colnames(hb_exp.matrix), value = T)
tumor.hb.exp.matrix <- hb_exp.matrix[,names.use]

# get dimensions after removal of normal cells
print(dim(tumor.hb.exp.matrix))

# adapt barcodes
colnames(tumor.hb.exp.matrix) <- str_split_fixed(colnames(tumor.hb.exp.matrix), "-", 2)[,1]

# melt df/matrix with exp values
melted.tumor.hb.exp.matrix <- melt(tumor.hb.exp.matrix)
colnames(melted.tumor.hb.exp.matrix) <- c("hgnc_symbol", "barcode", "exp_value")

# list of genes in matrix
genes.of.interest <- data.frame(hgnc_symbol = rownames(tumor.hb.exp.matrix))

# use biomart to get hgnc_symbols
mart = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
hb.ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
             "hgnc_symbol", as.character(genes.of.interest$hgnc_symbol), mart)

hb.ann$chromosome_name <- paste0("chr", hb.ann$chromosome_name)

# remove duplicated gene names
hb.ann <- distinct(hb.ann, ensembl_gene_id, .keep_all = TRUE)
hb.ann <- hb.ann[hb.ann$chromosome_name %in% paste0("chr", c(1:22)),]

# remove duplicated gene names
hb.ann <- distinct(hb.ann, hgnc_symbol, .keep_all = TRUE)
rownames(hb.ann) <- hb.ann$hgnc_symbol
colnames(hb.ann) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id")

# make GRange object
hb.GRange <- makeGRangesFromDataFrame(hb.ann, keep.extra.columns = T)

# add length column to metadata
hb.GRange$length <- width(hb.GRange)

#######################################################################################                                                                                                              
##                     READ IN AND PREPROCESS CONICSMAT                                                                                       
#######################################################################################
# read in CONICSmat gene expression matrix
CONICS_exp.matrix <- readRDS("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/CONICSmat.matrix.rds")

# get dimensions
print(dim(CONICS_exp.matrix))

# melt df/matrix with exp values
melted.tumor.CONICS.exp.matrix <- melt(CONICS_exp.matrix)
colnames(melted.tumor.CONICS.exp.matrix) <- c("hgnc_symbol", "barcode", "exp_value")

# list of genes in matrix
genes.of.interest <- data.frame(hgnc_symbol = rownames(CONICS_exp.matrix))

# use biomart to get hgnc_symbols
CONICSmat.ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
                "hgnc_symbol", as.character(genes.of.interest$hgnc_symbol), mart)
CONICSmat.ann$chromosome_name <- paste0("chr", CONICSmat.ann$chromosome_name)

# remove duplicated gene names
CONICSmat.ann <- distinct(CONICSmat.ann, ensembl_gene_id, .keep_all = TRUE)
CONICSmat.ann <- CONICSmat.ann[CONICSmat.ann$chromosome_name %in% paste0("chr", c(1:22)),]

# remove duplicated gene names
CONICSmat.ann <- distinct(CONICSmat.ann, hgnc_symbol, .keep_all = TRUE)
rownames(CONICSmat.ann) <- CONICSmat.ann$hgnc_symbol
colnames(CONICSmat.ann) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id")

# convert to GRange object
CONICS.GRange <- makeGRangesFromDataFrame(CONICSmat.ann, keep.extra.columns = T)

# add length column to metadata
CONICS.GRange$length <- width(CONICS.GRange)

#######################################################################################                                                                                                                     
##                      READ IN AND PREPROCESS INFERCNV                                                                                          
#######################################################################################

# get gene expression matrix
inferCNV_exp.matrix <- as.matrix(read.table("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/infercnv.median_filtered.observations.txt", header = T, sep = " "))

# get dimensions
print(dim(inferCNV_exp.matrix))

# adapt barcodes
colnames(inferCNV_exp.matrix) <- str_split_fixed(colnames(inferCNV_exp.matrix), "\\.1_pos", 2)[,1]

# melt df/matrix with exp values
melted.inferCNV.exp.matrix <- melt(inferCNV_exp.matrix)
colnames(melted.inferCNV.exp.matrix) <- c("hgnc_symbol", "barcode", "exp_value")
melted.inferCNV.exp.matrix$exp_value <- as.numeric(as.character(melted.inferCNV.exp.matrix$exp_value)) - 1

# list of genes in matrix
genes.of.interest <- data.frame(hgnc_symbol = rownames(inferCNV_exp.matrix))

# use biomart to get hgnc_symbols
inferCNV.ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
                       "hgnc_symbol", as.character(genes.of.interest$hgnc_symbol), mart)
inferCNV.ann$chromosome_name <- paste0("chr", inferCNV.ann$chromosome_name)

# remove duplicated gene names
inferCNV.ann <- distinct(inferCNV.ann, ensembl_gene_id, .keep_all = TRUE)
inferCNV.ann <- inferCNV.ann[inferCNV.ann$chromosome_name %in% paste0("chr", c(1:22)),]

# remove duplicated gene names
inferCNV.ann <- distinct(inferCNV.ann, hgnc_symbol, .keep_all = TRUE)
rownames(inferCNV.ann) <- inferCNV.ann$hgnc_symbol
colnames(inferCNV.ann) <- c("hgnc_symbol", "chr", "start", "end", "ensembl_gene_id")

# make data frame a Grange object
inferCNV.GRange <- makeGRangesFromDataFrame(inferCNV.ann, keep.extra.columns = T)

# add length column to metadata
inferCNV.GRange$length <- width(inferCNV.GRange)

#######################################################################################                                                                                                                     
##                 CALCULATE LOGR AVERAGES PER GENE FROM WGS
#######################################################################################

# find the overlapping regions of the WGS segments with the genes
inferCNV.GR.merge <- mergeByOverlaps(inferCNV.GRange, CNV.GRange)
CONICS.GR.merge <- mergeByOverlaps(CONICS.GRange, CNV.GRange)
hb.GR.merge <- mergeByOverlaps(hb.GRange, CNV.GRange)

# calculate mean logR for each gene (hb)
hb.logR <- data.frame("hgnc_symbol" = rownames(hb.GR.merge),
                      "logR" = hb.GR.merge$LogR,
                      "gene_length" = hb.GR.merge$length,
                      stringsAsFactors = F)

hb.logR <- aggregate(hb.logR[,c(2,3)], list(hb.logR$hgnc_symbol), mean)
colnames(hb.logR) <- c("hgnc_symbol", "mean_logR", "gene_length")

# calculate mean logR for each gene (inferCNV)
inferCNV.logR <- data.frame("hgnc_symbol" = rownames(inferCNV.GR.merge),
                            "logR" = inferCNV.GR.merge$LogR,
                            "gene_length" = inferCNV.GR.merge$length,
                            stringsAsFactors = F)

inferCNV.logR <- aggregate(inferCNV.logR[,c(2,3)], list(inferCNV.logR$hgnc_symbol), mean)
colnames(inferCNV.logR) <- c("hgnc_symbol", "mean_logR", "gene_length")

# calculate mean logR for each gene (CONICSmat)
CONICS.logR <- data.frame("hgnc_symbol" = rownames(CONICS.GR.merge),
                          "logR" = CONICS.GR.merge$LogR,
                          "gene_length" = CONICS.GR.merge$length,
                          stringsAsFactors = F)

CONICS.logR <- aggregate(CONICS.logR[,c(2,3)], list(CONICS.logR$hgnc_symbol), mean)
colnames(CONICS.logR) <- c("hgnc_symbol", "mean_logR", "gene_length")

#######################################################################################                                                                                                                     
##                      PLOT EXPRESSION/LOGR VALUES PER GENE
#######################################################################################
# plot WGS logR
plotWGS <- ggplot(wgs.cnv.table, aes(LogR)) + 
  geom_density(alpha=0.2) + 
  labs(fill="Stats", 
       title = "WGS",
       x="LogR",
       y="Density") + 
  xlim(-3, +3) +
  theme_bw()

# plot exp_value distribution 
plothb1 <- ggplot(melted.tumor.hb.exp.matrix, aes(exp_value)) + 
  geom_density(alpha=0.2) + 
  labs(fill="Stats",
       title = "HoneyBADGER",
       x="Normalized Expression",
       y="Density") + 
  xlim(-3, +3) +
  theme_classic() + theme(plot.title=element_text(size=20, 
                                                  face="bold", 
                                                  hjust=0.5,
                                                  lineheight=1.2),  # title
                          axis.title.x=element_text(size=15, face = "bold"),  # X axis title
                          axis.title.y=element_text(size=15, face = "bold"),  # Y axis title
                          axis.text.x=element_text(size=10, 
                                                   angle = 30,
                                                   vjust=.5),  # X axis text
                          axis.text.y=element_text(size=10))  # Y axis text

plotCONICS1 <- ggplot(melted.tumor.CONICS.exp.matrix, aes(exp_value)) + 
  geom_density(alpha=0.2) + 
  labs(fill="Stats", 
       title = "CONICSmat",
       x="Normalized Expression",
       y="Density") + 
  xlim(-3, +3) +
  theme_classic() + theme(plot.title=element_text(size=20, 
                                                  face="bold", 
                                                  hjust=0.5,
                                                  lineheight=1.2),  # title
                          axis.title.x=element_text(size=15, face = "bold"),  # X axis title
                          axis.title.y=element_text(size=15, face = "bold"),  # Y axis title
                          axis.text.x=element_text(size=10, 
                                                   angle = 30,
                                                   vjust=.5),  # X axis text
                          axis.text.y=element_text(size=10))  # Y axis text

plotinferCNV1 <- ggplot(melted.inferCNV.exp.matrix, aes(exp_value)) + 
  geom_density(alpha=0.2) + 
  labs(fill="Stats", 
       title = "inferCNV",
       x="Normalized Expression",
       y="Density") + 
  xlim(-3, +3) +
  theme_classic() + theme(plot.title=element_text(size=20, 
                                                  face="bold", 
                                                  hjust=0.5,
                                                  lineheight=1.2),  # title
                          axis.title.x=element_text(size=15, face = "bold"),  # X axis title
                          axis.title.y=element_text(size=15, face = "bold"),  # Y axis title
                          axis.text.x=element_text(size=10, 
                                                   angle = 30,
                                                   vjust=.5),  # X axis text
                          axis.text.y=element_text(size=10))  # Y axis text

# save patchwork plot
pdf("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/exp_limitedX_axis_density.pdf", width = 15, height = 20)
print(plotWGS / plothb1 / plotCONICS1 / plotinferCNV1)
dev.off()

pdf("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/exp_density_scRNAtools.pdf", width = 30, height = 10)
print(plothb1 | plotCONICS1 | plotinferCNV1)
dev.off()

#######################################################################################                                                                                                                     
##                     MERGE WGS WITH THE INDIVIDUAL TOOL - HoneyBADGER                                                                                      
#######################################################################################
# merge the individual melted dfs with the logR dfs
hb.final.df <- merge(melted.tumor.hb.exp.matrix, hb.logR, by = "hgnc_symbol")

# add gene coordinates
hb.GRange.df <- as.data.frame(hb.GRange)
hb.GRange.df$hgnc_symbol <- rownames(hb.GRange.df)
hb.final.df <- merge(hb.final.df, hb.GRange.df, by.x = "hgnc_symbol")
colnames(hb.final.df) <- c("hgnc_symbol", "barcode", "exp_value", "mean_logR",
                           "gene_length", "chr", "start", "end", "width", "strand", "length")

# remove unnecessary columns
hb.final.df[,c(9:11)] <- NULL

# add coordinates column
hb.final.df$coordinates <- paste(hb.final.df$chr, hb.final.df$start, hb.final.df$end, sep = ":")

# check chromsome order
chrOrder <- paste0("chr", c(1:22))
hb.final.df$chr <- factor(hb.final.df$chr, levels = chrOrder)

#######################################################################################                                                                                                                     
##                     MERGE WGS WITH THE INDIVIDUAL TOOL - CONICSmat                                                                                     
#######################################################################################
# merge the individual melted dfs with the logR dfs
CONICS.final.df <- merge(melted.tumor.CONICS.exp.matrix, CONICS.logR, by = "hgnc_symbol")

# add gene coordinates
CONICS.GRange.df <- as.data.frame(CONICS.GRange)
CONICS.GRange.df$hgnc_symbol <- rownames(CONICS.GRange.df)
CONICS.final.df <- merge(CONICS.final.df, CONICS.GRange.df, by.x = "hgnc_symbol")
colnames(CONICS.final.df) <- c("hgnc_symbol", "barcode", "exp_value", "mean_logR",
                               "gene_length", "chr", "start", "end", "width", "strand", "length")

# remove unnecessary columns
CONICS.final.df[,c(9:11)] <- NULL

# add coordinates column
CONICS.final.df$coordinates <- paste(CONICS.final.df$chr, CONICS.final.df$start, CONICS.final.df$end, sep = ":")

# check chromsome order
chrOrder <- paste0("chr", c(1:22))
CONICS.final.df$chr <- factor(CONICS.final.df$chr, levels = chrOrder)

#######################################################################################                                                                                                                     
##                     MERGE WGS WITH THE INDIVIDUAL TOOL - inferCNV                                                                                      
#######################################################################################
# merge the individual melted dfs with the logR dfs
inferCNV.final.df <- merge(melted.inferCNV.exp.matrix, inferCNV.logR, by = "hgnc_symbol")
inferCNV.final.df <- inferCNV.final.df[,c(1:3, 7:8)]

# add gene coordinates
inferCNV.GRange.df <- as.data.frame(inferCNV.GRange)
inferCNV.GRange.df$hgnc_symbol <- rownames(inferCNV.GRange.df)
inferCNV.final.df <- merge(inferCNV.final.df, inferCNV.GRange.df, by.x = "hgnc_symbol")
colnames(inferCNV.final.df) <- c("hgnc_symbol", "barcode", "exp_value", "mean_logR",
                                 "gene_length", "chr", "start", "end", "width", "strand", "length")

# remove unnecessary columns
inferCNV.final.df[,c(9:11)] <- NULL

# add coordinates column
inferCNV.final.df$coordinates <- paste(inferCNV.final.df$chr, inferCNV.final.df$start, inferCNV.final.df$end, sep = ":")

# check chromsome order
chrOrder <- paste0("chr", c(1:22))
inferCNV.final.df$chr <- factor(inferCNV.final.df$chr, levels = chrOrder)

#######################################################################################                                                                                                                     
##            WRITE A DATAFRAME CONTAINING ALL THE INFORMATION PER METHOD                                                                                      
#######################################################################################

## HoneyBADGER
# calculate z-score
hb.z.score <- (hb.final.df$exp_value - mean(hb.final.df$exp_value))/sd(hb.final.df$exp_value)
hb.final.df$Z_score <- hb.z.score

# write the data.frame with all information to a file
write.table(hb.final.df, "/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/hb_complete_df.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# calculate the max, min, mean number of genes per chromosome
hb.gene.per.chr <- table(unique(hb.final.df[c(1,6)]))
hb.gene.per.chr <- colSums(hb.gene.per.chr)
hb.mean.gene.per.chr <- mean(hb.gene.per.chr)
hb.min.gene.per.chr <- min(hb.gene.per.chr)
hb.max.gene.per.chr <- max(hb.gene.per.chr)

## CONICSmat
# calculate z-score
CONICS.z.score <- (CONICS.final.df$exp_value - mean(CONICS.final.df$exp_value))/sd(CONICS.final.df$exp_value)
CONICS.final.df$Z_score <- CONICS.z.score

# write the data.frame with all information to file
write.table(CONICS.final.df, "/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/CONICS_complete_df.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# calculate the max, min, mean number of genes per chromosome
CONICS.gene.per.chr <- table(unique(CONICS.final.df[c(1,6)]))
CONICS.gene.per.chr <- colSums(CONICS.gene.per.chr)
CONICS.mean.gene.per.chr <- mean(CONICS.gene.per.chr)
CONICS.min.gene.per.chr <- min(CONICS.gene.per.chr)
CONICS.max.gene.per.chr <- max(CONICS.gene.per.chr)

## inferCNV
# calculate z-score
inferCNV.z.score <- (inferCNV.final.df$exp_value - mean(inferCNV.final.df$exp_value))/sd(inferCNV.final.df$exp_value)
inferCNV.final.df$Z_score <- inferCNV.z.score

# write the data.frame with all information to file
write.table(inferCNV.final.df, "/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/inferCNV_complete_df.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# calculate the max, min, mean number of genes per chromosome
inferCNV.gene.per.chr <- table(unique(inferCNV.final.df[c(1,6)]))
inferCNV.gene.per.chr <- colSums(inferCNV.gene.per.chr)
inferCNV.mean.gene.per.chr <- mean(inferCNV.gene.per.chr)
inferCNV.min.gene.per.chr <- min(inferCNV.gene.per.chr)
inferCNV.max.gene.per.chr <- max(inferCNV.gene.per.chr)

#######################################################################################                                                                                                                     
##            Subset the individual dataframes to genes which are altered                                                                                      
#######################################################################################

# make a dataframe with the counts per chromosome
chr.gene.df <- as.data.frame(inferCNV.gene.per.chr)
colnames(chr.gene.df) <- "inferCNV"

# add info for HoneyBadger and CONICSmat
chr.gene.df$HoneyBADGER <- hb.gene.per.chr
chr.gene.df$CONICSmat <- CONICS.gene.per.chr
chr.gene.df$chr <- rownames(chr.gene.df)

# melt dataframe to plot
melted.chr.gene.df <- melt(chr.gene.df)

# check chromosome order
chrOrder <- paste0("chr", c(1:22))
melted.chr.gene.df$chr <- factor(melted.chr.gene.df$chr, levels = chrOrder)
melted.chr.gene.df$variable <- factor(melted.chr.gene.df$variable, levels = c("CONICSmat", "inferCNV", "HoneyBADGER"))

# plot number of chromosomes per gene
ggplot(melted.chr.gene.df, aes(x=chr, y=value, fill = variable)) + 
  geom_bar(stat="identity", position = position_dodge2(preserve = c("single"))) + 
  labs(title="#genes per chr", 
       y = "#genes", x = "Chromosomes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust=0.6)) +
  scale_fill_manual(values = viridis(3), name = "Algorithm") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "top")
ggsave("/Users/mjprzy/Documents/Master_thesis/Data/correlation_test/numberGenes_Methods.pdf", dpi = 600, width = 8.5, height = 7.5)

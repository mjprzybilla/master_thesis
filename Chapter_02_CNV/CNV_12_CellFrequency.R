#############################################################################################################################
##                                                                                                                      
##  DIRECTLY ASSESS THE FREQUENCY OF ALTERATIONS FOR THE ESSENTIAL CHROMOSOMES KNOWN FROM WGS
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
list.of.packages <- c("tidyverse", "Seurat", "biomaRt", "RColorBrewer", "GenomicRanges", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# define the input directory containing the seurat objects
w.dir <- "/labs/ccurtis2/mjprzy/correlation_test/"

# get sample name
sample.tmp <- "P2C2R2T2"
message(sample.tmp)

###########################################################################
#                   READ IN DATA OF INTEREST FROM inferCNV
###########################################################################
## Read in infercnv data for 6077 C5 R2 T2

# read in the bayesNet output for these samples
BayesNet.matrix <- as.matrix(read.table(paste0("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/Sequencing8_C5T2R2_cr3_WT/results/BayesNetOutput.HMMi6.rand_trees.hmm_mode-subclusters/infercnv.NormalProbabilities.PostFiltering.observations.txt"), header = T))

# get dimensions
print(dim(BayesNet.matrix))

# adapt barcodes
colnames(BayesNet.matrix) <- str_split_fixed(colnames(BayesNet.matrix), "_pos", 2)[,1]

# list of genes in matrix
genes.of.interest <- data.frame(hgnc_symbol = rownames(BayesNet.matrix))

# use biomart to get hgnc_symbols
mart = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl", mirror = "uswest")
ann <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
             "hgnc_symbol", as.character(genes.of.interest$hgnc_symbol), mart)

ann$chromosome_name <- paste0("chr", ann$chromosome_name)

# remove duplicated gene names
ann <- distinct(ann, ensembl_gene_id, .keep_all = TRUE)

###########################################################################
#                   READ IN DATA OF INTEREST FROM CONICS
###########################################################################

# read in CONICSmat gene expression matrix
CONICS.matrix <- readRDS("/labs/ccurtis2/mjprzy/correlation_test/P2C2R2T2/CONICSmat.binMatrix.rds")

# get dimensions
print(dim(CONICS.matrix))

# subset to WGS dimensions
CONICS.matrix <- CONICS.matrix[, colnames(CONICS.matrix) %in% c("del_3p", "amp_3q", "del_4p", "del_4q", "del_9p", "amp_9q", "del_13p", "del_13q", "amp_20q")]

# only consider complete cases
CONICS.matrix <- CONICS.matrix[complete.cases(CONICS.matrix),]
CONICS.matrix <- CONICS.matrix[-grep("WT", rownames(CONICS.matrix)), ]
rownames(CONICS.matrix) <- str_split_fixed(rownames(CONICS.matrix), "-", 2)[,1]

###########################################################################
#                   READ IN DATA OF INTEREST FROM SCDNA
###########################################################################

# read in ploidy matrix from scDNA (construct_ploidyMatrix_scDNA.R)
scDNA.matrix <- as.matrix(read.table("/labs/ccurtis2/mjprzy/correlation_test/cnv.cells.mat.txt", header = T, stringsAsFactors = F))

# get dimensions
print(dim(scDNA.matrix))

# get the segments and make GRange from it 
seg.coordinates <- data.frame(chr = paste0("chr", str_split_fixed(scDNA.matrix[,1], ":",2)[,1]), 
                              start = str_split_fixed(str_split_fixed(scDNA.matrix[,1], ":",2)[,2], "-", 2)[,1], 
                              end = str_split_fixed(str_split_fixed(scDNA.matrix[,1], ":",2)[,2], "-", 2)[,2])
seg.coordinates.GRange <- makeGRangesFromDataFrame(seg.coordinates, keep.extra.columns = T)

# remove id column
rownames(scDNA.matrix) <- paste(seg.coordinates$chr, seg.coordinates$start, seg.coordinates$end, sep = "-")
scDNA.matrix <- scDNA.matrix[,-1]

###########################################################################
#                    PROCESS THE inferCNV OUTPUT
###########################################################################
# read in the chromosome arm positions
chr.arm.pos <- read.table("/labs/ccurtis2/mjprzy/correlation_test/chromosome_arm_positions_GRCH38.txt", header = T)
chr.arm.pos$Chrom <- paste0("chr", chr.arm.pos$Chrom)
chr.arm.pos.GRange <- makeGRangesFromDataFrame(chr.arm.pos, keep.extra.columns = T)

# melt df/matrix with exp values
melted.BayesNet.matrix <- reshape2::melt(BayesNet.matrix)
colnames(melted.BayesNet.matrix) <- c("hgnc_symbol", "barcode", "post_prob")
melted.BayesNet.matrix <- merge(melted.BayesNet.matrix, ann, by = "hgnc_symbol")
colnames(melted.BayesNet.matrix) <- c("hgnc_symbol", "barcode", "post_prob", "chr", "start", "end", "ensembl_gene_id")

# make GRange object for infercnv genes
inferCNV.genes <- makeGRangesFromDataFrame(melted.BayesNet.matrix, keep.extra.columns = T)

# then merge each dataframe information by overlap
# find the overlapping regions of the WGS segments with the genes
inferCNV.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, inferCNV.genes)

# subset it to altered chromosomes
subset.inferCNV.GR.merge <- inferCNV.GR.merge[inferCNV.GR.merge$Idf %in% c("3p", "3q", "4p", "4q", "9p", "9q", "13p", "13q", "20q"),]

# get the essential columns and make a new dataframe
chr.arm.inferCNV <- data.frame("chr_arm" = subset.inferCNV.GR.merge$Idf,
                            "Cell_barcode" = subset.inferCNV.GR.merge$barcode,
                            "post_prob" = subset.inferCNV.GR.merge$post_prob,
                            stringsAsFactors = F)

# make the post_prob binary
chr.arm.inferCNV$post_prob[chr.arm.inferCNV$post_prob >= 0.5] <- 1
chr.arm.inferCNV$post_prob[chr.arm.inferCNV$post_prob < 0.5] <- 0

# convert each row
chr.arm.inferCNV$chr_arm <- as.character(chr.arm.inferCNV$chr_arm)
chr.arm.inferCNV$Cell_barcode <- as.character(chr.arm.inferCNV$Cell_barcode)

sd <- sd(chr.arm.inferCNV$post_prob)
# reshape the melted dataframe into an matrix 
inferCNV.bin.matrix <- reshape2::acast(chr.arm.inferCNV, Cell_barcode ~ chr_arm, value.var = "post_prob", fun.aggregate = mean)

# make binary again if more than 50% were assigned as a alteration before
inferCNV.bin.matrix[which(inferCNV.bin.matrix >= 0.5,arr.ind = T)] <- 1
inferCNV.bin.matrix[which(inferCNV.bin.matrix < 0.5,arr.ind = T)] <- 0

# reorder columns
inferCNV.bin.matrix <- inferCNV.bin.matrix[,c(3:ncol(inferCNV.bin.matrix),1,2)]
colnames(inferCNV.bin.matrix) <- colnames(CONICS.matrix)

###########################################################################
#                    PROCESS THE scCNV OUTPUT
###########################################################################
# melt df/matrix with exp values
melted.scDNA.matrix <- reshape2::melt(scDNA.matrix)
colnames(melted.scDNA.matrix) <- c("seg_coordinates", "cell_id", "copy_number")
melted.scDNA.matrix <- melted.scDNA.matrix %>% separate(seg_coordinates, c("chr", "start", "end"), "-")

# make GRange object for infercnv genes
scDNA.genes <- makeGRangesFromDataFrame(melted.scDNA.matrix, keep.extra.columns = T)

# then merge each dataframe information by overlap
# find the overlapping regions of the WGS segments with the genes
scDNA.GR.merge <- mergeByOverlaps(chr.arm.pos.GRange, scDNA.genes)

# subset it to altered chromosomes
subset.scDNA.GR.merge <- scDNA.GR.merge[scDNA.GR.merge$Idf %in% c("3p", "3q", "4p", "4q", "9p", "9q", "13p", "13q", "20q"),]

# get the essential columns and make a new dataframe
chr.arm.scDNA <- data.frame("chr_arm" = scDNA.GR.merge$Idf,
                               "cell_id" = scDNA.GR.merge$cell_id,
                               "copy_number" = scDNA.GR.merge$copy_number,
                               stringsAsFactors = F)

# convert each row
chr.arm.scDNA$chr_arm <- as.character(chr.arm.scDNA$chr_arm)
chr.arm.scDNA$cell_id <- as.character(chr.arm.scDNA$cell_id)
gsub("  ", "", chr.arm.scDNA$copy_number)

test <- gsub("  ", "", as.character(chr.arm.scDNA$copy_number))
test <- gsub(" ", "", test)
test <- as.numeric(test)
chr.arm.scDNA$copy_number <- as.numeric(test)

# make the copy number binary
chr.arm.scDNA$copy_number[chr.arm.scDNA$copy_number < 2 | chr.arm.scDNA$copy_number > 2 ] <- 1
chr.arm.scDNA$copy_number[chr.arm.scDNA$copy_number == 2] <- 0

# reshape the melted dataframe into an matrix 
scDNA.arm.matrix <- reshape2::acast(chr.arm.scDNA, cell_id ~ chr_arm, value.var = "copy_number", fun.aggregate = mean)

# reorder columns
scDNA.arm.matrix<- scDNA.arm.matrix[, colnames(scDNA.arm.matrix) %in% c("3p", "3q", "4p", "4q", "9p", "9q", "13q", "20q")]
scDNA.arm.matrix <- scDNA.arm.matrix[,c(4:ncol(scDNA.arm.matrix), 1:3)]
colnames(scDNA.arm.matrix) <- colnames(CONICS.matrix)

# make the copy number binary
scDNA.arm.matrix[which(scDNA.arm.matrix < 0.5, arr.ind = T)] <- 0
scDNA.arm.matrix[which(scDNA.arm.matrix >= 0.5, arr.ind = T)] <- 1

###########################################################################
#                 CALCULATE THE FRACTIONS AND PLOT THEM
###########################################################################
# calculate percentage of cell carrying each alteration
CONICS.cell.percentage <- round(colSums(CONICS.matrix)/nrow(CONICS.matrix), 3)

# calculate percentage of cells carrying each alteration
inferCNV.cell.percentage <- round(colSums(inferCNV.bin.matrix)/nrow(inferCNV.bin.matrix), 3)

# calculate percentage of cells carrying each alteration
scDNA.cell.percentage <- round(colSums(scDNA.arm.matrix)/nrow(scDNA.arm.matrix), 3)

# combine the frequencies into a dataframe
cell.perc.df <- rbind(data.frame(chr_arm = names(CONICS.cell.percentage), cell_percentage = CONICS.cell.percentage*100, tool = "CONICSmat"), 
                      data.frame(chr_arm = names(inferCNV.cell.percentage), cell_percentage = inferCNV.cell.percentage*100, tool = "inferCNV"), 
                      data.frame(chr_arm = names(scDNA.cell.percentage), cell_percentage = scDNA.cell.percentage*100, tool = "scCNV"))

# determine order
arm.order <- colnames(CONICS.matrix)
cell.perc.df$chr_arm <- factor(cell.perc.df$chr_arm, levels = arm.order)

# plot the distinct representations
ggplot(cell.perc.df, aes(x=chr_arm, y=cell_percentage, fill = tool)) + 
  geom_bar(stat="identity", position = position_dodge2(preserve = "single")) + 
  labs(y = "Fraction of cells [%]", x = "Copy Number Alterations") +
  theme_classic() + 
  scale_fill_manual(values = viridis(3), name = "Method") +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold" ),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold" ),
        axis.title = element_text(colour = "black", size = 16, face = "bold" ),
        plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 12, face = "bold",),
        legend.text = element_text(colour="black", size=10, face="bold"),
        legend.position = "top")
ggsave(paste0(w.dir, "subclonal_comparison.pdf"), height = 5, width = 15, dpi = 600)



#############################################################################################################################
##                                                                                                                      
##  Analysis of individual subclones in scRNA-seq with ECB information from Sequencing 18
##                                                                                                                      
##  Date: 23 July 2020                                                                                                                    
##  
##  Author: Moritz Przybilla and Kasper Karlsson                                                                                                                 
##                                                                                                                      
#############################################################################################################################

# clear workspace
rm(list = ls())
set.seed(14) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("tidyverse", "patchwork", "Seurat", "Matrix", "biomaRt", "scater", "DoubletFinder", "viridis",
                      "GO.db", "HTSanalyzeR2", "org.Hs.eg.db", "KEGGREST", "igraph", "tidyr", "stats", "reshape2", "ggplot2",
                      "forcats","heatmap3", "gplots", "ggpubr", "dplyr", "ComplexHeatmap", "ggplotify", "karyoploteR", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                      "regioneR", "RIPSeeker", "Repitools")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

###########################################################################
#                               FUNCTIONS
###########################################################################

# set up functions which are used 
`%notin%` <- Negate(`%in%`)

### FUNCTION TO GET GSEA ENRICHMENT FOR SEURAT OUTPUT FROM FIND MARKERS

geneSetAnalysis <- function(dfile) {
  dfile_subs <- subset(dfile,p_val <= 0.05) ### SUBSET DEGS FILE TO PVAL <= 0.05
  phenotype <- as.vector(dfile_subs$avg_logFC) ### USE LOG FOLD CHANGE AS PHENOTYPE VECTOR
  names(phenotype) <- rownames(dfile_subs)
  
  ## specify the gene sets type you want to analyze
  
  # HALLMARK
  MSig_H <- MSigDBGeneSets(species = "Hs", collection = "H", subcategory = NULL) # Hallmarks!
  ListGSC <- list(MSig_H=MSig_H)
  
  
  ## iniate a *GSCA* object
  gsca <- GSCA(listOfGeneSetCollections=ListGSC, 
               geneList=phenotype)
  
  ## preprocess
  gsca1 <- preprocess(gsca, species="Hs", initialIDs="SYMBOL",
                      keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                      orderAbsValue=FALSE)
  
  
  ## analysis
  if (requireNamespace("doParallel", quietly=TRUE)) {
    doParallel::registerDoParallel(cores=4)
  }  ## support parallel calculation using multiple cores
  
  
  gsca2 <- analyze(gsca1, 
                   para=list(pValueCutoff=0.05, pAdjustMethod="BH",
                             nPermutations=10000, minGeneSetSize=1,
                             exponent=1), 
                   doGSOA = FALSE)
  return(getResult(gsca2)$GSEA.results$MSig_H)
  message("Successfully ran GSEA.")
}

####### FUNCTION TO GET PROPORTION OF GENES PER CHROMOSOME ARE THAT ARE DIFFERENTIALLY EXPRESSED
getAveAltUp <- function(rg,propup){
  
  ### READ IN DEGS FROM SEURAT FIND.MARKERS FILE
  degs <- read.table(paste0(sample.tmp, "/DEG/markers_",rg,"_vs_2c.txt"), sep="\t", header = T, stringsAsFactors = T)
  
  ### MERGE WITH GENCODE GENE LOCATION ANNOTATION
  degs$genes <- rownames(degs)
  degsComb <- merge(gencode,degs,by="genes")
  
  ### ADD CHR ARM ANNOATION
  degsComb$chrarm <- paste0(degsComb$chr,degsComb$arm)
  
  ### ADD DIRECTION (UP OR DOWN REGULATED)
  degsComb$direction = ifelse(degsComb$avg_logFC>0,1,-1)
  
  ### GET NR OF GENES PER CHR ARM THAT ARE UP OR DOWN REGULATED
  degsCombTable <- table(degsComb$chrarm,degsComb$direction)
  degsCombDF <- as.data.frame.matrix(degsCombTable)
  colnames(degsCombDF) <- c("downR","upR")
  degsCombDF$chrarm <- rownames(degsCombDF)
  
  ### MERGE WITH THE PREVIOUSLY CREATED DF TO ADD CHR WITH ZERO COUNTS
  degsCombDF_all <- merge(degsCombDF,chrdf0,by="chrarm",all.y=TRUE)
  degsCombDF_all$value <- NULL
  degsCombDF_all[is.na(degsCombDF_all)] = 0
  
  ### ADD NR GENES EXPRESSED PER CHR ARM
  degsCombDF_all_nrGenes <- merge(degsCombDF_all,genes_per_arm,by="chrarm")
  
  ### CALCULATE PROPORTION THAT IS DIFFERENTIALLY EXPRESSED
  degsCombDF_all_nrGenes$upProp <- degsCombDF_all_nrGenes$upR/degsCombDF_all_nrGenes$nrGenes
  propup <- append(propup,degsCombDF_all_nrGenes$upProp)
  return (propup)
}

### SAME AS ABOVE BUT FOR DOWNREGULATED GENES
getAveAltDown <- function(rg,propdown){
  
  ### READ IN DEGS FROM SEURAT FIND.MARKERS FILE
  degs <- read.table(paste0(sample.tmp, "/DEG/markers_",rg,"_vs_2c.txt"), sep="\t", header = T, stringsAsFactors = T)
  
  ### MERGE WITH GENCODE GENE LOCATION ANNOTATION
  degs$genes <- rownames(degs)
  degsComb <- merge(gencode,degs,by="genes")
  
  ### ADD CHR ARM ANNOATION
  degsComb$chrarm <- paste0(degsComb$chr,degsComb$arm)
  
  ### ADD DIRECTION (UP OR DOWN REGULATED)
  degsComb$direction = ifelse(degsComb$avg_logFC>0,1,-1)
  
  
  ### GET NR OF GENES PER CHR ARM THAT ARE UP OR DOWN REGULATED
  degsCombTable <- table(degsComb$chrarm,degsComb$direction)
  degsCombDF <- as.data.frame.matrix(degsCombTable)
  colnames(degsCombDF) <- c("downR","upR")
  degsCombDF$chrarm <- rownames(degsCombDF)
  
  ### MERGE WITH THE PREVIOUSLY CREATED DF TO ADD CHR WITH ZERO COUNTS
  degsCombDF_all <- merge(degsCombDF,chrdf0,by="chrarm",all.y=TRUE)
  degsCombDF_all$value <- NULL
  degsCombDF_all[is.na(degsCombDF_all)] = 0
  
  ### ADD NR GENES EXPRESSED PER CHR ARM
  degsCombDF_all_nrGenes <- merge(degsCombDF_all,genes_per_arm,by="chrarm")
  
  ### CALCULATE PROPORTION THAT IS DIFFERENTIALLY EXPRESSED
  degsCombDF_all_nrGenes$downProp <- degsCombDF_all_nrGenes$downR/degsCombDF_all_nrGenes$nrGenes
  propdown <- append(propdown,degsCombDF_all_nrGenes$downProp)
  return (propdown)
}

### FUNCTION TO PLOT PROPORTION DEGS PER CHR THAT ARE UP OR DOWN REGEULATED
plotDEGs <- function(rg){
  
  ### READ IN DEGS FROM SEURAT FIND.MARKERS FILE
  degs <- read.table(paste0(sample.tmp, "/DEG/markers_",rg,"_vs_2c.txt"), sep="\t", header = T, stringsAsFactors = T)
  
  
  ### MERGE WITH GENCODE GENE LOCATION ANNOTATION
  degs$genes <- rownames(degs)
  degsComb <- merge(gencode,degs,by="genes")
  
  ### ADD CHR ARM ANNOATION
  degsComb$chrarm <- paste0(degsComb$chr,degsComb$arm)
  
  ### ADD DIRECTION (UP OR DOWN REGULATED)
  degsComb$direction = ifelse(degsComb$avg_logFC>0,1,-1)
  
  ### GET NR OF GENES PER CHR ARM THAT ARE UP OR DOWN REGULATED
  degsCombTable <- table(degsComb$chrarm,degsComb$direction)
  degsCombDF <- as.data.frame.matrix(degsCombTable)
  colnames(degsCombDF) <- c("downR","upR")
  degsCombDF$chrarm <- rownames(degsCombDF)
  
  ### MERGE WITH THE PREVIOUSLY CREATED DF TO ADD CHR WITH ZERO COUNTS
  degsCombDF_all <- merge(degsCombDF,chrdf0,by="chrarm",all.y=TRUE)
  degsCombDF_all$value <- NULL
  degsCombDF_all[is.na(degsCombDF_all)] = 0
  
  ### ADD NR GENES EXPRESSED PER CHR ARM
  degsCombDF_all_nrGenes <- merge(degsCombDF_all,genes_per_arm,by="chrarm")
  
  ### CALCULATE PROPORTION THAT IS DIFFERENTIALLY EXPRESSED
  degsCombDF_all_nrGenes$downProp <- degsCombDF_all_nrGenes$downR/degsCombDF_all_nrGenes$nrGenes
  degsCombDF_all_nrGenes$upProp <- degsCombDF_all_nrGenes$upR/degsCombDF_all_nrGenes$nrGenes
  
  ### EXTRACT CHROMOSOME ARMS OF INTEREST
  degsCombDF_all_nrGenes_subset <- subset(degsCombDF_all_nrGenes,chrarm %in% keeps)
  
  ### ORDER CHROMOSOME ARMS AS IN THE LIST "KEEPS"
  degsCombDF_all_nrGenes_subset$chrarm <- factor(degsCombDF_all_nrGenes_subset$chrarm,levels = keeps)
  
  ### GGPLOT TO PLOT UP REGUALATED DEG PROPORTION  
  up <- ggplot(data=degsCombDF_all_nrGenes_subset, aes(x=chrarm, y=upProp)) +
    geom_bar(stat="identity",width = 0.8,fill="darkred") +
    theme(axis.line = element_line(colour = "black"),legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.text.y = element_text(size=10, face = "bold", color = "black"),
          axis.text.x = element_text(size=10, face = "bold", color = "black", hjust = 0.5,vjust=-2))+
    scale_y_continuous(limits = c(0, 0.25),labels = scales::percent)+
    geom_hline(yintercept=median(propup), linetype="dashed", color = "black", size=0.5)
  
  ### GGPLOT TO PLOT DOWN REGUALATED DEG PROPORTION  
  down <- ggplot(data=degsCombDF_all_nrGenes_subset, aes(x=chrarm, y=downProp)) +
    geom_bar(stat="identity",width = 0.8,fill="darkblue") +
    theme(axis.line = element_line(colour = "black"),legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          panel.background = element_rect(fill = 'white'),
          axis.text.y = element_text(size=10, face = "bold", color = "black"))+
    scale_x_discrete(position = "top") +
    scale_y_reverse(limits = c(0.25, 0),labels = scales::percent)+
    geom_hline(yintercept=median(propdown), linetype="dashed", color = "black", size=0.5)
  
  ### PRINT COMBINED PLOT TO FILE
  
  
  pdf(paste0(o.dir, "/", sample.tmp, "/DEG/markers_",rg,"_vs_2c.pdf"),width=9.5,height=1.5,pointsize=0.1)
  print (ggarrange(up,down, ncol = 1, nrow = 2))
  dev.off()
}

###########################################################################
#                       READ IN THE DATA OF INTEREST
###########################################################################

# define the cellranger input directory of the sample of interest
c.dir <- "/labs/ccurtis2/mjprzy/infercnv_gastric/Sequencing18_D3_APA4_Early/outs/"

# define output directory
o.dir <- "/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze"

# create o.dir
dir.create(o.dir)
setwd(o.dir)

# get sample.ids
sample.tmp <- unique(str_split_fixed(list.files(c.dir, full.names = T), "/", 8)[,6])

# which sample? 
message(sample.tmp)
dir.create(paste0(o.dir, "/", sample.tmp, "/DEG"))

# read seurat objects
seurat.obj <- readRDS(file = paste0(o.dir, "/" , sample.tmp, "/", sample.tmp, "_seurat_obj.rds"))

# change the cluster identities from RTAM to SCT
Idents(seurat.obj) <- "SCT_snn_res.0.5"

# read in Genecode data
gencode <- read.csv("/labs/ccurtis2/mjprzy/scRNA_analysis/KarlssonEtal_Figure3/gencode_v34_simple_Havana.txt", sep="\t", header=FALSE)
colnames(gencode) <- c("chr","start","end","arm","genes")
gencode$chrarm <- paste0(gencode$chr,gencode$arm)

###########################################################################
#         GET THE BARCODES FOR 0A AND 0B AND ADAPT SEURAT OBJECT
###########################################################################

# read in the copy number clone file from the copy number assignement in inferCNV
cnv.clone.files <- list.files("/labs/ccurtis2/mjprzy/infercnv_gastric/freeze/arm_level", pattern = "_cnv_clones.txt", full.names = T)
cnv.clone.files <- cnv.clone.files[grep("4230_AP_R1_T2", cnv.clone.files)]

# read file
clone.file <- read.table(cnv.clone.files[1], header = T, stringsAsFactors = F)
clone.file$Cell_barcodes <- gsub("\\.", "-", clone.file$Cell_barcodes)

# get the metadata from the seurat.obj
seurat.metadata <- seurat.obj@meta.data

# subset to barcode 0
clone0.barcodes <- seurat.metadata[seurat.metadata$ECB_RG == "0", "Cell_Barcode"]
clone2.barcodes <- seurat.metadata[seurat.metadata$ECB_RG == "2", "Cell_Barcode"]

# subset the clone file to subclone 0
clone0.file <- clone.file[clone.file$Cell_barcodes %in% clone0.barcodes, ]
clone0a.barcodes <- clone0.file[grep("9q+", clone0.file$cnv_id), "Cell_barcodes"]
clone0b.barcodes <- clone0.file[-grep("9q+", clone0.file$cnv_id), "Cell_barcodes"]

# subset the clone file to subclone 2
clone2.file <- clone.file[clone.file$Cell_barcodes %in% clone2.barcodes, ]
clone2.file[clone2.file$clone_id == "clone39", "cnv_id"] <- "3p-:16p+"
clone2c.barcodes <- clone2.file[clone2.file$clone_id == "clone22", "Cell_barcodes"]
clone2b.barcodes <- clone2.file[clone2.file$clone_id == "clone19", "Cell_barcodes"]
clone2a.barcodes <- clone2.file[clone2.file$clone_id == "clone18", "Cell_barcodes"]

# write barcodes to file
write.table(data.frame(clone2c.barcodes), paste0(sample.tmp, "/DEG/clone2c_barcodes.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# change the identity of the seurat object to subclones
Idents(seurat.obj) <- "ECB_RG"

# change the identity of the subset of cells belonging to 0a to 0a
seurat.obj <- SetIdent(seurat.obj, cells = clone0a.barcodes, value = "0a")
seurat.obj <- SetIdent(seurat.obj, cells = clone0b.barcodes, value = "0b")

# change the identity of the subset of cells belonging to 0a to 0a
seurat.obj <- SetIdent(seurat.obj, cells = clone2c.barcodes, value = "2c")
seurat.obj <- SetIdent(seurat.obj, cells = clone2b.barcodes, value = "2b")
seurat.obj <- SetIdent(seurat.obj, cells = clone2a.barcodes, value = "2a")

# get the valid ECB clones
ECB.clones <- table(Idents(seurat.obj)) > 20
valid.clones <- names(ECB.clones[ECB.clones == T])

# manually add 0a
valid.clones <- factor(valid.clones, levels = c("0a", "0b", "2a", "2b", "2c", 0:10000))
valid.clones <- valid.clones[order(valid.clones)]

# remove clone 1
valid.clones <- valid.clones[-c(5,6,8)]

###########################################################################
#             PERFORM DIFFERENTIAL GENE EXPRESSION ANALYSIS
###########################################################################

# iterate over each clone and FindMarkers
for (subclone in valid.clones) {
  print(subclone)
  
  # run DGE analysis
  subclone.markers <- FindMarkers(seurat.obj, ident.1= subclone, ident.2="2c", logfc.threshold = 0.25)
  
  # write all the markers to a file
  write.table(subclone.markers, paste0(o.dir, "/" , sample.tmp, "/DEG/markers_subclone", subclone,"_vs_2c.txt") , quote=FALSE, sep="\t", col.names = T, row.names = T)
}

###########################################################################
#           GET THE NUMBER OF GENES WHICH IS EXPRESSED PER ARM
###########################################################################
# set identities to the same
seurat.obj$subclones <- Idents(seurat.obj)
Idents(seurat.obj) <- sample.tmp

# calculate average gene expression per clone id for Aziz
SCT.pseudobulk <- AverageExpression(seurat.obj, assay = "SCT")
SCT.pseudobulk <- SCT.pseudobulk$SCT

### ORDER IN DECREASING ORDER
SCT.pseudobulk <- SCT.pseudobulk[order(SCT.pseudobulk$Sequencing18_D3_APA4_Early, decreasing = TRUE),, drop = FALSE]
SCT.pseudobulk$genes <- rownames(SCT.pseudobulk)

### MERGE WITH GENCODE TO GET CHROMOSOME ARM LOCATION PER GENE
SCT.pseudobulk.comb <- merge(gencode,SCT.pseudobulk,by="genes")
SCT.pseudobulk.comb$chrarm <- paste0(SCT.pseudobulk.comb$chr,SCT.pseudobulk.comb$arm)

### PRINT TO FILE 
write.table(SCT.pseudobulk,paste0(sample.tmp, "/DEG/SCT_average_exp_all_genes.txt"), quote=FALSE, sep="\t")
write.table(table(SCT.pseudobulk.comb$chrarm), paste0(sample.tmp, "/DEG/SCT_exp_genes_per_chrarm.txt"), quote=FALSE, sep="\t")

###########################################################################
#           GET THE NUMBER OF GENES WHICH IS EXPRESSED PER ARM
###########################################################################

# read in files of interest
marker.files <- list.files(paste0("/labs/ccurtis2/mjprzy/scRNA_analysis/hashECB_data_freeze/", sample.tmp, "/DEG"), pattern="markers_", full.names=T)

# get subclone id
subclone.id <- str_split_fixed(basename(marker.files), "_", 3)[,2]

# READ IN ALL FILES
marker.files <- lapply(marker.files, read.table, sep="\t", header=TRUE, stringsAsFactors = F)

# GET GENE SET ENRICHMENT FOR EACH DEG FILE
GSEA.results <- lapply(marker.files, try(geneSetAnalysis))

### PRINT TO FILE
for (i in 1:length(GSEA.results)){
  
  # add the hallmark column to the file
  results <- GSEA.results[[i]]
  results$hallmark <- rownames(results)
  
  # add a subclone column
  results$subclone <- subclone.id[i]
  
  # exchange the files
  GSEA.results[[i]] <- results
  
  write.table(results, paste0(sample.tmp, "/DEG/", subclone.id[i], "_Hallmark_GSEA_vs_2c.txt"),sep="\t", quote=FALSE, row.names = TRUE, col.names=TRUE)
}

# combine all the files together
all <- bind_rows(GSEA.results)
all$Leading.Edge <- NULL

# replace score according to down or up-regulation
all[all$Observed.score < 0, "Observed.score"] <- -1
all[all$Observed.score > 0, "Observed.score"] <- 1

# replace 0 is in the pvalue column by 0.00000000001
all[all$Pvalue > 0.05, "Pvalue"] <- 1
all[all$Pvalue <= 0.05 & all$Pvalue > 0.001, "Pvalue"] <- 0.05
all[all$Pvalue <= 0.001 & all$Pvalue > 0.00001, "Pvalue"] <- 0.001
all[all$Pvalue == 0, "Pvalue"] <- 0.00001

# filter the dataframe according to significance
# all <- all[all$Pvalue < 0.05,]

# reorder file to make it consistent
all <- all[,c(5,4,1:3)]
colnames(all) <- c("sample","hallmark","direction","pVal","adjp")
unique(all$sample)

# make a pvalue column with - and + according to up and down-regulation
all$pVal_dir <- all$pVal*all$direction

# Remove for plot non-essential columns
all[,c("direction","pVal","adjp")] <- NULL

# tranpose with Hallmarks as columns instead
all.wide <- all %>% spread(hallmark, pVal_dir,fill=1) # Go from tall to wide format
rownames(all.wide) <- all.wide$sample # Make samples name as row name

# remove sample column
all.wide$sample <- NULL

# Transpose
all.wide.t <- as.data.frame(t(all.wide))

# For plotting purpose add a column for the RG we're comparing with. Seq8: RG1, Seq18: RG2 (actually RG2c, but here called RG2)
all.wide.t$subclone2c <- 1

### KEEP AND ORDER GENE SETS YOU WANT TO PLOT

### GENE SET TO PLOT FOR SEQ 18
keeps <- c("HALLMARK_HYPOXIA","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_ALLOGRAFT_REJECTION",
           "HALLMARK_E2F_TARGETS","HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

# subset the matrix to the pathways of interest
sub.all.wide.t <- all.wide.t[keeps,]

### CONVERT TO MATRIX TO COMPLY WITH HEATMAP PROGRAM
allm <- data.matrix(sub.all.wide.t)

# subset and order allm
allm <- allm[,c("subclone0a", "subclone0b", "subclone2a", "subclone2b",  "subclone2c", "subclone7", "subclone8", "subclone3", "subclone4", "subclone1")]

### ADD IN BREAKS AND COLORS (FOR PVALUE)
colors = c(-1,-0.05,-0.001,-0.00001,0,0.00001,0.001,0.05,1)
my_palette <- c("white","#b3cde0","#005b96","#011f4b","black","#a70000","#ff0000","#ffbaba","white")

### for merging marker files
for (i in 1:length(marker.files)){
  
  # add the hallmark column to the file
  file <- marker.files[[i]]
  file$genes <- rownames(file)
  
  # add a subclone column
  file$subclone <- subclone.id[i]
  
  # exchange the files
  marker.files[[i]] <- file
}

# combine all the files together
marker.file.df <- bind_rows(marker.files)
marker.file.df <- marker.file.df[marker.file.df$p_val < 0.05,]

# get the number of differential expressed genes per subclone
marker.df <- data.frame(table(marker.file.df$subclone))
rownames(marker.df) <- marker.df$Var1
marker.df$Var1 <- NULL

# manually add subclone 1
marker.df[nrow(marker.df)+1,] <- 0
rownames(marker.df)[nrow(marker.df)] <- "subclone2c"

# convert to matrix and reorder
marker.matrix <- as.matrix(marker.df)
marker.matrix <- as.matrix(marker.matrix[colnames(allm),])

# get colors for each bar plot
ECB.colors <- distinct(seurat.obj@meta.data[seurat.obj$subclones %in% str_split_fixed(colnames(allm), "subclone", 2)[,2], c("subclones", "RG_color")])
ECB.colors$subclones <- paste0("subclone", ECB.colors$subclones)
rownames(ECB.colors) <- ECB.colors$subclones
ECB.colors <- ECB.colors[colnames(allm),]

# generate the first heatmap based on the score matrix
ht1 <- Heatmap(allm, 
               name = "score", 
               col = circlize::colorRamp2(colors, my_palette), 
               cluster_rows = F, 
               cluster_columns = FALSE,
               show_row_names = F, 
               row_names_side = "left",
               show_row_dend = FALSE,
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 9), 
               rect_gp = gpar(col = "black", lwd = 1), 
               column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               heatmap_legend_param = list(color_bar = "discrete",
                                           at = colors,
                                           title = ""),
               top_annotation = HeatmapAnnotation(column_barplot = anno_barplot(marker.matrix, 
                                                                                border = TRUE, 
                                                                                height = unit(1.5, "cm"),
                                                                                annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                                                                                bar_width = 0.75, 
                                                                                gp = gpar(col = "white", fill = ECB.colors$RG_color, fontsize = 10, fontface = "bold"), 
                                                                                axis_param = list(at = c(0, 100, 200, 300),
                                                                                                  labels = c("0", "100", "200", "300")),
                                                                                width = unit(2, "cm")),
                                                  show_annotation_name = F), 
               border = T)

pdf(paste0(sample.tmp, "/DEG/GSEA_heatmap_pval_ordered_hallmarks.pdf"),width=6,height=5,pointsize=0.1)
print(ht1)
dev.off()

###########################################################################
#         GENERATE INDIVIDUAL GSEA PLOTS HERE FOR 0A VS. OTHERS
###########################################################################
# subset the GSEA.results to 0a vs. 1
clone0.markers <- read.table(paste0(sample.tmp, "/DEG/markers_subclone0a_vs_2c.txt"), stringsAsFactors = F, header = T)
clone0.markers <- subset(clone0.markers, p_val <= 0.05)

phenotype <- as.vector(clone0.markers$avg_logFC)
names(phenotype) <- rownames(clone0.markers)
names(phenotype)

#print (phenotype)
## specify the gene sets type you want to analyze
MSig_H <- MSigDBGeneSets(species = "Hs", collection = "H", subcategory = NULL) # Hallmarks!
ListGSC <- list(MSig_H=MSig_H)

## iniate a *GSCA* object
gsca <- GSCA(listOfGeneSetCollections=ListGSC, 
             geneList=phenotype)

## preprocess
gsca1 <- preprocess(gsca, species="Hs", initialIDs="SYMBOL",
                    keepMultipleMappings=TRUE, duplicateRemoverMethod="max",
                    orderAbsValue=FALSE)

## analysis
if (requireNamespace("doParallel", quietly=TRUE)) {doParallel::registerDoParallel(cores=4)}  ## support parallel calculation using multiple cores

# perform analysis here
gsca2 <- analyze(gsca1, para=list(pValueCutoff=0.051, pAdjustMethod="BH", nPermutations=10000, minGeneSetSize=1, exponent=1), doGSOA = FALSE)

## append gene sets terms
gsca3 <- appendGSTerms(gsca2, msigdbGSCs=c("MSig_H"))

## draw GSEA plot for a specific gene set
topGS <- getTopGeneSets(gsca3, resultName="GSEA.results",
                        gscs=c("MSig_H"),allSig=TRUE)#, allSig=TRUE)

## IL6
pdf(paste0(sample.tmp, "/DEG/GSEA_EnrichmentPlot_subclone0a_vs_2c_IL6.pdf"),width=5,height=4,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_IL6_JAK_STAT3_SIGNALING", main.title = "IL6 JAK STAT3 Signaling")
dev.off()

## G2M
pdf(paste0(sample.tmp, "/DEG/GSEA_EnrichmentPlot_subclone0a_vs_2c_Hypoxia.pdf"),width=5,height=4,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_HYPOXIA", main.title = "Hypoxia")
dev.off()

## TNFA
pdf(paste0(sample.tmp, "/DEG/GSEA_EnrichmentPlot_subclone0a_vs_2c_Angiogenesis.pdf"),width=5,height=4,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_ANGIOGENESIS", main.title = "Angiogenesis")
dev.off()

## MYC
pdf(paste0(sample.tmp, "/DEG/GSEA_EnrichmentPlot_subclone0a_vs_2c_TNFA.pdf"),width=5,height=4,pointsize=0.1)
viewGSEA(gsca3, gscName="MSig_H", gsName="HALLMARK_TNFA_SIGNALING_VIA_NFKB", main.title = "TNFa signaling via NFKß")
dev.off()

# ## view enrichment Map
# viewEnrichMap(gsca3, gscs=c("MSig_H"),
#               ntop = 5, gsNameType = "id")

###########################################################################
#         GENERATE INDIVIDUAL GSEA PLOTS HERE FOR 0A VS. OTHERS
###########################################################################

# number of expressed genes in the underlying matrix
genes_per_arm <- as.data.frame(table(SCT.pseudobulk.comb$chrarm))
colnames(genes_per_arm) <- c("chrarm","nrGenes")

### CREATE DATA FRAME WITH CHR ARM AND ZERO (COUNTS TO BE USED LATER)
chrarm <- unique(genes_per_arm$chrarm)
chrdf0 <- as.data.frame(chrarm)
chrdf0$value <- 0

### Chromosomes you want to plot
keeps <- c("chr3p","chr3q", "chr9p","chr9q", "chr11p","chr11q", "chr15p","chr15q","chr20p","chr20q")

### CALCULATE MEDIAN PROPORTION ACROSS ALL CHROMOSOMES FOR THE RGS OF INTEREST
propup <- c()
propdown <- c()
rgs <- subclone.id
for (rg in rgs){
  print (rg)
  propup <- getAveAltUp(rg,propup)
  propdown <- getAveAltDown(rg,propdown)
  #print (propup)
  print (propdown)
}
mean(propup)
median(propup)
mean(propdown)
median(propdown)


### SPECIFY RG OF INTEREST
rgs <- subclone.id
for (rg in rgs){
  print (rg)
  plotDEGs(rg)
}


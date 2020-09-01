#############################################################################################################################
##                                                                                                                      
##  GENERATE SCATTERPLOTS FROM SHALLOW WGS RESULTS PROCESSED USING QDNASEQ
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
list.of.packages <- c("reshape2", "optparse", "BSgenome", "RColorBrewer", "ggplot2", "scales", "DescTools", "dendextend", "tidyverse", 
                      "Matrix", "devtools", "Matrix.utils", "matrixStats", "readr", "magrittr", "fishplot", "Signac", "BiocManager", 
                      "biomaRt", "httr", "ComplexHeatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages, repos = "http://cran.us.r-project.org")
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

# add function
`%notin%` <- Negate(`%in%`)

#####################################################################################
# READ IN DATA OF INTEREST FROM SHALLOW WGS
#####################################################################################
# set output directory
o.dir <- "/labs/ccurtis2/mjprzy/correlation_test/"

# set input directory and collect the files
wgs.files <- list.files("/labs/ccurtis2/kasper/Clonal_evolution/WGS/wgs_rest/Mapping_hg19", pattern = "sorted_CNV_Smooth_Segment_per_bin.txt|sorted_bin50kb_CNV_Smooth_Segment_per_bin.txt", recursive = T, full.names = T)
wgs.files <- wgs.files[grep("6077_P_C5_Rep2", wgs.files)]
wgs.files2 <- list.files("/labs/cjkuo/Kasper/WGS/WGS6/Mapping_hg19", pattern = "sorted_CNV_Smooth_Segment_per_bin.txt|sorted_bin50kb_CNV_Smooth_Segment_per_bin.txt", recursive = T, full.names = T)
wgs.files2 <- wgs.files2[grep("ECB_EARLY_6077_P_C5_R2_T12_sorted_CNV", wgs.files2)]

# combine lists and only get the P2C5R2 files
all.wgs.files <- c(wgs.files, wgs.files2)

# set the sample ids
sample.ids <- c("P2C2R2T3", "P2C2R2T8", "P2C2R2T12")
timepoints <- c("T3", "T8", "T12")

all.wgs.files <- lapply(all.wgs.files, read.table, stringsAsFactors = F, header = T)

i <- 1
for (i in 1:length(sample.ids)){
  
  # add a sample id column
  table <- all.wgs.files[[i]]
  colnames(table)[ncol(table)] <- "logR"
  table$sample_id <- timepoints[i]
  all.wgs.files[[i]] <- table
  
}

# bind them together
wgs.table <- bind_rows(all.wgs.files)
  
# transform the data to get approx cnv states
wgs.logR<- 2+(wgs.table$logR)
wgs.cn.states <- wgs.logR

# check length of each sample
table(wgs.table$sample_id)

# set up a dataframe with all information needed
wgs.data <- data.frame(cn_state=wgs.cn.states, chr=wgs.table$chromosome, coordinates=c(1:49857,1:49857, 1:49857), sample_id = wgs.table$sample_id)

# Upperlimit the data to 8 copies
wgs.data$cn_state[which(wgs.data$cn_state>8)] <- 8
wgs.data$cn_state[which(wgs.data$cn_state < 0)] <- 0

#This is needed to plot the chromosomes in the right order
wgs.data$order <- c(1:49857,1:49857, 1:49857)
wgs.data$sample_id <- factor(wgs.data$sample_id, levels = c("T3", "T8", "T12"))

# create a ggplot for this
plotbulk<- ggplot(wgs.data, aes(coordinates, cn_state)) +
  geom_point(aes(colour = cn_state)) +
  scale_colour_gradient2(low = "darkblue", high = "darkred", mid ="lightgray", midpoint = 2) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        strip.text = element_text(face="bold", size=20, colour = "black",),
        strip.background = element_rect(fill="white", colour="black", size=1), 
        axis.text = element_text(colour = "black", size = 20, face = "bold" ),
        axis.title = element_text(colour = "black", size = 25, face = "bold" ),
        plot.title = element_text(colour = "black", size = 25, face = "bold", hjust = 0.5),
        legend.title = element_text(color = "black", size = 22, face = "bold",),
        legend.text = element_text(colour="black", size=20, face="bold"),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing.x=unit(0, "lines"), panel.border = element_rect(linetype =3))+
  facet_grid(sample_id ~ reorder(chr,order), scales="free", space="free") + 
  ggExtra::removeGrid() + labs(x="Chromosomes", y="", color = "Copy Number")+
  scale_x_continuous(expand = c(0.01, 0.01)) 

pdf(paste0(o.dir, "P2C2R2_WGS_Scatterplot.pdf"), width = 28, height = 12)
print(plotbulk + NoLegend())
dev.off()




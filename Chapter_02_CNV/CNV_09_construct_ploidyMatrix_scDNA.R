#########################################################################################
###  
###   Create a bins x cell matrix from 10X scCNV data for D2R2 replicate 2.
###  
###   Author: Hana Susak, adapted by Moritz Przybilla
###   Date: 16/05/2020
###   Notes: This script reads in the sample names, from 10x CNV. It further needs a path
###          to the bed files from CNV incorporating the segment with CN calls. 
###          First, noisy cells in the data are removed. Afterwards, the minimum segments
###          are determined and intersected across cells. Lastly, a matrix with minimum
###          segments x cells is constructed.
###
#########################################################################################

rm(list = ls()) 
set.seed(16011985) # set the seed of random number generator 

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("reshape2", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE))) # ignore all "simple" diagnostic messages (warnings or errors)

samples <- trimws(unlist(strsplit(c('6077_P_C5_ECB_R2_T2'), split=',' )))
samples.seq.type <- trimws(unlist(strsplit( c('10x'), split=',' )))
samples.cnv.path <- trimws(unlist(strsplit( c('/labs/ccurtis2/mjprzy/scDNA_6077_C2R2/node_cnv_calls.bed'), split=',' )))
out.dir <- paste0('/labs/ccurtis2/mjprzy/correlation_test', '/')

###################################################################################
# load cells CNV calls and remove noisy cells
###################################################################################

info.per.cell.all <- data.frame() #only 10x data
node.cnv.calls.all <- data.frame()
# i <- 1
for (i in 1:length(samples)){
  print(samples[i])
  
  # assign tmp variable for the current sample with all information on seq and cnv path
  sample.tmp <- samples[i]
  sample.seq.type.tmp <- samples.seq.type[i]
  sample.cnv.path.tmp <- samples.cnv.path[i]
  
  node.cnv.calls <- read.table( sample.cnv.path.tmp )
  if(sample.seq.type.tmp=='10x'){
    # add original colnames present in the bed file
    colnames(node.cnv.calls) <- c('chrom', 'start', 'end',	'id',	'copy_number',	'event_confidence')
  }
  # rename the ids to sampleName_previousID
  node.cnv.calls$sample <- sample.tmp
  node.cnv.calls$id <- paste(sample.tmp,node.cnv.calls$id, sep='_')
  
  # 10x only! Calculate and filter noisy cells
  if(sample.seq.type.tmp=='10x'){
    # read in cells.info provided by 10x and get the number of cells
    # df with total_num_reads frac_bases_R1_Q30 frac_bases_R2_Q30 correct_bc_rate frac_non_cell_barcode etc
    cells.info <- read.table( gsub('node_cnv_calls.bed','summary.csv',sample.cnv.path.tmp), sep = ',', header = T)
    num.cells <- as.numeric(cells.info['num_cells'])
    rm(cells.info)
    
    # read in summary metrics for all cells
    # df with barcode cell_id total_num_reads num_unmapped_reads  etc.
    info.per.cell <- read.table(gsub('node_cnv_calls.bed','per_cell_summary_metrics.csv',sample.cnv.path.tmp) , header = T, sep = ",", stringsAsFactors = F ) 
    node.cnv.calls <-  node.cnv.calls[node.cnv.calls$id %in% paste(sample.tmp, info.per.cell$cell_id, sep='_'),] # keep info on cells which are in node.cnv.calls
    
    ## calculating new noisy cells
    num.seg <- list()
    # j <- 1
    for (j in 1:(num.cells)){ 
      # print(j)
      node.cnv.calls.tmp <- node.cnv.calls[node.cnv.calls$id==paste(sample.tmp, (j-1), sep='_') ,]
      num.seg[j] <- nrow(node.cnv.calls.tmp)
    } 
    rm(node.cnv.calls.tmp)
    rm(j)
    # list with number of segments per cell in the noisy cells
    num.seg <- unlist(num.seg)
    
    info.per.cell$num_seg <- num.seg
    info.per.cell$new_noisy <- 0
    rm(num.seg)
    
    # calculate standard devition on number of segments
    sd.num_seg <- mean(info.per.cell$num_seg) + sd(info.per.cell$num_seg)
    sd2.num_seg <- mean(info.per.cell$num_seg) + 2*sd(info.per.cell$num_seg)
    
    ### cell is noisy if high_dimapd and number of segments greater then mean+sd
    info.per.cell[info.per.cell$is_high_dimapd == 1  & info.per.cell$num_seg > sd.num_seg ,'new_noisy'] <- 1
    
    ### cell is noisy if ploidy_confidence low (between 0 and 2) and number of segments greater then mean+sd
    info.per.cell[(info.per.cell$ploidy_confidence >=0 & info.per.cell$ploidy_confidence <=2 )  & info.per.cell$num_seg > sd.num_seg ,'new_noisy'] <- 1
    
    ### cell is noisy if number of segments greater then mean+2*sd
    info.per.cell[info.per.cell$num_seg > sd2.num_seg ,'new_noisy'] <- 1
    
    rm(sd.num_seg)
    rm(sd2.num_seg)
    
    ###################
    info.per.cell$id <- paste(sample.tmp, info.per.cell$cell_id, sep='_')
    info.per.cell$sample <- sample.tmp
    # remove the new noisy cells and keep only those which are assigned as 0
    node.cnv.calls <- node.cnv.calls[node.cnv.calls$id %in% info.per.cell[info.per.cell$new_noisy==0, 'id'], ]
    # print(dim(info.per.cell))
    # print(dim(node.cnv.calls))
    
    info.per.cell.all <- rbind.data.frame(info.per.cell.all, info.per.cell )
    
  } else { # here, code for removing noisy cells in G&T and basically everthin else which is not 10x CNV is missing
    ## HANA!! add code to remove noisy/control cells
    noisy.cells <- c('C5', 'B3', 'B12', 'C12')
    control.cells <- c('A1', 'H11', 'H12', 'H13')
    
    node.cnv.calls <- node.cnv.calls[! node.cnv.calls$id %in% c( paste(sample.tmp,noisy.cells, sep='_') ),]
    node.cnv.calls <- node.cnv.calls[! node.cnv.calls$id %in% c( paste(sample.tmp,control.cells, sep='_') ),]
    
  }
  
  node.cnv.calls.all <- rbind.data.frame(node.cnv.calls.all, node.cnv.calls )
  node.cnv.calls.all$chrom <-  gsub('chr','',node.cnv.calls.all$chrom )
  
  rm(sample.tmp)
  rm(sample.seq.type.tmp)
  rm(sample.cnv.path.tmp)
}
rm(i)
rm(info.per.cell)
rm(node.cnv.calls)

# Find bed files with all intersects
setwd(out.dir)

# get cell ids of all cells which we have bed file information for
cells.ids <- unique(node.cnv.calls.all$id)
cells.ids <- sort(cells.ids)
cells.cnvs <- node.cnv.calls.all[node.cnv.calls.all$id %in% cells.ids,  ] #? Hana:simplify # only cells which are not noisy are kept 
cells.cnvs <- cells.cnvs[with(cells.cnvs,order(id,chrom,start,end)), ] #? Hana:simplify # ordered bed files 
num.cells <- length(unique(cells.cnvs$id)) # total number of cells

###################################################################################
### Creating/Loading minimum segments intersects between all cells
###################################################################################
# setwd("/home/przybilm/cells_ploidy_matrix_test/")
setwd(out.dir)

if(!file.exists('min.segments.bed')){
  for (i in 1:(num.cells)){ 
    # i <- 2
    # make a dataframe for the first cell, get its information from the cnv bedfile and sort the file
    if(i==1){
      cell.id.tmp <- cells.ids[i]
      cells.cnvs.tmp <- cells.cnvs[cells.cnvs$id==cell.id.tmp, ]
      write.table(cells.cnvs.tmp, file = 'cellA.bed', col.names=F, row.names = F, quote = F, sep = '\t')
      system(paste0('sort -V -k1,1 -k2,2n cellA.bed > cellA.sorted.bed')) #-V -k1,1 -k2,2
      while (!file.exists('cellA.sorted.bed')) {
        Sys.sleep(1)
      }
      system(paste0('rm cellA.bed'))
      
    } else {
      cell.id.tmp <- cells.ids[i]
      cells.cnvs.tmp <- cells.cnvs[cells.cnvs$id==cell.id.tmp, ]
      write.table(cells.cnvs.tmp, file = 'cellB.bed', col.names=F, row.names = F, quote = F, sep = '\t')
      while (!file.exists('cellB.bed')) {
        Sys.sleep(1)
      }
      system(paste0('sort -V -k1,1 -k2,2n cellB.bed > cellB.sorted.bed'))
      while (!file.exists('cellB.sorted.bed')) {
        Sys.sleep(1)
      }
      system(paste0('rm cellB.bed'))
      
      # intersect cellA and cellB files 
      system(paste0('/scg/apps/software/bedtools2/2.27.1/bin/intersectBed -sorted -wb -a cellA.sorted.bed -b cellB.sorted.bed > cellAnew.bed'))
      while (!file.exists('cellAnew.bed')) {
        Sys.sleep(1)
      }
      system(paste0('rm cellA.sorted.bed'))
      while (file.exists('cellA.sorted.bed')) {
        Sys.sleep(1)
      }
      system(paste0('rm cellB.sorted.bed'))
      
      #new.intersect <- read.table('cellAnew.bed', stringsAsFactors = F)
      system(paste0('cut -f1-3 cellAnew.bed > cellAnew2.bed')) # could go also 1-7 I think # get the bed file with min segments
      system(paste0('rm cellAnew.bed'))
      
      # sort the min segment file
      system(paste0('sort -V -k1,1 -k2,2n cellAnew2.bed > cellA.sorted.bed'))
      while (!file.exists('cellA.sorted.bed')) {
        Sys.sleep(1)
      }
      system(paste0('rm cellAnew2.bed'))
    }
    #print(cell.id.tmp)
    rm(cell.id.tmp)
    rm(cells.cnvs.tmp)
  }
  
  # the intersection between all cells is the the min.segments.bed
  inter.bed <- read.table('cellA.sorted.bed')
  inter.bed <- inter.bed[,c(1:3)]
  write.table(inter.bed, file = 'min.segments.bed', col.names=F, row.names = F, quote = F, sep = '\t')
  while (!file.exists('min.segments.bed')) {
    Sys.sleep(1)
  }
  system(paste0('rm cellA.sorted.bed'))
  rm(i)
  
} else {
  # read final segements
  inter.bed <- read.table('min.segments.bed') # load bed file with intersected regions across cells (header #chrom, start, end)
  
}

###################################################################################
### get cells cnv calls for each min segment
###################################################################################

fin.cnv.df <- data.frame(id.pos=paste0(inter.bed$V1,':', inter.bed$V2, '-', inter.bed$V3, sep='')) # merge the chrom:start-end together
fin.cnv.conf.df <- data.frame(id.pos=paste0(inter.bed$V1,':', inter.bed$V2, '-', inter.bed$V3, sep=''))
order.cnv <- as.character(fin.cnv.df$id.pos)
fin.cnv.df$id.pos <- factor(fin.cnv.df$id.pos , levels = order.cnv)
fin.cnv.conf.df$id.pos <- factor(fin.cnv.conf.df$id.pos , levels = order.cnv)

if(!file.exists('cnv.cells.mat.txt')){
  # i <- 1
  # write a bed file for every cell and sort it 
  for (i in 1:(num.cells)){
    cell.id.tmp <- cells.ids[i]
    print(cell.id.tmp)
    cells.cnvs.tmp <- cells.cnvs[cells.cnvs$id==cell.id.tmp, ]
    write.table(cells.cnvs.tmp, file = 'cell.bed', col.names=F, row.names = F, quote = F, sep = '\t')
    system(paste0('sort -V -k1,1 -k2,2n cell.bed > cell.sorted.bed'))
    while (!file.exists('cell.sorted.bed')) {
      Sys.sleep(1)
    }
    system(paste0('rm cell.bed'))
    
    # intersect the min.segments.bed with each individual cell
    # problem here is the use of the & in G&T
    cell.id.tmp <- gsub("&","n", cell.id.tmp)
    system(paste0('/scg/apps/software/bedtools2/2.27.1/bin/intersectBed -sorted -wb -a min.segments.bed -b cell.sorted.bed > cell_', cell.id.tmp, '.bed'))
    while (!file.exists(paste0('cell_',cell.id.tmp,'.bed')) ){
      Sys.sleep(1)
    }
    system(paste0('rm cell.sorted.bed'))
    
    # read in the cell specific bed
    cell.seg <- read.table(paste0('cell_',cell.id.tmp,'.bed'))
    cell.seg$id.pos <- paste0(cell.seg$V1,':', cell.seg$V2, '-', cell.seg$V3, sep='')
    cell.seg <- cell.seg[,-c(1:3)] # remove the first three columns as they are summarized in id.pos
    colnames(cell.seg) <- c('chr', 'start','end','cell_id', 'copy_number', 'event_confidence', 'sample','id.pos')
    
    # merge the fin.cnv.df data frame with the segments together with the copy number information on the segments
    fin.cnv.df <- merge(x=fin.cnv.df, y=cell.seg[,c('copy_number','id.pos')],  by='id.pos')
    fin.cnv.df <- fin.cnv.df[order(fin.cnv.df$id.pos), ]
    colnames(fin.cnv.df)[i+1] <- paste0('cell_',cell.id.tmp)
    
    # merge the fin.cnv.conf.df data frame with the segments together with the copy number information on the segments
    fin.cnv.conf.df <- merge(x=fin.cnv.conf.df, y=cell.seg[,c('event_confidence','id.pos')],  by='id.pos')
    fin.cnv.conf.df <- fin.cnv.conf.df[order(fin.cnv.conf.df$id.pos),]
    colnames(fin.cnv.conf.df)[i+1] <- paste0('cell_',cell.id.tmp)
    
    rm(cell.seg)
    rm(cell.id.tmp)
    rm(cells.cnvs.tmp)
  }
  system(paste0('rm cell*bed'))
  
  # write the matrix for segments x cells
  write.table(fin.cnv.df, file = 'cnv.cells.mat.txt', col.names=T, row.names = F, quote = F, sep = '\t')
  write.table(fin.cnv.conf.df, file = 'cnv_conf.cells.mat.txt', col.names=T, row.names = F, quote = F, sep = '\t')
  
} else {
  ### load data if already exist
  fin.cnv.df <- read.table('cnv.cells.mat.txt', header = T, stringsAsFactors = F)
  fin.cnv.conf.df <- read.table('cnv_conf.cells.mat.txt', header = T, stringsAsFactors = F)
}

# other information on the cells and their copy number
write.table(info.per.cell.all, file = 'per_cell_summary_metrics.tsv', col.names=T, row.names = F, quote = F, sep = '\t')
write.table(node.cnv.calls.all, file = 'node_cnv_calls.bed', col.names=T, row.names = F, quote = F, sep = '\t')


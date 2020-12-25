source("./utilityFunc/rnaPCAPlot.R")
source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/rnaMiArPCAPlot.R")
source("./utilityFunc/plottingPalettes.R")

library(RColorBrewer)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(ggsci)
options(stringsAsFactors = F)
library(optparse)

option_list <- list( 
  make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-o", "--output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
##NOTE: These plots only show the result after normalization of size factor(account for differences in sequencing depth)
##NOTE: We recommand you run this PCA/Heatmap only when you want to see within-group/between-group variability across multiple groups 
##NOTE: If you want to see the detail DE genes, run Heatmap/Volcano between two groups using "group_by_group"/[expr group] options
if(data_type == "RNAseq"){
  rnaPCAPlot(dataProc = dataProc_RNA, 
             fileLoc = opt$output)
}
if(data_type == "MiAr"){
  rnaMiArPCAPlot(dataProc = dataProc_RNA, 
                 fileLoc = opt$output)
}
##!!!!!DEV: rnaHeatmapPlot input differ in different situation
rnaHeatmapPlot(DEAresult = dataProc_RNA, showallgroups = T, type = data_type, 
               fileLoc = opt$output)

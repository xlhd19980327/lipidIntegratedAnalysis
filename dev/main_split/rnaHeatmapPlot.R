source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")

library(DESeq2)
library(limma)
library(tidyverse)
library(ggrepel)
library(ggpubr)
options(stringsAsFactors = F)
library(RColorBrewer)
library(ggsci)
library(optparse)

option_list <- list( 
  make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-v", "--show_heatmap_top"), action="store", default = 75), 
  make_option(c("-w", "--heatmap_output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
rnaHeatmapPlot(DEAresult, showtop = opt$show_heatmap_top, 
               type = data_type, fileLoc = opt$heatmap_output)
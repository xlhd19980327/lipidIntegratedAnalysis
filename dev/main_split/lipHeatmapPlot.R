source("./utilityFunc/lipHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")

library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(ggpubr)
options(stringsAsFactors = F)
library(RColorBrewer)
library(ggsci)
library(optparse)

option_list <- list( 
  make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-y", "--heatmap_output"), action="store"), 
  make_option(c("-e", "--top_number"), action="store", default = NA, type = "integer")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
msg.vec <- c()
lipHeatmapPlot(dataSet = dataSet, mSet = mSet, 
               fileLoc = paste0(opt$heatmap_output, "MARresults/"), 
               topnum = opt$top_number)

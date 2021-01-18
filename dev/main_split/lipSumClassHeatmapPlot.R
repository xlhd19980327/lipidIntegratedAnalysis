source("./utilityFunc/lipSumClassHeatmapPlot.R")
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
  make_option(c("-u", "--lipidclass_output"), action="store"), 
  make_option(c("-w", "--ignore_subclass"), action="store", default = T), 
  make_option(c("-z", "--show_detail"), action="store", default = F)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
lipSumClassHeatmapPlot(dataSet = dataSet, mSet = mSet, 
               fileLoc = paste0(opt$lipidclass_output, "headgroup/"), 
               ignore = opt$ignore_subclass, show_detail = opt$show_detail)

source("./utilityFunc/lipVolcanoPlot.R")
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
  make_option(c("-s", "--show_lipid_class"), action="store", default = F), 
  make_option(c("-p", "--volcano_output"), action="store"), 
  make_option(c("-b", "--paired"), action="store", default = F),
  make_option(c("-x", "--p_type"), action="store", default = "raw"), 
  make_option(c("-j", "--fc_thresh"), action="store", default = 2), 
  make_option(c("-k", "--p_thresh"), action="store", default = 0.1), 
  make_option(c("-m", "--show_volcano_top"), action="store", default = 10),
  make_option(c("-w", "--ignore_subclass"), action="store", default = T)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
lipVolcanoPlot(dataSet = dataSet, mSet = mSet, showLipClass = opt$show_lipid_class,
               fileLoc = paste0(opt$volcano_output, "MARresults/"), paired = opt$paired, pval.type = opt$p_type, 
               fcthresh = opt$fc_thresh, pthresh = opt$p_thresh, showtop = opt$show_volcano_top)

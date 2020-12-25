source("./utilityFunc/rnaVolcanoPlot.R")
source("./utilityFunc/plottingPalettes.R")

library(tidyverse)
library(ggrepel)
library(ggpubr)
options(stringsAsFactors = F)
library(RColorBrewer)
library(ggsci)
library(optparse)

option_list <- list( 
  make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-s", "--volcano_output"), action="store"), 
  make_option(c("-f", "--fc_thresh"), action="store", default = 2.0), 
  make_option(c("-p", "--p_thresh"), action="store", default = 0.1), 
  make_option(c("-u", "--show_volcano_top"), action="store", default = 20)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
rnaVolcanoPlot(DEAresult = DEAresult, type = data_type, fcthresh = opt$fc_thresh, 
               pthresh = opt$p_thresh, showtop = opt$show_volcano_top, 
               fileLoc = opt$volcano_output)

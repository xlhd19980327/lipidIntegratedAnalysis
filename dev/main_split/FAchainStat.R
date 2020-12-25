source("./utilityFunc/FAchainStat.R")
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
  make_option(c("-v", "--fachains_output"), action="store"), 
  make_option(c("-g", "--plot_type"), action="store", default = "FA_info"),
  make_option(c("-w", "--ignore_subclass"), action="store", default = T)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
FAchainStat(dataSet = dataSet, mSet = mSet, 
            fileLoc = paste0(opt$fachains_output, "FAchainVisual/"), 
            plotInfo = opt$plot_type, ignore = opt$ignore_subclass)
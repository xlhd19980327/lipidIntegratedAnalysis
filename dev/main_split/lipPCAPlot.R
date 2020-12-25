source("./utilityFunc/lipPCAPlot.R")
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
  make_option(c("-q", "--pca_output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
lipPCAPlot(dataSet = dataSet, mSet = mSet, 
           fileLoc = paste0(opt$pca_output, "MARresults/"))


library(dplyr)
library(tibble)
library(optparse)
options(stringsAsFactors = F)
source("./utilityFunc/LIONenrich_tar.R")

option_list <- list( 
  make_option(c("-i", "--data_file"), action="store"),
  make_option(c("-j", "--row"), action="store", type = "integer"), 
  #make_option(c("-k", "--colum"), action="store", type = "integer"), 
  
  make_option(c("-o", "--output_loc"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dataloc <- paste0(opt$data_file, 
                  dir(opt$data_file, pattern = paste0("lipids_", opt$row, "\\.csv")))
lipid_target_list <- read.csv(dataloc)$x
lipid_background_list <- row.names(read.csv(paste0(opt$data_file, "correlation.csv"), row.names = 1))
LIONenrich_tar(lipid_target_list = lipid_target_list, lipid_background_list = lipid_background_list, 
               fileLoc = opt$output_loc, reg = "none")

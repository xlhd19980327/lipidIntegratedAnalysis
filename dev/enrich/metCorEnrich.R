library(MetaboAnalystR)
library(optparse)
options(stringsAsFactors = F)
source("./utilityFunc/metEnrichFunc.R")

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
metEnrichFunc(csvfile = dataloc, filename = "none", 
              fileLoc = opt$output_loc)
library(optparse)
options(stringsAsFactors = F)
source("./utilityFunc/plotCircos.R")

option_list <- list( 
  make_option(c("-r", "--godata_file"), action="store"),
  make_option(c("-i", "--data_file"), action="store"),
  make_option(c("-j", "--row"), action="store", type = "integer"), 
  make_option(c("-k", "--colum"), action="store", type = "integer"), 
  
  make_option(c("-t", "--thresh_cor"), action="store", default = 0.6), 
  make_option(c("-n", "--go_topnum"), action="store", default = 20), 
  make_option(c("-o", "--outfileLoc"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

invisible(plotCircos(gofileLoc = opt$godata_file, corfileLoc = opt$data_file, 
           k = opt$colum, j = opt$row, thresh1 = opt$thresh_cor, 
           topnum = opt$go_topnum, outfileLoc = opt$outfileLoc))


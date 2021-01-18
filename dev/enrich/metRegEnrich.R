source("./utilityFunc/metEnrichFunc.R")

library(MetaboAnalystR)
library(optparse)

option_list <- list( 
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-o", "--metenrich_output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

up <- read.csv(paste0(opt$input_file, "up.csv"))$x
down <- read.csv(paste0(opt$input_file, "down.csv"))$x
metEnrichFunc(csvfile = paste0(opt$input_file, "up.csv"), filename = "up", 
                           fileLoc = opt$metenrich_output)
metEnrichFunc(csvfile = paste0(opt$input_file, "down.csv"), filename = "down", 
                           fileLoc = opt$metenrich_output)

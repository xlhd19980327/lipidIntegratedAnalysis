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
if(length(up) == 0){
  cat("No UP regulation metabolites found! Try reduce the threshold of Fold change/p.value!\n")
}else{
  metEnrichFunc(csvfile = paste0(opt$input_file, "up.csv"), filename = "up", 
                           fileLoc = opt$metenrich_output)
}
if(length(down) == 0){
  cat("No DOWN regulation metabolites found! Try reduce the threshold of Fold change/p.value!\n")
}else{
  metEnrichFunc(csvfile = paste0(opt$input_file, "down.csv"), filename = "down", 
                fileLoc = opt$metenrich_output)
}


source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/MARpreproc.R")

options(stringsAsFactors = F)
library(optparse)
library(tidyverse)
library(MetaboAnalystR)

option_list <- list( 
  make_option(c("-a", "--analysis_option"), action="store"),
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  make_option(c("-c", "--control_group"), action="store", default = ""), 
  make_option(c("-t", "--data_type"), action="store"), 
  make_option(c("-f", "--feature_field"), action="store", default = NA, type = "character"), 
  make_option(c("-l", "--lipid_odd_chain_deletion"), action="store", default = F), 
  make_option(c("-o", "--tidy_output"), action="store"), 
  make_option(c("-n", "--NA_string"), action="store", default = NULL, type = "character"), 
  
  make_option(c("-p", "--output_temp"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

analOpt <- opt$analysis_option
dataSet <- readingLipidData(datafile = opt$input_file, sampleList = opt$description_file, 
                            controlGrp = opt$control_group, dataType = opt$data_type, 
                            lipField = opt$feature_field, delOddChainOpt = opt$lipid_odd_chain_deletion,
                            na.char = opt$NA_string)
if(analOpt == "all_together"){
  cat("All groups will be analyzed together\n")
  mSet <- MARpreproc(dataSet = dataSet, fileLoc = opt$tidy_output)
}else{
  prepDataSet <- function(x, dataset = dataSet){
    ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
    ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
    dataset$data <- dataset$data[, ind]
    dataset$groupsLevel <- dataset$groupsLevel[ind2]
    dataset$allgroups <- dataset$allgroups[ind]
    return(dataset)
  }
  
  cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
  dataSet <- prepDataSet(analOpt)
  mSet <- MARpreproc(dataSet = dataSet, fileLoc = opt$tidy_output)
}
save(mSet, dataSet, file = paste0(opt$output_temp, "data.RData"))

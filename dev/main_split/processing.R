source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/MARpreproc.R")

options(stringsAsFactors = F)
library(optparse)
library(tidyverse)
library(MetaboAnalystR)

option_list <- list( 
  #make_option(c("-a", "--analysis_option"), action="store"),
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  make_option(c("-c", "--control_group"), action="store", default = ""), 
  make_option(c("-t", "--data_type"), action="store"), 
  #make_option(c("-f", "--feature_field"), action="store", default = NA, type = "character"), 
  make_option(c("-l", "--lipid_odd_chain_deletion"), action="store", default = F), 
  make_option(c("-o", "--tidy_output"), action="store"), 
  #make_option(c("-n", "--NA_string"), action="store", default = NULL, type = "character"), 
  make_option(c("-e", "--na_percent"), action="store", default = 67.0, type = "double"), 
  
  make_option(c("-p", "--output_temp"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

analOpt <- opt$analysis_option
dataSet <- readingLipidData(datafile = opt$input_file, sampleList = opt$description_file, 
                            controlGrp = opt$control_group, dataType = opt$data_type,
                            delOddChainOpt = opt$lipid_odd_chain_deletion)
#if(analOpt == "all_together"){
#  cat("All groups will be analyzed together\n")
#  if(length(dataSet[["groupsLevel"]]) != 2){
#    write.csv(FALSE, paste0(opt$output_temp, "show_enrich.csv"), row.names = F)
#  }else{
#    write.csv(TRUE, paste0(opt$output_temp, "show_enrich.csv"), row.names = F)
#  }
#  mSet <- MARpreproc(dataSet = dataSet, fileLoc = opt$tidy_output)
#}else{
#  write.csv(TRUE, paste0(opt$output_temp, "show_enrich.csv"), row.names = F)
#  prepDataSet <- function(x, dataset = dataSet){
#    ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
#    ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
#    dataset$data <- dataset$data[, ind]
#    dataset$groupsLevel <- dataset$groupsLevel[ind2]
#    dataset$allgroups <- dataset$allgroups[ind]
#    return(dataset)
#  }
#  
#  cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
#  dataSet <- prepDataSet(analOpt)
  mSet <- MARpreproc(dataSet = dataSet, fileLoc = opt$tidy_output, 
                     perc = opt$na_percent)
#}
data_type <- dataSet$dataType
data_proc <- t(mSet[["dataSet"]][["proc"]])
dataSet$lipidName <- rownames(data_proc)
rownames(data_proc) <- NULL
data_proc <- as.data.frame(data_proc)
dataSet$data <- data_proc
save(mSet, dataSet, data_type, file = paste0(opt$output_temp, "data.RData"))

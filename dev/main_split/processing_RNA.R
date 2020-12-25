source("./utilityFunc/readingRNAData.R")
source("./utilityFunc/DESeq2preproc.R")
source("./utilityFunc/limmaPreproc.R")
source("./utilityFunc/DEAnalysis.R")
source("./utilityFunc/DEAnalysis_limma.R")

library(tidyverse)
options(stringsAsFactors = F)
library(optparse)
library(DESeq2)
library(limma)
library(apeglm)

option_list <- list( 
  make_option(c("-a", "--analysis_option"), action="store"),
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  make_option(c("-c", "--control_group"), action="store", default = ""), 
  make_option(c("-o", "--tidy_output"), action="store"), 
  make_option(c("-n", "--norm"), action="store", default = T), 
  
  make_option(c("-t", "--data_type"), action="store", default = "RNAseq"), 
  make_option(c("-p", "--output_temp"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dataSet_RNA <- readingRNAData(datafile = opt$input_file, 
                              sampleList = opt$description_file, 
                              controlGrp = opt$control_group)
if(opt$data_type == "RNAseq"){
  dataProc_RNA <- DESeq2preproc(dataSet = dataSet_RNA, 
                                fileLoc = opt$tidy_output)
  DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = opt$analysis_option, 
                          fileLoc = opt$tidy_output)
}
if(opt$data_type == "MiAr"){
  dataProc_RNA <- limmaPreproc(dataSet = dataSet_RNA, norm = opt$norm, 
                               fileLoc = opt$tidy_output)
  DEAresult <- DEAnalysis_limma(dataProc = dataProc_RNA, experGrp = opt$analysis_option, 
                                fileLoc = opt$tidy_output)
}
sampleList <- opt$description_file
sampleInfo <- read.csv(sampleList)
groupsLevel <- unique(sampleInfo$conditions)
if(length(groupsLevel) > 2){
  show_variability <- T
}else{
  show_variability <- F
}
data_type <- opt$data_type
write.csv(show_variability, file = paste0(opt$output_temp, "show_variability.csv"), row.names = F)
save(dataProc_RNA, DEAresult, data_type, file = paste0(opt$output_temp, "data.RData"))

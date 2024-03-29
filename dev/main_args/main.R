## Use this script to test function procedure
source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/MARpreproc.R")
source("./utilityFunc/lipVolcanoPlot.R")
source("./utilityFunc/lipPCAPlot.R")
source("./utilityFunc/lipHeatmapPlot.R")
source("./utilityFunc/headgroupStat.R")
source("./utilityFunc/FAchainStat.R")
source("./utilityFunc/plottingPalettes.R")
source("./utilityFunc/statFAChains.R")
source("./utilityFunc/statFAChains_pathAna.R")

library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(ggpubr)
options(stringsAsFactors = F)

#Color
library(RColorBrewer)
library(ggsci)

#Shell interaction
#library(getopt)
library(optparse)

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
  make_option(c("-s", "--show_lipid_class"), action="store", default = F), 
  make_option(c("-p", "--volcano_output"), action="store"), 
  make_option(c("-q", "--pca_output"), action="store"), 
  make_option(c("-y", "--heatmap_output"), action="store"), 
  make_option(c("-e", "--top_number"), action="store", default = NA, type = "integer"), 
  make_option(c("-u", "--lipidclass_output"), action="store"), 
  make_option(c("-v", "--fachains_output"), action="store"), 
  make_option(c("-g", "--plot_type"), action="store", default = "FA_info"),
  make_option(c("-b", "--paired"), action="store", default = F),
  make_option(c("-x", "--p_type"), action="store", default = "raw"), 
  make_option(c("-j", "--fc_thresh"), action="store", default = 2), 
  make_option(c("-k", "--p_thresh"), action="store", default = 0.1), 
  make_option(c("-m", "--show_volcano_top"), action="store", default = 10),
  make_option(c("-w", "--ignore_subclass"), action="store", default = T)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}
##!!!Client options: can be "all_together"/"group_by_group"/[expr group]
analOpt <- opt$analysis_option
dataSet <- readingLipidData(datafile = opt$input_file, sampleList = opt$description_file, 
                            controlGrp = opt$control_group, dataType = opt$data_type, 
                            lipField = opt$feature_field, delOddChainOpt = opt$lipid_odd_chain_deletion,
                            na.char = opt$NA_string)
if(analOpt == "all_together"){
  cat("All groups will be analyzed together\n")
  mSet <- MARpreproc(dataSet = dataSet, fileLoc = opt$tidy_output)
  lipPCAPlot(dataSet = dataSet, mSet = mSet, 
             fileLoc = opt$pca_output)
  lipHeatmapPlot(dataSet = dataSet, mSet = mSet, 
                 fileLoc = opt$heatmap_output, 
                 topnum = opt$top_number)
  if(opt$data_type %in% c("LipidSearch", "MS_DIAL")){
    headgroupStat(dataSet = dataSet, mSet = mSet, 
                fileLoc = opt$lipidclass_output, ignore = opt$ignore_subclass)
    FAchainStat(dataSet = dataSet, mSet = mSet, 
                fileLoc = opt$fachains_output, 
                plotInfo = opt$plot_type, ignore = opt$ignore_subclass)
  }
}else if(analOpt == "group_by_group"){
  cat("Group-by-group analysis mode\n")
  for(i in dataSet$groupsLevel[dataSet$groupsLevel != dataSet$controlGrp]){
    dataset <- prepDataSet(i)
    
    mSet <- MARpreproc(dataSet = dataset, fileLoc = opt$tidy_output)
    lipVolcanoPlot(dataSet = dataset, mSet = mSet, showLipClass = opt$show_lipid_class,
                   fileLoc = opt$volcano_output, paired = opt$paired, pval.type = opt$p_type, 
                   fcthresh = opt$fc_thresh, pthresh = opt$p_thresh, showtop = opt$show_volcano_top)
    lipPCAPlot(dataSet = dataset, mSet = mSet, 
               fileLoc = opt$pca_output)
    lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
                   fileLoc = opt$heatmap_output, 
                   topnum = opt$top_number)
    if(opt$data_type %in% c("LipidSearch", "MS_DIAL")){
      headgroupStat(dataSet = dataset, mSet = mSet, 
                    fileLoc = opt$lipidclass_output, ignore = opt$ignore_subclass)
      FAchainStat(dataSet = dataset, mSet = mSet, 
                  fileLoc = opt$fachains_output, 
                  plotInfo = opt$plot_type, ignore = opt$ignore_subclass)
    }
  }
}else{
  cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
  dataset <- prepDataSet(analOpt)
  
  mSet <- MARpreproc(dataSet = dataset, fileLoc = opt$tidy_output)
  lipVolcanoPlot(dataSet = dataset, mSet = mSet, showLipClass = opt$show_lipid_class,
                 fileLoc = opt$volcano_output, paired = opt$paired, pval.type = opt$p_type, 
                 fcthresh = opt$fc_thresh, pthresh = opt$p_thresh, showtop = opt$show_volcano_top)
  lipPCAPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = opt$pca_output)
  lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = opt$heatmap_output, 
                 topnum = opt$top_number)
  if(opt$data_type %in% c("LipidSearch", "MS_DIAL")){
    headgroupStat(dataSet = dataset, mSet = mSet, 
                  fileLoc = opt$lipidclass_output, ignore = opt$ignore_subclass)
    FAchainStat(dataSet = dataset, mSet = mSet, 
                fileLoc = opt$fachains_output, 
                plotInfo = opt$plot_type, ignore = opt$ignore_subclass)
  }
}



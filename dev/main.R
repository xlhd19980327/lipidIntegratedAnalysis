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
source("./utilityFunc/lipOPLSDAPlot.R")
source("./utilityFunc/lipSumClassHeatmapPlot.R")
source("./utilityFunc/lipSumChainHeatmapPlot.R")

outputLoc <- "./testData/zly_20210918data/output/"
prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}
##!!!Client options: can be "all_together"/"group_by_group"/[expr group]
analOpt <- "preSens"
dataSet <- readingLipidData(datafile = "./testData/zly_20210918data/input/lipid.csv",
                            sampleList = "./testData/zly_20210918data/input/sampleList.csv", 
                            controlGrp = "preRes", dataType = "Lipids", delOddChainOpt = T)
if(analOpt == "group_by_group"){
  cat("Group-by-group analysis mode\n")
  for(i in dataSet$groupsLevel[dataSet$groupsLevel != dataSet$controlGrp]){
    dataset <- prepDataSet(i)
    
    mSet <- MARpreproc(dataSet = dataset, fileLoc = outputLoc)
    data_type <- dataset$dataType
    data_proc <- t(mSet[["dataSet"]][["proc"]])
    dataset$lipidName <- rownames(data_proc)
    rownames(data_proc) <- NULL
    data_proc <- as.data.frame(data_proc)
    dataset$data <- data_proc
    lipVolcanoPlot(dataSet = dataset, mSet = mSet, showLipClass = T,
                   fileLoc = paste0(outputLoc, "MARresults/"))
    lipPCAPlot(dataSet = dataset, mSet = mSet, 
               fileLoc = paste0(outputLoc, "MARresults/"))
    lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
                   fileLoc = paste0(outputLoc, "MARresults/"), 
                   topnum = 75)
    
    headgroupStat(dataSet = dataset, mSet = mSet, ignore = T, 
                  fileLoc = paste0(outputLoc, "headgroup/"))
    FAchainStat(dataSet = dataset, mSet = mSet, ignore = T, 
                fileLoc = paste0(outputLoc, "FAchainVisual/"), 
                plotInfo = "FA_info")
  }
}else if(analOpt == "all_together"){
  cat("All groups will be analyzed together\n")
  mSet <- MARpreproc(dataSet = dataSet, fileLoc = outputLoc, perc = 2/3*100)
  data_type <- dataSet$dataType
  data_proc <- t(mSet[["dataSet"]][["proc"]])
  dataSet$lipidName <- rownames(data_proc)
  rownames(data_proc) <- NULL
  data_proc <- as.data.frame(data_proc)
  dataSet$data <- data_proc
  lipPCAPlot(dataSet = dataSet, mSet = mSet, 
             fileLoc = paste0(outputLoc, "MARresults/"))
  lipOPLSDAPlot(dataSet = dataSet, mSet = mSet, 
                fileLoc = paste0(outputLoc, "MARresults/"))
  lipHeatmapPlot(dataSet = dataSet, mSet = mSet, 
                 fileLoc = paste0(outputLoc, "MARresults/"), 
                 topnum = 75)
  
  headgroupStat(dataSet = dataSet, mSet = mSet, 
                fileLoc = paste0(outputLoc, "headgroup/"))
  lipSumClassHeatmapPlot(dataSet = dataSet, mSet = mSet, show_detail = F,  
                         fileLoc = paste0(outputLoc, "headgroup/"))
  FAchainStat(dataSet = dataSet, mSet = mSet, 
              fileLoc = paste0(outputLoc, "FAchainVisual/"), 
              plotInfo = "all_info", topnum = 25)
  lipSumChainHeatmapPlot(dataSet = dataSet, mSet = mSet, 
                         fileLoc = paste0(outputLoc, "FAchainVisual/"), 
                         plotInfo = "FA_info", show_detail = F)
}else{
  cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
  dataset <- prepDataSet(analOpt)
  
  mSet <- MARpreproc(dataSet = dataset, fileLoc = outputLoc, perc = 2/3*100)
  data_type <- dataset$dataType
  data_proc <- t(mSet[["dataSet"]][["proc"]])
  dataset$lipidName <- rownames(data_proc)
  rownames(data_proc) <- NULL
  data_proc <- as.data.frame(data_proc)
  dataset$data <- data_proc
  lipVolcanoPlot(dataSet = dataset, mSet = mSet, showLipClass = T,
                 fileLoc = paste0(outputLoc, "MARresults/"), ignore = T, showtop = 5)
  lipPCAPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = paste0(outputLoc, "MARresults/"))
  lipOPLSDAPlot(dataSet = dataset, mSet = mSet, 
                fileLoc = paste0(outputLoc, "MARresults/"))
  lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = paste0(outputLoc, "MARresults/"), 
                 topnum = 75)
  
  headgroupStat(dataSet = dataset, mSet = mSet, ignore = T, 
                fileLoc = paste0(outputLoc, "headgroup/"))
  lipSumClassHeatmapPlot(dataSet = dataset, mSet = mSet, show_detail = F,  
                         fileLoc = paste0(outputLoc, "headgroup/"))
  FAchainStat(dataSet = dataset, mSet = mSet, ignore = T, 
              fileLoc = paste0(outputLoc, "FAchainVisual/"), 
              plotInfo = "all_info")
}

##FAchainStat output: lipid_subclass_integStat, use it to do our stat
lipsample <- "COLD"
dataset <- prepDataSet(lipsample)
mSet <- MARpreproc(dataSet = dataset)
stat_res <- FAchainStat(dataSet = dataset, mSet = mSet, 
                        #fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/FAchainVisual/", 
                        #"plotInfo", "stat" can not modify
                        plotInfo = "FA_info", stat = T)
regStat_gene <- statFAChains(lipid_subclass_tidyStat = stat_res, 
             fileLoc = "./branch/benchmark/output/", 
             lipsample = lipsample, spe = "mmu")
regStat_path <- statFAChains_pathAna(lipid_subclass_tidyStat = stat_res, 
                                     fileLoc = "./branch/benchmark/output/", 
                                     lipsample = lipsample, spe = "mmu")

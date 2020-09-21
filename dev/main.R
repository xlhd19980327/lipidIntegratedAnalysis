## Use this script to test function procedure
source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/MARpreproc.R")
source("./utilityFunc/lipVolcanoPlot.R")
source("./utilityFunc/lipPCAPlot.R")
source("./utilityFunc/lipHeatmapPlot.R")
source("./utilityFunc/headgroupStat.R")
source("./utilityFunc/FAchainStat.R")
source("./utilityFunc/plottingPalettes.R")

prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}
analOpt <- "all_together"
dataSet <- readingLipidData(datafile = "./testData/zsy_DGATinhibitors/HeLaData/input/data_tidy_testFormat.csv", 
                            controlGrp = "OA", dataType = "MS_DIAL", delOddChainOpt = T,
                            fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/")
if(analOpt == "group_by_group"){
  cat("Group-by-group analysis mode\n")
  for(i in dataSet$groupsLevel[dataSet$groupsLevel != dataSet$controlGrp]){
    dataset <- prepDataSet(i)
    
    mSet <- MARpreproc(dataSet = dataset)
    lipVolcanoPlot(dataSet = dataset, mSet = mSet, 
                   fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/")
    lipPCAPlot(dataSet = dataset, mSet = mSet, 
               fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/")
    lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
                   fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/", 
                   topnum = 75)
    
    headgroupStat(dataSet = dataset, mSet = mSet, 
                  fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/headgroup/")
    FAchainStat(dataSet = dataset, mSet = mSet, 
                fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/FAchainVisual/", 
                plotInfo = "FA_info")
  }
}else if(analOpt == "all_together"){
  cat("All groups will be analyzed together\n")
  mSet <- MARpreproc(dataSet = dataSet)
  lipPCAPlot(dataSet = dataSet, mSet = mSet, 
             fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/")
  lipHeatmapPlot(dataSet = dataSet, mSet = mSet, 
                 fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/", 
                 topnum = 75)
  
  headgroupStat(dataSet = dataSet, mSet = mSet, 
                fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/headgroup/")
  FAchainStat(dataSet = dataSet, mSet = mSet, 
              fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/FAchainVisual/", 
              plotInfo = "FA_info")
}else{
  cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
  dataset <- prepDataSet(analOpt)
  
  mSet <- MARpreproc(dataSet = dataset)
  lipVolcanoPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/")
  lipPCAPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/")
  lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
                 fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/MARresults/", 
                 topnum = 75)
  
  headgroupStat(dataSet = dataset, mSet = mSet, 
                fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/headgroup/")
  FAchainStat(dataSet = dataset, mSet = mSet, 
              fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/FAchainVisual/", 
              plotInfo = "all_info")

}

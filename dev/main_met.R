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


prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}

##!!!Client options: can be "all_together"/"group_by_group"/[expr group]
analOpt <- "all_together"
##!!!!!DEV: delOddChainOpt should always be F in the Metabolites mode
dataSet <- readingLipidData(datafile = "./branch/benchmark/input/HANLipidMediator_imm.CSV", 
                            controlGrp = "", dataType = "Metabolites", delOddChainOpt = F,
                            lipField = "LipidMediator",
                            fileLoc = "./branch/benchmark/output/")

cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
dataset <- prepDataSet(analOpt)
mSet <- MARpreproc(dataSet = dataset)
##!!!!!DEV: showLipClass should always be F in the Metabolites mode
lipVolcanoPlot(dataSet = dataset, mSet = mSet, showLipClass = F, 
               fileLoc = "./testData/CerebrospinalFluid_multiomics/output/MARresults/")
lipPCAPlot(dataSet = dataset, mSet = mSet, 
           fileLoc = "./testData/CerebrospinalFluid_multiomics/output/MARresults/")
lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
               fileLoc = "./testData/CerebrospinalFluid_multiomics/output/MARresults/", 
               topnum = 75)

if(analOpt == "all_together"){
  cat("All groups will be analyzed together\n")
  mSet <- MARpreproc(dataSet = dataSet)
  lipPCAPlot(dataSet = dataSet, mSet = mSet, 
             fileLoc = "./branch/benchmark/output/MARresults/")
  lipHeatmapPlot(dataSet = dataSet, mSet = mSet, 
                 fileLoc = "./branch/benchmark/output/MARresults/", 
                 topnum = 75)
  
  headgroupStat(dataSet = dataSet, mSet = mSet, 
                fileLoc = "./testData/CerebrospinalFluid_multiomics/output/headgroup/")
  FAchainStat(dataSet = dataSet, mSet = mSet, 
              fileLoc = "./testData/CerebrospinalFluid_multiomics/output/FAchainVisual/", 
              plotInfo = "FA_info")
}

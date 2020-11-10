source("./branch/benchmark/code_extra/readingLipidData_add.R")
source("./branch/benchmark/code_extra/tidyFormat_HCC.R")
source("./utilityFunc/MARpreproc.R")
source("./branch/benchmark/code_extra/lipVolcanoPlot_add.R")
source("./utilityFunc/lipPCAPlot.R")
source("./utilityFunc/lipHeatmapPlot.R")
source("./branch/benchmark/code_extra/headgroupStat_add.R")
source("./branch/benchmark/code_extra/FAchainStat_add.R")
source("./utilityFunc/plottingPalettes.R")
source("./utilityFunc/statFAChains.R")

prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}
analOpt <- "Fasn_neg_ mice injected with PLASMIDS"
dataSet <- readingLipidData(datafile = "./branch/benchmark/input/lipid_tidy.CSV", 
                            controlGrp = "Fasn_pos_  mice injected with PLASMIDS", 
                            dataType = "HCC", delOddChainOpt = F, lipField = "LipidName",
                            fileLoc = "./branch/benchmark/output/", 
                            na.char = "")
dataset <- prepDataSet(analOpt)

mSet <- MARpreproc(dataSet = dataset)
lipVolcanoPlot(dataSet = dataset, mSet = mSet, 
               fileLoc = "./branch/benchmark/output/MARresults/")
lipPCAPlot(dataSet = dataset, mSet = mSet, 
           fileLoc = "./branch/benchmark/output/MARresults/")
lipHeatmapPlot(dataSet = dataset, mSet = mSet, 
               fileLoc = "./branch/benchmark/output/MARresults/", 
               topnum = 75)

headgroupStat(dataSet = dataset, mSet = mSet, 
              fileLoc = "./branch/benchmark/output/headgroup/")
FAchainStat(dataSet = dataset, mSet = mSet, 
            fileLoc = "./branch/benchmark/output/FAchainVisual/", 
            plotInfo = "FA_info")












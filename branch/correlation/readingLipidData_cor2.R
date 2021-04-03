readingLipidData_cor <- function(datafile, sampleList, dataType, delOddChainOpt, 
                             perc, normOpt){
  source("./utilityFunc/readingLipidData.R")
  source("./utilityFunc/MARpreproc.R")
  dataSet <- readingLipidData(datafile = datafile, sampleList = sampleList, 
                              dataType = dataType, delOddChainOpt = delOddChainOpt, 
                              cor = T)
  mSet <- MARpreproc(dataSet = dataSet, fileLoc = NA, perc = perc, 
                     normOpt = normOpt)
  return(mSet[["dataSet"]][["norm"]])
}
  
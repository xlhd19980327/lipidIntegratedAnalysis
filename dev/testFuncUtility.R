## Test readingLipidData
dataSet <- readingLipidData(datafile = "~/temp/Area_0_2020951653_Serum.CSV", 
                            controlGrp = "", dataType = "MS_DIAL", delOddChainOpt = T, 
                            fileLoc = "~/temp/")
## Test MARpreproc
mSet <- MARpreproc(dataSet = dataSet)

## Test FAchainStat
dataSet = dataSet
mSet = mSet
fileLoc = "~/temp/"
plotInfo = "FA_info"

write.csv(lipid_subclass_stat_output, 
          paste0(fileLoc, "Serum.csv"), 
          row.names = F) 

## Test statFAChains
lipid_subclass_tidyStat = stat_res
lipsample = lipsample
spe = "hsa"

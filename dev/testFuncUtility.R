## Test readingLipidData
dataSet <- readingLipidData(datafile = "./testData/otherlipidData/20201020_lipid.CSV", 
                            controlGrp = "", dataType = "MS_DIAL", delOddChainOpt = T, 
                            fileLoc = "~/temp/")
datafile = "./testData/cold_induced/input/lipid_tidy.CSV"
controlGrp = "COLD"
dataType = "LipidSearch"
delOddChainOpt = T
lipField = NA
fileLoc = "./testData/cold_induced/output/"

## Test MARpreproc
mSet <- MARpreproc(dataSet = dataSet)

## Test headgroupStat
dataSet = dataset
mSet = mSet

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

## Test readingRNAData
datafile = "./testData/cold_induced/input/rna.csv"
sampleList = "./testData/cold_induced/input/sampleList.CSV"
controlGrp = "RT"

datafile = "./testData/cold_induced/input/rna_genesymbol.CSV"
sampleList = "./testData/cold_induced/input/sampleList.CSV"
controlGrp = "RT"

datafile = "./testData/RNAseq_test/input/GSE148729_Caco2_polyA_readcounts.tsv"
sampleList = "./testData/RNAseq_test/input/fileSample_batch.csv"
controlGrp = ""

## Test rnaHeatmapPlot
DEAresult = dataProc_RNA
showtop = 75
showallgroups = T
fileLoc = "./testData/cold_induced/output/"

## Test DESeq2preproc
dataSet = dataSet_RNA

## Test statFAChains
lipid_subclass_tidyStat = stat_res
fileLoc = "./testData/cold_induced/output/"
lipsample = lipsample
spe = "mmu"



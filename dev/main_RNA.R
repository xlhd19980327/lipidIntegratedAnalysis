## Use this script to test function procedure
source("./utilityFunc/readingRNAData.R")
source("./utilityFunc/DESeq2preproc.R")
source("./utilityFunc/rnaPCAPlot.R")
source("./utilityFunc/DEAnalysis.R")
source("./utilityFunc/rnaVolcanoPlot.R")
source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")

dataSet_RNA <- readingRNAData(datafile = "./testData/RNAseq_test/input/GSE148729_Caco2_polyA_readcounts.tsv", 
                              sampleList = "./testData/RNAseq_test/input/fileSample_batch.csv", 
                              controlGrp = "untr")
dataProc_RNA <- DESeq2preproc(dataSet = dataSet_RNA, 
                              fileLoc = "./testData/RNAseq_test/output/")
rnaHeatmapPlot(DEAresult = dataProc_RNA, showallgroups = T, 
               fileLoc = "./testData/RNAseq_test/output/")
rnaPCAPlot(dataProc = dataProc_RNA, 
           fileLoc = "./testData/RNAseq_test/output/")
analOpt = "group_by_group"
if(analOpt == "group_by_group"){
  for(i in dataSet_RNA$groupsLevel[dataSet_RNA$groupsLevel != dataSet_RNA$controlGrp]){
    DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = i, 
                            fileLoc = "./testData/RNAseq_test/output/")
    rnaVolcanoPlot(DEAresult, 
                   fileLoc = "./testData/RNAseq_test/output/")
    rnaHeatmapPlot(DEAresult, 
                   fileLoc = "./testData/RNAseq_test/output/")
  }
}else{
  DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = analOpt, 
                          fileLoc = "./testData/RNAseq_test/output/")
  rnaVolcanoPlot(DEAresult, 
                 fileLoc = "./testData/RNAseq_test/output/")
  rnaHeatmapPlot(DEAresult, 
                 fileLoc = "./testData/RNAseq_test/output/")
}

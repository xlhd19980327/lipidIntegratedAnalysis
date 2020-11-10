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
                              controlGrp = "")
dataProc_RNA <- DESeq2preproc(dataSet = dataSet_RNA, 
                              fileLoc = "./testData/RNAseq_test/output/")
## These plots only show the result after normalization of size factor(account for differences in sequencing depth)
## We recommand you run this PCA/Heatmap only when you want to see within-group/between-group variability across multiple groups 
## If you want to see the detail DE genes, run Heatmap/Volcano between two groups in the later section
rnaPCAPlot(dataProc = dataProc_RNA, 
           fileLoc = "./testData/RNAseq_test/output/")
##!!!!!DEV: rnaHeatmapPlot input differ in different situation
rnaHeatmapPlot(DEAresult = dataProc_RNA, showallgroups = T, 
               fileLoc = "./testData/RNAseq_test/output/")
analOpt <- "mock"
if(analOpt == "group_by_group"){
  for(i in dataSet_RNA$groupsLevel[dataSet_RNA$groupsLevel != dataSet_RNA$controlGrp]){
    DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = i, 
                            fileLoc = "./testData/cold_induced/output/")
    rnaVolcanoPlot(DEAresult, 
                   fileLoc = "./testData/cold_induced/output/")
    rnaHeatmapPlot(DEAresult, 
                   fileLoc = "./testData/cold_induced/output/")
  }
}else{
  DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = analOpt, 
                          fileLoc = "./testData/RNAseq_test/output/")
  rnaVolcanoPlot(DEAresult, 
                 fileLoc = "./testData/RNAseq_test/output/")
  rnaHeatmapPlot(DEAresult, 
                 fileLoc = "./testData/RNAseq_test/output/")
}

source("./utilityFunc/readingRNAData.R")
source("./utilityFunc/limmaPreproc.R")
source("./utilityFunc/rnaMiArPCAPlot.R")
source("./utilityFunc/DEAnalysis_limma.R")
source("./utilityFunc/rnaVolcanoPlot.R")
source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")

dataSet_RNA <- readingRNAData(datafile = "./branch/benchmark/input/gene_tidy.CSV", 
                              sampleList = "./branch/benchmark/input/sampleList.CSV", 
                              controlGrp = "HCC_FASNpos")
dataProc_RNA <- limmaPreproc(dataSet = dataSet_RNA, norm = T, 
                             fileLoc = "./branch/benchmark/output/")
## These plots only show the result after normalization
## We recommand you run this PCA/Heatmap only when you want to see within-group/between-group variability across multiple groups 
## If you want to see the detail DE genes, run Heatmap/Volcano between two groups in the later section
rnaMiArPCAPlot(dataProc = dataProc_RNA, 
               fileLoc = "./branch/benchmark/output/")
##!!!!!DEV: rnaHeatmapPlot input differ in different situation
rnaHeatmapPlot(DEAresult = dataProc_RNA, showallgroups = T, type = "MiAr", 
               fileLoc = "./branch/benchmark/output/")
analOpt = "HCC_FASNneg"
if(analOpt == "group_by_group"){
  for(i in dataSet_RNA$groupsLevel[dataSet_RNA$groupsLevel != dataSet_RNA$controlGrp]){
    DEAresult <- DEAnalysis_limma(dataProc = dataProc_RNA, experGrp = i, 
                            fileLoc = "./testData/cold_induced/output/")
    rnaVolcanoPlot(DEAresult, 
                   fileLoc = "./testData/cold_induced/output/")
    rnaHeatmapPlot(DEAresult, 
                   fileLoc = "./testData/cold_induced/output/")
  }
}else{
  DEAresult <- DEAnalysis_limma(dataProc = dataProc_RNA, experGrp = analOpt, 
                          fileLoc = "./branch/benchmark/output/")
  rnaVolcanoPlot(DEAresult, type = "MiAr",
                 fileLoc = "./branch/benchmark/output/")
  rnaHeatmapPlot(DEAresult, type = "MiAr",
                 fileLoc = "./testData/cold_induced/output/")
}

## Use this script to test function procedure
source("./utilityFunc/readingRNAData.R")
source("./utilityFunc/DESeq2preproc.R")
source("./utilityFunc/rnaPCAPlot.R")
source("./utilityFunc/DEAnalysis.R")
source("./utilityFunc/rnaVolcanoPlot.R")
source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")

outputLoc <- "./branch/benchmark/output/"
dataSet_RNA <- readingRNAData(datafile = "./branch/benchmark/input/HANgene_tidy_geneid_allgroups.CSV", 
                              sampleList = "./branch/benchmark/input/HANsampleList_allgroups.CSV", 
                              controlGrp = "")
dataProc_RNA <- DESeq2preproc(dataSet = dataSet_RNA, 
                              fileLoc = outputLoc)

##!!!Client options: can be "group_by_group"/"only_show_variability"/[expr group]
analOpt <- "Ly6ChighD4"
if(analOpt == "group_by_group"){
  for(i in dataSet_RNA$groupsLevel[dataSet_RNA$groupsLevel != dataSet_RNA$controlGrp]){
    DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = i, 
                            fileLoc = outputLoc)
    rnaVolcanoPlot(DEAresult, 
                   fileLoc = outputLoc)
    rnaHeatmapPlot(DEAresult, 
                   fileLoc = outputLoc)
  }
}else if(analOpt == "only_show_variability"){
  ##NOTE: These plots only show the result after normalization of size factor(account for differences in sequencing depth)
  ##NOTE: We recommand you run this PCA/Heatmap only when you want to see within-group/between-group variability across multiple groups 
  ##NOTE: If you want to see the detail DE genes, run Heatmap/Volcano between two groups using "group_by_group"/[expr group] options
  rnaPCAPlot(dataProc = dataProc_RNA, 
             fileLoc = outputLoc)
  ##!!!!!DEV: rnaHeatmapPlot input differ in different situation
  rnaHeatmapPlot(DEAresult = dataProc_RNA, showallgroups = T, 
                 fileLoc = outputLoc)
}else{
  DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = analOpt, 
                          fileLoc = outputLoc)
  rnaVolcanoPlot(DEAresult, 
                 fileLoc = outputLoc)
  rnaHeatmapPlot(DEAresult, 
                 fileLoc = outputLoc)
}

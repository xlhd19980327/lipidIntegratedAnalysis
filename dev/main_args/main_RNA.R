## Use this script to test function procedure
source("./utilityFunc/readingRNAData.R")
source("./utilityFunc/DESeq2preproc.R")
source("./utilityFunc/rnaPCAPlot.R")
source("./utilityFunc/DEAnalysis.R")
source("./utilityFunc/rnaVolcanoPlot.R")
source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")

source("./utilityFunc/limmaPreproc.R")
source("./utilityFunc/rnaMiArPCAPlot.R")
source("./utilityFunc/DEAnalysis_limma.R")

library(RColorBrewer)

#RNA-seq
library(DESeq2)
library(limma)
library(apeglm)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(ggsci)
options(stringsAsFactors = F)

#Shell interaction
#library(getopt)
library(optparse)

option_list <- list( 
  make_option(c("-a", "--analysis_option"), action="store"),
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  make_option(c("-c", "--control_group"), action="store", default = ""), 
  make_option(c("-o", "--tidy_output"), action="store"), 
  make_option(c("-n", "--norm"), action="store", default = T), 
  make_option(c("-q", "--dea_putput"), action="store"), 
  make_option(c("-t", "--data_type"), action="store", default = "RNAseq"), 
  make_option(c("-s", "--volcano_output"), action="store"), 
  make_option(c("-f", "--fc_thresh"), action="store", default = 2.0), 
  make_option(c("-p", "--p_thresh"), action="store", default = 0.1), 
  make_option(c("-u", "--show_volcano_top"), action="store", default = 20), 
  make_option(c("-v", "--show_heatmap_top"), action="store", default = 75), 
  make_option(c("-w", "--heatmap_output"), action="store"), 
  make_option(c("-r", "--pca_output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dataSet_RNA <- readingRNAData(datafile = opt$input_file, 
                              sampleList = opt$description_file, 
                              controlGrp = opt$control_group)
if(opt$data_type == "RNAseq"){
  dataProc_RNA <- DESeq2preproc(dataSet = dataSet_RNA, 
                              fileLoc = opt$tidy_output)
}
if(opt$data_type == "MiAr"){
  dataProc_RNA <- limmaPreproc(dataSet = dataSet_RNA, norm = opt$norm, 
                               fileLoc = opt$tidy_output)
}

##!!!Client options: can be "group_by_group"/"only_show_variability"/[expr group]
analOpt <- opt$analysis_option
if(analOpt == "group_by_group"){
  for(i in dataSet_RNA$groupsLevel[dataSet_RNA$groupsLevel != dataSet_RNA$controlGrp]){
    if(opt$data_type == "RNAseq"){
      DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = i, 
                            fileLoc = opt$dea_putput)
    }
    if(opt$data_type == "MiAr"){
      DEAresult <- DEAnalysis_limma(dataProc = dataProc_RNA, experGrp = i, 
                                    fileLoc = opt$dea_putput)
    }
    rnaVolcanoPlot(DEAresult = DEAresult, type = opt$data_type, fcthresh = opt$fc_thresh, 
                   pthresh = opt$p_thresh, showtop = opt$show_volcano_top, 
                   fileLoc = opt$volcano_output)
    rnaHeatmapPlot(DEAresult, showtop = opt$show_heatmap_top, 
                   type = opt$data_type, fileLoc = opt$heatmap_output)
  }
}else if(analOpt == "only_show_variability"){
  ##NOTE: These plots only show the result after normalization of size factor(account for differences in sequencing depth)
  ##NOTE: We recommand you run this PCA/Heatmap only when you want to see within-group/between-group variability across multiple groups 
  ##NOTE: If you want to see the detail DE genes, run Heatmap/Volcano between two groups using "group_by_group"/[expr group] options
  if(opt$data_type == "RNAseq"){
    rnaPCAPlot(dataProc = dataProc_RNA, 
             fileLoc = opt$pca_output)
  }
  if(opt$data_type == "MiAr"){
    rnaMiArPCAPlot(dataProc = dataProc_RNA, 
                   fileLoc = opt$pca_output)
  }
  ##!!!!!DEV: rnaHeatmapPlot input differ in different situation
  rnaHeatmapPlot(DEAresult = dataProc_RNA, showallgroups = T, type = opt$data_type, 
                 fileLoc = opt$heatmap_output)
}else{
  if(opt$data_type == "RNAseq"){
    DEAresult <- DEAnalysis(dataProc = dataProc_RNA, experGrp = analOpt, 
                          fileLoc = opt$dea_putput)
  }
  if(opt$data_type == "MiAr"){
    DEAresult <- DEAnalysis_limma(dataProc = dataProc_RNA, experGrp = analOpt, 
                                  fileLoc = opt$dea_putput)
  }
  rnaVolcanoPlot(DEAresult = DEAresult, type = opt$data_type, fcthresh = opt$fc_thresh, 
                 pthresh = opt$p_thresh, showtop = opt$show_volcano_top, 
                 fileLoc = opt$volcano_output)
  rnaHeatmapPlot(DEAresult, showtop = opt$show_heatmap_top, 
                 type = opt$data_type, fileLoc = opt$heatmap_output)
}

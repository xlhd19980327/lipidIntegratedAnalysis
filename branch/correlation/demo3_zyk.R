## Do correlation
library(tidyverse)
library(dplyr)
library(pheatmap)
options(stringsAsFactors = F)
source("./branch/correlation/readingLipidData_cor.R")
source("./branch/correlation/readingRNAData_cor.R")
kg <- 10
kl <- 4

## Data input
lipiddataSet <- readingLipidData(datafile = "./branch/benchmark/input/HANLipidMediator_imm_forcor.CSV", 
                                 controlGrp = "", dataType = "Metabolites", delOddChainOpt = F,
                                 lipField = "LipidMediator", sampleList = "./branch/benchmark/input/HANsampleList_lipmid.csv", 
                                 fileLoc = "./branch/benchmark/output/")
genedataSet <- readingRNAData(datafile = "./branch/benchmark/input/HANgene_tidy.CSV", 
                              sampleList = "./branch/benchmark/input/HANsampleList.CSV")

## Tidy the input files
alllipidnames <- lipiddataSet$lipidName
lipid_data <- as.data.frame(lipiddataSet$data)
rownames(lipid_data) <- paste0("L", 1:nrow(lipid_data)) 
allgenenames <- rownames(genedataSet$data)
gene_data <- as.data.frame(genedataSet$data)
rownames(gene_data) <- paste0("G", 1:nrow(gene_data))

##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

### Waiting time(~15+ min)
correlation <- cor(t(lipid_data), t(gene_data), use='pairwise.complete.obs', method ="spearman")
correlation[is.na(correlation)] <- 0
correlation <- correlation[,which(colSums(abs(correlation))>0)]
correlation <- correlation[which(rowSums(abs(correlation))>0),]
##rearrange lipid and gene ~5min
list <- pheatmap(correlation)
lipidrank <- list[["tree_row"]][["order"]]
generank <- list[["tree_col"]][["order"]]
correlation <- correlation[,rank(generank)]
correlation <- correlation[rank(lipidrank),]

##Figure out the size of submatrix
rownames(correlation) <- alllipidnames[as.integer(gsub("^L([0-9]+)", "\\1", 
                                                       rownames(correlation)))]
colnames(correlation) <- allgenenames[as.integer(gsub("^G([0-9]+)", "\\1", 
                                                      colnames(correlation)))]
list <- pheatmap(correlation)
gene_cluster_num <- cutree(list$tree_col, k=kg)
lipid_cluster_num <- cutree(list$tree_row, k=kl)
gene_cluster <- split(as.data.frame(t(correlation)), gene_cluster_num)
lipid_cluster <- split(as.data.frame(correlation), lipid_cluster_num)

##extract submatrix
correlation_data <- list()
gene_lipid <- list()
initial <- 1
for(i in 1:length(gene_cluster)){
  gene_size <- dim(gene_cluster[[i]][1])
  sub1 <- as.data.frame(correlation[,initial:(initial-1+gene_size)])
  sub_initial <- 1
  for(j in 1:length(lipid_cluster)){
    lipid_size <- dim(lipid_cluster[[j]][1])
    sub2 <- as.data.frame(sub1[sub_initial:(sub_initial-1+lipid_size),])
    name <- paste0(j, sep = ',', i)
    correlation_data[[name]] <- sub2
    lipid_name <- rownames(sub2)
    gene_name <- colnames(sub2)
    label <- list(lipid_name, gene_name)
    gene_lipid[[name]] <- label
    sub_initial <- sub_initial+lipid_size
 }
  initial <- initial+gene_size
}












## Do correlation
library(tidyverse)
library(dplyr)
library(pheatmap)
options(stringsAsFactors = F)
source("./branch/correlation/readingRNAData_cor.R")
source("./branch/correlation/readingLipidData_cor.R")


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

##rearrange lipid and gene ~5min
correlation <- correlation[,which(colSums(abs(correlation))>0)]
correlation <- correlation[which(rowSums(abs(correlation))>0),]
list <- pheatmap(correlation)
lipidrank <- list[["tree_row"]][["order"]]
generank <- list[["tree_col"]][["order"]]
correlation <- correlation[,rank(generank)]
correlation <- correlation[rank(lipidrank),]

##Figure out the size of submatrix
gene_length <- dim(correlation)[2]
lipid_length <- dim(correlation)[1]
gene_size <- ceiling(gene_length/100)
lipid_size <- ceiling(lipid_length/100)

##extract submatrix
correlation_data <- list()
gene_lipid <- list()
for(i in 1:gene_length){
  if((i-1+gene_size) <= gene_length){
    sub1 <- as.data.frame(correlation[,i:(i-1+gene_size)])
  }
  else{
    sub1 <- as.data.frame(correlation[,i:gene_length])
  }
  for(j in 1:lipid_length){
    if((j-1+lipid_size) <= lipid_length){
      sub2 <- as.data.frame(sub1[j:(j-1+lipid_size),])
    }
    else{
      sub2 <- as.data.frame(sub1[j:lipid_length,])
    }
    name <- paste0(j, sep = ',', i)
    correlation_data[[name]] <- sub2
    lipid_name <- rownames(sub2)
    gene_name <- colnames(sub2)
    label <- list(lipid_name, gene_name)
    gene_lipid[[name]] <- label
    j <- j+lipid_size
 }
  i <- i+gene_size
}












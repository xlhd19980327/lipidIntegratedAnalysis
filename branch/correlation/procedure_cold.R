library(tidyverse)
library(pheatmap)
options(stringsAsFactors = F)
source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/readingRNAData.R")

h_gene=18
h_lipid=10

## Data input
lipiddataSet <- readingLipidData(datafile = "./testData/cold_induced/input/lipid_tidy_forcor.CSV", 
                                 dataType = "LipidSearch", delOddChainOpt = T,
                                 fileLoc = NA, 
                                 na.char = "######")
genedataSet <- readingRNAData(datafile = "./testData/cold_induced/input/rna_genesymbol.CSV", 
                              sampleList = "./testData/cold_induced/input/sampleList.CSV")

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

# cluster analysis
sort_function <- function(x){
  n <- nrow(x)
  l <- max(round(0.03*n), 1)
  x[(1:l),]
}

# ——gene
lipidsum <- colSums(abs(correlation), na.rm = T, dims = 1)
gene_rank <- order(lipidsum, decreasing = T)
correlation <- correlation[, gene_rank]
### Waiting time(~??15+ min, rupt)
list <- pheatmap(correlation)
gene_cluster <- cutree(list$tree_col,h=h_gene)
gene_group_list <- split(data.frame(t(correlation)), gene_cluster)
group_gene_sort_result<- lapply(gene_group_list, sort_function)
genefinal <- do.call(rbind, group_gene_sort_result)
genenames <- lapply(group_gene_sort_result, rownames)
genenames <- do.call(c, genenames)
rownames(genefinal) <- genenames

#——lipid
genesum <- rowSums(abs(genefinal), na.rm = T, dims = 1)
lipid_rank <- order(genesum, decreasing = T)
genefinal <- genefinal[lipid_rank, ]
### Waiting time(Very Fast)
list2 <- pheatmap(genefinal)
lipid_cluster <- cutree(list2$tree_col,h=h_lipid)
lipid_group_list <- split(data.frame(t(genefinal)), lipid_cluster)
group_lipid_sort_result<- lapply(lipid_group_list, sort_function)
allfinal <- do.call(rbind, group_lipid_sort_result)
lipidnames <- lapply(group_lipid_sort_result, rownames)
lipidnames <- do.call(c, lipidnames)
rownames(allfinal) <- lipidnames

#Adjust names of the lipids and genes
finalgenenames <- allgenenames[as.integer(gsub("^G([0-9]+)", "\\1", colnames(allfinal)))]
colnames(allfinal) <- finalgenenames
finallipidnames <- alllipidnames[as.integer(gsub("^L([0-9]+)", "\\1", rownames(allfinal)))]
rownames(allfinal) <- finallipidnames

pdf(file = 'refine.pdf',width = 50, height = 50)
pheatmap(allfinal)
dev.off()
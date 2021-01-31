## Do correlation
library(pheatmap)
library(tidyverse)
library(MetaboAnalystR)
library(DESeq2)
library(limma)
options(stringsAsFactors = F)
source("./branch/correlation/readingLipidData_cor2.R")
source("./branch/correlation/readingRNAData_cor.R")
kg <- 6
kl <- 7

## Data input
#lipid_data <- readingLipidData(datafile = "./testData/SVFmultiomics_210118/input/lipids.csv", 
#                               controlGrp = "", dataType = "LipidSearch", delOddChainOpt = T, na.char = "",
#                               lipField = "LipidIon", sampleList = "./testData/SVFmultiomics_210118/input/sampleList.csv")
lipid_data <- readingLipidData(datafile = "./testData/SVFmultiomics_210118/input/lipids.csv", 
                               dataType = "Lipids", delOddChainOpt = T, perc = 2/3*100, 
                               sampleList = "./testData/SVFmultiomics_210118/input/sampleList.csv")
gene_data <- readingRNAData(datafile = "./testData/SVFmultiomics_210118/input/RNAseq_genesymbol.csv", 
                              sampleList = "./testData/SVFmultiomics_210118/input/sampleList.csv",
                              type = "RNAseq")

##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

### Waiting time(~2+ min, larger data may require more)
correlation <- cor(lipid_data, t(gene_data), method ="spearman")
correlation[is.na(correlation)] <- 0
#correlation <- correlation[,which(apply(correlation, 2, FUN = sd) != 0)]
#correlation <- correlation[which(apply(correlation, 1, FUN = sd) != 0),]
cutoff <- function(x){
  x<-abs(x)
  x<-max(x)
}

#70% quantile
max_list <- apply(correlation,2,max)
value <- quantile(max_list,0.7)
#cor thresh
value <- 0.8

correlation <- correlation[,which(apply(correlation, 2, cutoff)>value)]## set a parameter for cutoff, default =0.8
#correlation <- correlation[which(apply(correlation, 1, cutoff)>0.8),]## set a parameter for cutoff

##Figure out the size of submatrix
result_lipid <- dist(correlation, method = "euclidean")
list_lipid <- hclust(d = result_lipid)
result_gene <- dist(t(correlation), method = "euclidean")
list_gene <- hclust(d = result_gene)
gene_cluster_num <- cutree(list_gene, k=kg)
lipid_cluster_num <- cutree(list_lipid, k=kl)
#gene_cluster <- split(as.data.frame(t(correlation)), gene_cluster_num)
#lipid_cluster <- split(as.data.frame(correlation), lipid_cluster_num)
###randomly select gene and lipid from each cluster
#selected_gene <- data.frame()
#for(i in 1:kg){
#  n <- dim(gene_cluster[[i]])[1]
#  selected <- sample(n, ceiling(0.1*n))
#  selected_cluster_gene <- gene_cluster[[i]][selected,]
#  selected_gene <- rbind(selected_gene, selected_cluster_gene)
#}
#result_lipid.1 <- as.matrix(result_lipid)
#result_gene.1 <- as.matrix(result_gene)
#result_lipid.2 <- result_lipid.1[,colnames(selected_gene)]
#result_lipid.2 <- result_lipid.2[colnames(selected_gene),]
#result_lipid.dist <- as.dist(result_lipid.2)
#result_gene.2 <- result_gene.1[,rownames(selected_gene)]
#result_gene.2 <- result_gene.2[rownames(selected_gene),]
#result_gene.dist <- as.dist(result_gene.2)
#refine_pheatmap <- pheatmap(t(selected_gene), 
#                            clustering_distance_rows = result_lipid.dist, clustering_distance_cols = result_gene.dist, 
#                            cutree_col = kg, cutree_rows = kl)

### Waiting time(~2+ min, larger data may require more)
#png("~/temp/cor/correlationPlot.png", 
#    width = 720, height = 720)
pdf("~/temp/cor/correlationPlot.pdf")
list <- pheatmap(correlation,
                 cutree_col = kg, cutree_rows = kl)
dev.off()

cutgene <- cutree(list$tree_col, k = kg)[list[["tree_col"]][["order"]]]
cutlipid <- cutree(list$tree_row, k = kl)[list[["tree_row"]][["order"]]]
genegaps <- c(names(cutgene[1]), 
               names(which((cutgene[-1] - cutgene[-length(cutgene)]) != 0)))
lipidgaps <- c(names(cutlipid[1]), 
               names(which((cutlipid[-1] - cutlipid[-length(cutlipid)]) != 0)))
for(i in 1:kl){
  for(j in 1:kg){
    subdata <- correlation[names(cutlipid)[cutlipid == i], names(cutgene)[cutgene == j], drop = F]
    ith <- which(lipidgaps %in% rownames(subdata))
    jth <- which(genegaps %in% colnames(subdata))
    sublipids <- rownames(subdata)
    subgenes <- colnames(subdata)
    write.csv(subdata, paste0("~/temp/cor_sub/", "correlation_", ith, "_", jth, ".csv"))
    write.csv(sublipids, paste0("~/temp/cor_sub/", "lipids_", ith, ".csv"), row.names = F)
    write.csv(subgenes, paste0("~/temp/cor_sub/", "genes_", jth, ".csv"), row.names = F)
  }
}


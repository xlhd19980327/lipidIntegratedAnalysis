## Do correlation
library(pheatmap)
library(tidyverse)
library(MetaboAnalystR)
library(DESeq2)
library(limma)
library(dbscan)
library(MCL)
library(scales)
library(reshape)
library(ggplot2)

options(stringsAsFactors = F)
source("./branch/correlation/readingLipidData_cor.R")
source("./branch/correlation/readingRNAData_cor.R")
kg <- 3
kl <- 5
MinPts.l <- 5
MinPts.g <- 5

## Data input
lipid_data <- readingLipidData(datafile = "./testData/SVFmultiomics_210118/input/lipids.csv", 
                               controlGrp = "", dataType = "LipidSearch", delOddChainOpt = T, na.char = "",
                               lipField = "lipidName", sampleList = "./testData/SVFmultiomics_210118/input/sampleList.csv")
gene_data <- readingRNAData(datafile = "./testData/SVFmultiomics_210118/input/RNAseq_id.csv", 
                              sampleList = "./testData/SVFmultiomics_210118/input/sampleList.csv",
                              type = "RNAseq")

##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

### Waiting time(~2+ min, larger data may require more)
cutoff <- function(x){
  x<-abs(x)
  x<-max(x)
}
correlation <- cor(lipid_data, t(gene_data), method ="spearman")
correlation[is.na(correlation)] <- 0
correlation <- correlation[,which(apply(correlation, 2, FUN = sd) != 0)]
correlation <- correlation[which(apply(correlation, 1, FUN = sd) != 0),]
max_list <- apply(correlation,2,max)
value <- quantile(max_list,0.7)

correlation <- correlation[,which(apply(correlation, 2, cutoff)>value)]

##hierarchical clustering data preparation
result_lipid <- dist(correlation, method = "euclidean")
list_lipid <- hclust(d = result_lipid)
result_gene <- dist(t(correlation), method = "euclidean")
list_gene <- hclust(d = result_gene)
gene_cluster_num <- cutree(list_gene, k=kg)
lipid_cluster_num <- cutree(list_lipid, k=kl)
list <- pheatmap(t(correlation), fontsize_row = 4, clustering_distance_rows = result_gene, clustering_distance_cols = result_lipid,
                 cutree_col = kl, cutree_rows = kg)
#k-means clustering data preparation
kmeans.lipid.cluster <- kmeans(correlation, kl)
kmeans.gene.cluster <- kmeans(t(correlation), kg)

cluster.l <- data.frame(kmeans.lipid.cluster[["cluster"]])
cluster.l[,'lipid'] <- rownames(cluster.l)
cluster.l <- cluster.l[order(cluster.l$kmeans.lipid.cluster...cluster...),]
colnames(cluster.l) <- c('cluster.l', 'lipid')

cluster.g <- data.frame(kmeans.gene.cluster[["cluster"]])
cluster.g[,'gene'] <- rownames(cluster.g)
cluster.g <- cluster.g[order(cluster.g$kmeans.gene.cluster...cluster...),]
colnames(cluster.g) <- c('cluster.g', 'gene')

mat <- correlation[order(kmeans.lipid.cluster[["cluster"]]),order(kmeans.gene.cluster[["cluster"]])]
mat.m <- melt(mat)
colnames(mat.m) <- c('lipid', 'gene', 'value')
mat.m <- merge(mat.m, cluster.g, by = 'gene')
mat.m <- merge(mat.m, cluster.l, by = 'lipid')

p=ggplot(mat.m,aes(gene, lipid))+geom_tile(aes(fill=value))
p=p+scale_fill_gradientn(colours = c('#03A9F4','#81D4FA','#FFF9C4','#FFAB91','#FF5722'))
p=p+theme(text=element_text(size=4),axis.text.x=element_text(angle=90,vjust=0))
pic=p+facet_grid(cluster.l~cluster.g, scales='free', space='free')
pic

lipid_group <- split(cluster.l, cluster.l$cluster.l)
gene_group <- split(cluster.g, cluster.g$cluster.g)

#t <- pheatmap(mat, kmeans_k = 10, cutree_rows = 5, cutree_cols = 3)

#DBSCAN
Min <- min(correlation)
Max <- max(correlation)
Rank <- Max - Min
correlation <- (correlation-Min)/Rank

dbscan.lipid.cluster <- hdbscan(correlation, minPts = 2)
dbscan.gene.cluster <- hdbscan(t(correlation), minPts =2)

cluster.l <- data.frame(dbscan.lipid.cluster[["cluster"]])
cluster.l[,'lipid'] <- rownames(correlation)
cluster.l <- cluster.l[order(cluster.l$dbscan.lipid.cluster...cluster...),]
colnames(cluster.l) <- c('cluster.l', 'lipid')

cluster.g <- data.frame(dbscan.gene.cluster[["cluster"]])
cluster.g[,'gene'] <- colnames(correlation)
cluster.g <- cluster.g[order(cluster.g$dbscan.gene.cluster...cluster...),]
colnames(cluster.g) <- c('cluster.g', 'gene')

mat <- correlation[order(cluster.l$cluster.l),order(cluster.g$cluster.g)]
mat.m <- melt(mat)
colnames(mat.m) <- c('lipid', 'gene', 'value')
mat.m <- merge(mat.m, cluster.g, by = 'gene')
mat.m <- merge(mat.m, cluster.l, by = 'lipid')

p=ggplot(mat.m,aes(gene, lipid))+geom_tile(aes(fill=value))
p=p+scale_fill_gradientn(colours = c('#03A9F4','#81D4FA','#FFF9C4','#FFAB91','#FF5722'))
p=p+theme(text=element_text(size=4),axis.text.x=element_text(angle=90,vjust=0))
pic=p+facet_grid(cluster.l~cluster.g, scales='free', space='free')
pic

lipid_group <- split(cluster.l, cluster.l$cluster.l)
gene_group <- split(cluster.g, cluster.g$cluster.g)

#Markov

lipid_cor <- cor(t(correlation))
lipid_cor <- (lipid_cor - min(lipid_cor))/(max(lipid_cor) - min(lipid_cor))
lipid_cut <- quantile(lipid_cor,0.55)
lipid_cor[lipid_cor<=lipid_cut] <- 0

gene_cor <- cor(correlation)
gene_cor <- (gene_cor - min(gene_cor))/(max(gene_cor) - min(gene_cor))
gene_cut <- quantile(gene_cor,0.5)
gene_cor[gene_cor<=gene_cut] <- 0
mcl.lipid.cluster <- mcl(lipid_cor, addLoops= T, allow1 = T)
mcl.gene.cluster <- mcl(gene_cor, addLoops = T, allow1 = T)

cluster.l <- data.frame(mcl.lipid.cluster[["Cluster"]])
cluster.l[,'lipid'] <- rownames(correlation)
colnames(cluster.l) <- c('cluster.l', 'lipid')
cluster.l[,'raw.order'] <- as.numeric(rownames(cluster.l))
cluster.l <- arrange(cluster.l, cluster.l)
cluster.l[,'new.order'] <- order(cluster.l$cluster.l)
cluster.l <- arrange(cluster.l, raw.order)

cluster.g <- data.frame(mcl.gene.cluster[["Cluster"]])
cluster.g[,'gene'] <- colnames(correlation)
colnames(cluster.g) <- c('cluster.g', 'gene')
cluster.g[,'raw.order'] <- as.numeric(rownames(cluster.g))
cluster.g <- arrange(cluster.g, cluster.g)
cluster.g[,'new.order'] <- order(cluster.g$cluster.g)
cluster.g <- arrange(cluster.g, raw.order)

mat <- correlation[order(cluster.l$new.order),order(cluster.g$cluster.g)]

mat.m <- melt(mat)
colnames(mat.m) <- c('lipid', 'gene', 'value')
mat.m <- merge(mat.m, cluster.g, by = 'gene')
mat.m <- merge(mat.m, cluster.l, by = 'lipid')

p=ggplot(mat.m,aes(gene, lipid))+geom_tile(aes(fill=value))
p=p+scale_fill_gradientn(colours = c('#03A9F4','#81D4FA','#FFF9C4','#FFAB91','#FF5722'))
p=p+theme(text=element_text(size=4),axis.text.x=element_text(angle=90,vjust=0))
pic=p+facet_grid(cluster.l~cluster.g, scales='free', space='free')
pic

### Waiting time(~2min 30s, larger data may require more)
png("./correlationPlot.png", 
    width = 720, height = 720)
pic
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
    write.csv(subdata, paste0("~/temp/cor2/", "correlation_", ith, "_", jth, ".csv"))
    write.csv(sublipids, paste0("~/temp/cor2/", "lipids_", ith, ".csv"), row.names = F)
    write.csv(subgenes, paste0("~/temp/cor2/", "genes_", jth, ".csv"), row.names = F)
  }
}


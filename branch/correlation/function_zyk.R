heatmap_analysis <- function(workdirectory, gene, lipid, h_gene=18, h_lipid=10){
  library(tidyverse)
  library(pheatmap)
  library(dplyr)
  
  setwd(workdirectory)
  lipid <- read.csv(lipid, header = TRUE, na.strings = NA)
  gene <- read.table(gene, header = TRUE, na.strings = NA)
  gene <- gene[,-(2:3)]
  t1 <- t(data.frame(gene, row.names = 1))
  gene <- as.data.frame(t1,row.names=1)
  data.frame(lipid, row.names = 1)
  t2 <- t(data.frame(lipid,row.names=1))
  lipid <- as.data.frame(t2,row.names=F)
  lipid <- lipid[,-1]
  
  correlation <- cor(lipid, gene, use='pairwise.complete.obs', method ="spearman")
  #write.csv(correlation, file='correlation.csv')
  cor <- abs(correlation)
  colsum <- colSums(cor, na.rm = T, dims = 1)
  data <- data.frame(t(rbind(correlation,colsum)))
  correlation[is.na(correlation)] <- 0
  #pdf(file = 'heatmap.pdf',width = 500, height = 350)
  #pheatmap(correlation,treeheight_row = 10000, treeheight_col = 10000, row.names = 1)
  #dev.off()
  
  # cluster analysis
  sort_function <- function(x){
    n <- ncol(x)
    l <- round(0.03*n)
    x[(1:l),]
  }
  list <- pheatmap(correlation)
  # ——gene
  lipidsum <- colSums(abs(correlation), na.rm = T, dims = 1)
  # 获得分群结果
  gene_cluster <- cutree(list$tree_col,h=h_gene)
  #gene_cluster <- cutree(list$tree_col,h=18)
  gene_group <- t(data.frame(gene_cluster))
  #将基因分群结果与基因与脂质相关性、总相关性结合
  grouped_gene <- data.frame(t(rbind(correlation, lipidsum, gene_group)))
  #根据总相关性排序
  gene_rank <- order(grouped_gene$lipidsum, decreasing = T)
  ranked_grouped_gene <- grouped_gene[gene_rank,]
  #按照cluster划分数据
  gene_group_list <- split(ranked_grouped_gene,ranked_grouped_gene$gene_cluster)
  #对每一个cluster应用sortfunction
  group_gene_sort_result<- lapply(gene_group_list,sort_function)
  #将每个cluster里sort出的数据汇总
  genefinal <- do.call(rbind, group_gene_sort_result)
  n <- ncol(genefinal)
  genefinal <- genefinal[,-((n-1):n)]
  genefinal <- genefinal[complete.cases(genefinal),]
  
  #——lipid
  list2 <- pheatmap(genefinal)
  genesum <- colSums(abs(genefinal), na.rm = T, dims = 1)
  lipid_cluster <- cutree(list2$tree_col, h = h_lipid)
  #lipid_cluster <- cutree(list2$tree_col, h = 10)
  lipid_group <- t(data.frame(lipid_cluster))
  genesum <- t(data.frame(genesum))
  grouped_lipid <- data.frame(t(rbind(genefinal,lipid_group, genesum)))
  lipid_rank <- order(grouped_lipid$genesum, decreasing = T)
  ranked_grouped_lipid <- grouped_lipid[lipid_rank,]
  lipid_group_list <- split(ranked_grouped_lipid,ranked_grouped_lipid$lipid_cluster)
  grouped_lipid_sort_result <- lapply(lipid_group_list, sort_function)
  lipidfinal <- do.call(rbind, grouped_lipid_sort_result)
  lipidfinal <- lipidfinal[complete.cases(lipidfinal),]
  n <- ncol(lipidfinal)
  lipidfinal <- lipidfinal[,-((n-1):n)]
  
  pdf(file = 'refine.pdf',width = 50, height = 50)
  pheatmap(lipidfinal,treeheight_row = 800, treeheight_col = 800)
  dev.off()
}


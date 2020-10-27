library(tidyverse)
library(pheatmap)
options(stringsAsFactors = F)

input <- "./branch/correlation/input/"
output <- "./branch/correlation/output/"
gene = "200216MultiOmicsSVF_gene_expression.txt"
lipid = "lipiddata.csv"
h_gene=3
h_lipid=3

lipid <- read.csv(paste0(input, lipid), header = TRUE, na.strings = "NA")
gene <- read.table(paste0(input, gene), header = TRUE, na.strings = "NA")

## Tidy the input files
gene <- gene[,-(2:3)]
lipid_sub <- lipid[1:100, ]
gene_sub <- gene[1:500, ]
alllipidnames <- lipid_sub$lipidName
lipid_sub <- select(lipid_sub, -lipidName)
rownames(lipid_sub) <- paste0("L", 1:nrow(lipid_sub)) 
allgenenames <- gene_sub$st_gene_id
gene_sub <- select(gene_sub, -st_gene_id)
rownames(gene_sub) <- paste0("G", 1:nrow(gene_sub))

### Waiting time(~10 min)
correlation <- cor(t(lipid_sub), t(gene_sub), use='pairwise.complete.obs', method ="spearman")
#correlation[is.na(correlation)] <- 0

# cluster analysis
sort_function <- function(x){
  n <- nrow(x)
  l <- max(round(0.03*n), 1)
  x[(1:l),]
}

# ——gene
lipidsum <- colSums(abs(correlation), na.rm = T, dims = 1)
#根据总相关性排序
gene_rank <- order(lipidsum, decreasing = T)
correlation <- correlation[, gene_rank]
### Waiting time(~20 min)
list <- pheatmap(correlation)
# 获得分群结果
gene_cluster <- cutree(list$tree_col,h=3)

#按照cluster划分数据
gene_group_list <- split(data.frame(t(correlation)), gene_cluster)
#对每一个cluster应用sortfunction
group_gene_sort_result<- lapply(gene_group_list, sort_function)
#将每个cluster里sort出的数据汇总
genefinal <- do.call(rbind, group_gene_sort_result)
genenames <- sapply(group_gene_sort_result, rownames)
rownames(genefinal) <- genenames

#——lipid
#
genesum <- rowSums(abs(genefinal), na.rm = T, dims = 1)
#根据总相关性排序
lipid_rank <- order(genesum, decreasing = T)
genefinal <- genefinal[lipid_rank, ]
### Waiting time(~20 min)
list2 <- pheatmap(genefinal)
# 获得分群结果
lipid_cluster <- cutree(list2$tree_col,h=3)

#按照cluster划分数据
lipid_group_list <- split(data.frame(t(genefinal)), lipid_cluster)
#对每一个cluster应用sortfunction
group_lipid_sort_result<- lapply(lipid_group_list, sort_function)
#将每个cluster里sort出的数据汇总
allfinal <- do.call(rbind, group_lipid_sort_result)
lipidnames <- sapply(group_lipid_sort_result, rownames)
rownames(allfinal) <- lipidnames

#Adjust names of the lipids and genes
finalgenenames <- allgenenames[as.integer(gsub("^G([0-9]+)", "\\1", colnames(allfinal)))]
colnames(allfinal) <- finalgenenames
finallipidnames <- alllipidnames[as.integer(gsub("^L([0-9]+)", "\\1", rownames(allfinal)))]
rownames(allfinal) <- finallipidnames

pdf(file = 'refine.pdf',width = 50, height = 50)
pheatmap(allfinal)
dev.off()
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(pheatmap)

# 读入counts数据
data <- read.table("/Users/LCX/Desktop/RNA-seq_fw/quantification_counts.txt",
                  header=T, row.names=1, com='', 
                  quote='',check.names=F, sep="")
head(data)


# 撇掉在多于两个样本中count<1的值,排除极低表达基因的干扰
data <- data[rowSums(data)>2,]
head(data)


# 读入分组信息
# sample一般分为两列，第一列为sample名字，
# 第二列为conditions（比如treatment，control）
sample <- read.table("/Users/LCX/Desktop/RNA-seq_fw/SampleGroup.tsv",
                     header=T, row.names=1, com='',quote='', check.names=F,
                     sep="\t", colClasses="factor")
sample <- sample[match(colnames(data), rownames(sample)),drop=F]                     
sample


# 产生DESeq数据集
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                            colData = sample, 
                                            design= ~ conditions)
dds <- DESeq(ddsFullCountTable)
dds


# 获取标准化后的数据
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)


# log转换后的结果
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]




# PCA分析
pca_data <- plotPCA(rld, intgroup=c("conditions"), returnData=T)
percentVar <- round(100*attr(pca_data, "percentVar")) 
ggplot(pca_data, aes(PC1, PC2, color=conditions)) + 
  geom_point(size=3) + 
  ggtitle("PCA") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))




# 差异基因分析
# 定义变量
sampleA = "Treatment"
sampleB = "Control"

# constrast用于指定比较的两组的信息
# 输出的log2FoldChange为log2(SampleA/SampleB)
contrastV <- c("conditions", sampleA, sampleB)
res <- results(dds,  contrast=contrastV)
res


# 获得第一组数据均值
baseA <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleA]
if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
head(baseMeanA)

# 获得第二组数据均值
baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
head(baseMeanB)

# 结果组合
res <- cbind(baseMeanB, baseMeanA, as.data.frame(res))
head(res)

# 增加ID信息
res <- cbind(ID=rownames(res), as.data.frame(res))
res$baseMean <- rowMeans(cbind(baseB, baseA))

# 校正后p-value为NA的复制为1
res$padj[is.na(res$padj)] <- 1
res$pvalue[is.na(res$pvalue)] <- 1

# 按pvalue排序, 把差异大的基因放前面
res <- res[order(res$pvalue),]
head(res)


# 差异基因筛选，padj<0.1
res_de <- subset(res, res$padj<0.1, 
                 select=c('ID', sampleB,
                          sampleA, 'log2FoldChange', 'padj'))
# foldchang > 1
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)

# 差异基因ID
res_de_up_id = data.frame(ID=res_de_up$ID, 
                          type=paste(sampleB,"_higherThan_", 
                                     sampleA, sep="."))
res_de_dw_id = data.frame(ID=res_de_dw$ID, 
                          type=paste(sampleB,"_lowerThan_", 
                                     sampleA, sep="."))
de_id = rbind(res_de_up_id, res_de_dw_id)



# 火山图
cut_off_pvalue = 0.01  #统计显著性
cut_off_logFC = 1
res$change = ifelse(res$pvalue < cut_off_pvalue & abs(res$log2FoldChange) >= cut_off_logFC,
                        ifelse(res$log2FoldChange> cut_off_logFC ,'Up','Down'),
                        'Nodiff')
head(res)
this_tile <- paste0(
  '\nThe number of up gene is ',nrow(res[res$change =='Up',]),
  '\nThe number of down gene is ',nrow(res[res$change =='Down',])
)

p <- ggplot(
  # 数据、映射、颜色
  dataa, aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
  geom_point(alpha=1, size=dataa$Size) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) +
  ggtitle( this_tile )
p



# heatmap
# 按pvalue排序, 把差异大的基因放前面
res <- res[order(res$pvalue),]
head(res)


# 差异基因筛选，padj<0.1
data_de_up_top_id <- as.vector(head(data_de_up$ID,10))
data_de_dw_top_id <- as.vector(head(data_de_dw$ID,10))


red_de_top <- c(data_de_up_top_id, data_de_dw_top_id)
# red_de_top

red_de_top_expr <- normalized_counts[normalized_counts$ID %in% red_de_top,]
rownames(red_de_top_expr) <- red_de_top_expr$ID
red_de_top_expr$ID <- NULL
# head(red_de_top_expr)


pheatmap(red_de_top_expr,show_colnames =T,show_rownames = T,
         scale="row",
         color = colorRampPalette(c("blue", "white", "red"))(10),
         cluster_cols = T)






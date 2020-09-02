### NOTE: RNA-seq expression file should be exported as "tsv" file
### NOTE2: Clients should prepare a sample-list file indicate the experiment design info
### NOTE2: Sample-list file should be a "csv" file ###
### NOTE2: Sample-list file should have the following feature columns: ###
### NOTE2: samples(sample names), conditions(usually the experiment design grouping info)
### NOTE3: RNA-seq expression data should have columns names which are same as the sample-list file "samples" colum ###

###!!!Client options: input RNA-seq expression file
datafile <- "./testData/RNAseq_test/input/GSE148729_Caco2_polyA_readcounts.tsv"
###!!!Client options: input sample-list file
sampleList <- "./testData/RNAseq_test/input/fileSample.csv"
###!!!Client options: control group, default will use group of the first column sample / first sample in sample-list file
controlGrp <- "untr"
#!!!Client options: Client should choose which group to do the comparison(DEA) or do group by group comparison("group_by_group")
experGrp <- "S2"

##!!!!!DEV: set the directory that plot files locate
fileLoc <- "./testData/RNAseq_test/output/"

library(DESeq2)
library(apeglm)
library(ggplot2)
library(dplyr)
library(pheatmap)
options(stringsAsFactors = F)
data <- read.table(datafile,
                   header=T, row.names=1, com='', 
                   quote='',check.names=F, sep="")
sampleInfo <- read.csv(sampleList)
data <- select(data, sampleInfo$samples)
groupsLevel <- unique(sampleInfo$conditions)
sampleInfo$conditions <- factor(sampleInfo$conditions, 
                                levels = c(groupsLevel[groupsLevel == controlGrp], groupsLevel[groupsLevel != controlGrp]))

## Delete low gene abundance feature
data <- data[rowSums(data)>2, ]

## DESeq2 procedure
ddsFullCountTable <- tryCatch(DESeqDataSetFromMatrix(countData = data,
                                            colData = sampleInfo, 
                                            design= ~ conditions), 
              error = function(e){if(conditionMessage(e) == "some values in assay are not integers"){
                cat("Convert some float to integer!\n")
                data <- round(data)
                DESeqDataSetFromMatrix(countData = data,
                                       colData = sampleInfo, 
                                       design= ~ conditions)
              }else{
                stop(paste0("Error in read counts input: ", conditionMessage(e), 
                            ". PROGRAM EXIT!"))
              }})
dds <- DESeq(ddsFullCountTable)
#Normalization 
normalized_counts <- counts(dds, normalized=TRUE)
#Sorting by mad(Median Absolute Deviation) value(Genes with greater differences are ranked higher)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts, paste0(fileLoc, "normalized_counts_sorted.csv"))
#log2 transformation
rld <- rlog(dds, blind=FALSE)
#PCA plot
pca_data <- plotPCA(rld, intgroup = "conditions",returnData = T)
percentVar <- round(100*attr(pca_data, "percentVar"), 1) 
pcaScorePlot <- ggplot() +
  geom_point(pca_data, mapping =aes(x = PC1, y = PC2, color = conditions)) +
  scale_color_aaas() +
  scale_fill_aaas() + 
  theme_bw() + 
  labs(x = paste0("PC1", "(", percentVar[1], "%)"), 
       y = paste0("PC2", "(", percentVar[2], "%)"), 
       color = "group", 
       title = "PCA Score Plot") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
ggsave(paste0(fileLoc, "PCA_score_plot_all.pdf"), plot = pcaScorePlot, 
       device = "pdf", width = 9, height = 9)

## Differential expression analysis
#!!!Client options: Fold change threshold
fcthresh <- 2.0
#!!!Client options: Fold change p.value
pthresh <- 0.1
#!!!Client options: Show top gene labels in volcano or set "F" to do not show any
showtop <- 20
#!!!Client options: Show top gene in heatmap
showtop2 <- 75
res <- results(dds, contrast = c("conditions", experGrp, controlGrp))
if(experGrp != "group_by_group"){
  coefVar <- paste0("conditions_", experGrp, "_vs_", controlGrp)
  resLFC <- tryCatch(lfcShrink(dds, coef=coefVar, type="apeglm"), 
                     error = function(e){
                       cat("Using DESeq2 default method to calculate LFC\n")
                       return(lfcShrink(dds, coef=coefVar, type="normal"))
                     })
  ## Volcano
  volcano.data <- do.call(cbind, resLFC@listData) %>%
    as.data.frame() %>%
    mutate(
      log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
      padj = ifelse(is.na(padj), 1, padj),
      gene = resLFC@rownames, 
      up = ifelse(log2FoldChange > log(fcthresh, 2) & padj < pthresh, 1, 0),
      down = ifelse(log2FoldChange < -log(fcthresh, 2) & padj < pthresh, 2, 0),
      regState = sapply(up + down, function(x) switch(x, "upreg", "downreg")), 
      regState = sapply(regState, function(x) ifelse(is.null(x), 0, x))
    )
  volcano.data_reg <- subset(volcano.data, subset = regState != 0)
  volcano.data_unreg <- subset(volcano.data, subset = regState == 0)
  #!!!!!WARNING: use p_value as rank standard, may use others(eg. mad?)
  if(showtop != F){
    if(showtop > nrow(volcano.data_reg)){
      cat("Not enough significant feature! Show all the feature label.\n")
    }
    volcano.data_top <- volcano.data_reg[order(volcano.data_reg$padj), ][1:showtop, ]
  }else{
    volcano.data_top <- volcano.data_reg[order(volcano.data_reg$padj), ][0, ]
  }
  volcano.plot <- ggplot(mapping = aes(x = log2FoldChange, y = -log(padj, 10))) +
    geom_point(data = volcano.data_reg, mapping = aes(color = factor(regState))) + 
    geom_point(data = volcano.data_unreg, color = "gray") +
    geom_text(data = volcano.data_top, mapping = aes(label = gene), 
              hjust = 0, nudge_x = 0.05) +
    theme_bw() +
    geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
    labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "regulation", 
         title = paste0(experGrp, " vs ", controlGrp)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))+ 
    scale_color_aaas()
  volcano.plot
  ## DE gene output
  getMeanCounts <- function(counts = normalized_counts, group){
    data <- counts[, colData(dds)$conditions == group]
    countmean <- apply(data, 1, mean)
    countmean <- countmean[match(volcano.data_reg$gene, names(countmean))]
    return(countmean)
  }
  #use reg feature that volcano used
  sigDEgene <- volcano.data_reg %>%
    mutate(controlGrp = getMeanCounts(group = controlGrp), 
           experGrp = getMeanCounts(group = experGrp)) %>%
    select(c('gene', 'controlGrp', 'experGrp', 'log2FoldChange', 'padj')) %>%
    rename(controlGrp = controlGrp, experGrp = experGrp)
  write.csv(sigDEgene, paste0(fileLoc, "sigDEgeneStatistics.csv"), 
            row.names = F)
  ## Heatmap
  if(showtop2 > nrow(volcano.data)){
    cat("Not enough significant feature! Show all the feature label.\n")
  }
  #!!!!!WARNING: use p_value as rank standard, may use others(eg. mad?)
  heatmap.data_topgene <- volcano.data[order(volcano.data$padj), ][1:showtop2, ]$gene
  heatmap.data_top <- normalized_counts[match(heatmap.data_topgene, rownames(normalized_counts)), 
                                        colData(dds)$conditions %in% c(controlGrp, experGrp)]
  legend_annotation <- data.frame(conditions = sampleInfo$conditions)
  rownames(legend_annotation) <- sampleInfo$samples
  heatmap <- pheatmap::pheatmap(mat = heatmap.data_top, 
                     annotation = legend_annotation, 
                     fontsize = 8, 
                     fontsize_row = 8, 
                     cluster_rows = T, 
                     cluster_cols = F, 
                     scale = "row"
  )
  pdf(paste0(fileLoc, "heatmap_top", showtop2, ".pdf"), 
      width=(ncol(heatmap.data_top)*25+300)/72, height=(nrow(heatmap.data_top)*18+150)/72)
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  dev.off()
}
#!!!!!Control flow WARNING
if(experGrp == "group_by_group"){
  
}


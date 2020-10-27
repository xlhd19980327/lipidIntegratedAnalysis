### NOTE: RNA-seq expression file should be exported as "tsv"/"csv" file(Integrate "csv" inputs into final function already)
### NOTE2: Clients should prepare a sample-list file indicate the experiment design info
### NOTE2: Sample-list file should be a "csv" file ###
### NOTE2: Sample-list file should have the following feature columns: ###
### NOTE2: samples(sample names), conditions(usually the experiment design grouping info)
### NOTE2: Sample-list file may(optical) have the following feature columns: ###
### NOTE2: batch(batch info, for removing batch effect) ###
### NOTE2: if batch colum exisits, program will automatically remove batch effect, ###
### Developer can add a option about removing batch effect(Not achieve yet) 
### NOTE3: RNA-seq expression data should have gene info(gene_id/gene name) in the 1st column ###
### NOTE3: RNA-seq expression data should have columns names which are same as the sample-list file "samples" colum ###
### NOTE3: DESeq2 notes: RNA-seq expression data should use un-normalized counts(Raw counts) as input ###
### NOTE3: (ref: http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

###!!!Client options: input RNA-seq expression file
datafile <- "./testData/RNAseq_test/input/GSE148729_Caco2_polyA_readcounts.tsv"
###!!!Client options: input sample-list file
sampleList <- "./testData/RNAseq_test/input/fileSample_batch.csv"
###!!!Client options: control group, default will use group of the first sample in sample-list file
controlGrp <- "untr"
#!!!Client options: Client should choose which group to do the comparison(DEA) or do group by group comparison("group_by_group")
experGrp <- "S2"

##!!!!!DEV: set the directory that plot files locate
fileLoc <- "./testData/RNAseq_test/output/"

source("./utilityFunc/plottingPalettes.R")

library(DESeq2)
library(limma)
library(apeglm)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggrepel)
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
batremvOpt <- T
ddsFullCountTable <- tryCatch(DESeqDataSetFromMatrix(countData = data,
                                            colData = sampleInfo, 
                                            design= parse(text = ifelse(batremvOpt, "~ conditions + batch", "~ conditions"))), 
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
#Normalization(note: only providing counts scaled by size or normalization factors, 
#experiment design info only used in the other normalization methods)
normalized_counts <- counts(dds, normalized=TRUE)
#Sorting by mad(Median Absolute Deviation) value(Genes with greater 'overall' differences are ranked higher)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts, paste0(fileLoc, "normalized_counts_sorted.csv"))
#log2 transformation and normalized(may use "vst" instead)
rld <- rlog(dds, blind=FALSE)
if(batremvOpt){
  ## Using limma::removeBatchEffect to remove batch effect, it will change the original count matrix
  assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)
}
#PCA plot
pca_data <- plotPCA(rld, intgroup = ifelse(batremvOpt, c("conditions", "batch"), "conditions"),returnData = T)
ellipse_coor_all <- data.frame()
for(i in groupsLevel){
  groupVar <- var(subset(pca_data, subset = group == i)[, 1:2], na.rm = T)
  groupMean <- cbind(mean(subset(pca_data, subset = group == i)[, 1], na.rm = T), 
                     mean(subset(pca_data, subset = group == i)[, 2], na.rm = T))
  ellipse_coor <- as.data.frame(ellipse::ellipse(groupVar, centre = groupMean, 
                                                 level = 0.95, npoints = 100))
  colnames(ellipse_coor) <- c("x", "y")
  ellipse_coor_all <- rbind(ellipse_coor_all, 
                            cbind(ellipse_coor, group = i))
}
percentVar <- round(100*attr(pca_data, "percentVar"), 1) 
colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
colorpars_fill <- alpha(colorpars, 0.3)
names(colorpars_fill) <- names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
pcaScorePlot <- ggplot() +
  geom_polygon(ellipse_coor_all, mapping = aes(x = x, y = y, 
                                               color = factor(group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])), 
                                               fill = factor(group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))) +
  geom_point(pca_data, mapping =aes(x = PC1, y = PC2, color = group)) +
  scale_color_manual(values = colorpars) +
  scale_fill_manual(values = colorpars_fill) +
  theme_bw() + 
  labs(x = paste0("PC1", "(", round(percentVar[1], 1), "%)"), 
       y = paste0("PC2", "(", round(percentVar[2], 1), "%)"), 
       color = "group", 
       fill = "group",
       title = "PCA Score Plot") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
ggsave(paste0(fileLoc, "PCA_score_plot_all.pdf"), plot = pcaScorePlot, 
       device = "pdf", width = 9, height = 9)

## Differential expression analysis
## Usually will be two groups comparison
## This section contains DEA & Volcano-plot & Heatmap, using the following procedure(same as the DESeq DEA procedure):
## DEA: dds <- DESeq(ddsFullCountTable); res <- results(dds); res
## Volcano: lfcShrink(dds, ...)
## Heatmap: use volcano.data to rank genes, plot top genes with normalized_counts
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
  ## DEA output
  deaResult <- as.data.frame(res) %>%
    rownames_to_column(var = "Gene")
  write.csv(deaResult, paste0(fileLoc, "DEgeneStatistics.csv"), 
            row.names = F)

  
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
    geom_text_repel(data = volcano.data_top, mapping = aes(label = gene), 
              hjust = 0, nudge_x = 0.05) +
    theme_bw() +
    geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
    labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "regulation", 
         title = paste0(experGrp, " vs ", controlGrp)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))+ 
    scale_color_aaas()
  ggsave(paste0(fileLoc, "volcano_", experGrp, "_vs_", controlGrp, ".pdf"), plot = volcano.plot, 
         device = "pdf", width = 9, height = 9)
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
                     color = plottingPalettes(100, type = "continuous"),
                     #one-by-one usage
                     annotation_colors = list(conditions = colorpars[names(colorpars) %in% c(experGrp, controlGrp)]),
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


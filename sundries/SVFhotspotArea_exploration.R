# Get hotspot area and differentiation-related genes
library(clusterProfiler)
library(cowplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
#!!!Client options: species("hsa"/"mmu")
#spe <- opt$species
#!!!Client options: gene type(sugg: "ENSEMBL"/"SYMBOL")(ref: keytypes([spe db]))
#gene_type <- opt$gene_type
#!!!Client options: show ontology item number
#shownum <- opt$show_num
#!!!Client options: select go term catalog("BP"/"CC"/"MF"/"ALL")
#gocat <- opt$go_term
orgdb <- org.Mm.eg.db
goopt <- "BP"
genes <- read.csv("~/temp/genes_3.csv")$x
go <- enrichGO(gene = genes, OrgDb = orgdb, ont=goopt,pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = "SYMBOL")
barplot(go,showCategory=50,drop=T, font.size=15)
genes <- go@result[go@result$ID %in% "GO:0045444", ]$geneID
genes <- strsplit(genes, split = "/")[[1]]

# Get genes normalized counts
a <- read.csv("./testData/SVFmultiomics_210118/output/normalized_counts_sorted.csv")
b <- a[a$X %in% genes, ]
rownames(b) <- b$X
b <- b[, -1]

# Do heatmap
source("./utilityFunc/readingRNAData.R")
source("./utilityFunc/DESeq2preproc.R")
source("./utilityFunc/rnaPCAPlot.R")
source("./utilityFunc/DEAnalysis.R")
source("./utilityFunc/rnaVolcanoPlot.R")
source("./utilityFunc/rnaHeatmapPlot.R")
source("./utilityFunc/plottingPalettes.R")
library(multcomp)
#Color
library(RColorBrewer)
library(ggsci)
#RNA-seq
library(DESeq2)
library(limma)
library(apeglm)
library(tidyverse)
library(pheatmap)
library(ggrepel)
options(stringsAsFactors = F)
dataSet_RNA <- readingRNAData(datafile = "./testData/SVFmultiomics_210118/input/RNAseq_genesymbol.csv", 
                              sampleList = "./testData/SVFmultiomics_210118/input/sampleList.csv", 
                              controlGrp = "Day0")
dataProc_RNA <- DESeq2preproc(dataSet = dataSet_RNA, 
                              fileLoc = "~/temp/")
#heatmap preparation
DEAresult = dataProc_RNA
controlGrp <- DEAresult$dataSet$controlGrp
sampleInfo <- DEAresult$dataSet$sampleInfo
groupsLevel <- DEAresult$dataSet$groupsLevel
colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
legend_annotation <- data.frame(conditions = sampleInfo$conditions)
rownames(legend_annotation) <- sampleInfo$samples
collegend_annotation <- list(conditions = colorpars)
#heatmap
heatmap <- pheatmap::pheatmap(mat = b, 
                              annotation = legend_annotation, 
                              fontsize_col = 20, 
                              fontsize_row = 20, 
                              fontsize = 14,
                              annotation_colors = collegend_annotation,
                              cluster_rows = T, 
                              cluster_cols = F, 
                              scale = "row"
)
pdf(paste0("~/temp/heatmap_difReGenes.pdf"), 
    width=9, 
    height=10)
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
dev.off()
write.csv(b, "~/temp/normalized_counts_difReGenes.csv")

# Get lipid-genes corelation data(in that hotspot area), do heatmap 
cordata <- read.csv("./sundries/bigdata/Fig4_cor_data/correlation_7_3.csv", 
                    row.names = 1)
sub_cordata <- cordata[, colnames(cordata) %in% genes]
tmp <- sub_cordata[order(rownames(sub_cordata)), ]
write.csv(tmp, "~/temp/data_cor_gene-lip.csv")
x <- pheatmap::pheatmap(mat = tmp, 
                        cluster_rows = F, 
                        cluster_cols = F, 
                        scale = "row"
)
pdf(paste0("~/temp/heatmap_cor_gene-lip.pdf"), 
    width=10, height=25)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()

a <- read.csv("./branch/benchmark/input/HANlipid_orig.csv")
lipid <- paste0(a$class, "(", a$Compound, ")")
b <- subset(a, select = c(-class, -Compound))
d <- cbind(lipid, b)
write.csv(d, "./branch/benchmark/input/HANlipid_onlyforallclass.csv", row.names = F)



source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/MARpreproc.R")
source("./utilityFunc/headgroupStat.R")
outputLoc <- "./branch/benchmark/output/"
dataSet <- readingLipidData(datafile = "./branch/benchmark/input/HANlipid_onlyforallclass.csv",
                            sampleList = "./branch/benchmark/input/HANsampleList_lipid.CSV", 
                            controlGrp = "Day0", dataType = "LipidSearch", delOddChainOpt = F,
                            lipField = "lipid", #fileLoc = outputLoc, 
                            na.char = "")
mSet <- MARpreproc(dataSet = dataSet, fileLoc = outputLoc)
headgroupStat(dataSet = dataSet, mSet = mSet, ignore = T, 
              fileLoc = paste0(outputLoc, "headgroup/"))

## RNA heatmap, top, with multiple groups
DEAresult = dataProc_RNA
controlGrp <- DEAresult$dataSet$controlGrp
sampleInfo <- DEAresult$dataSet$sampleInfo
groupsLevel <- DEAresult$dataSet$groupsLevel

colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
legend_annotation <- data.frame(conditions = sampleInfo$conditions)
rownames(legend_annotation) <- sampleInfo$samples

if(type == "RNAseq"){
  normalized_counts <- DEAresult$normalized_counts
}
if(type == "MiAr"){
  normalized_counts <- DEAresult$eset
}

normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ][1:75, ]

##We really donot recommand client to use showtop
#heatmap.data_top <- normalized_counts[1:showtop, ]
heatmap.data_top <- normalized_counts
collegend_annotation <- list(conditions = colorpars)
heatmap <- pheatmap::pheatmap(mat = heatmap.data_top, 
                              annotation = legend_annotation, 
                              #ugly
                              #color = plottingPalettes(100, type = "continuous"),
                              #one-by-one usage
                              annotation_colors = collegend_annotation,
                              fontsize = 20,
                              cluster_rows = T, 
                              cluster_cols = T, 
                              scale = "row"
)
pdf(paste0("~/temp/heatmap_allgroups.pdf"), width=10, height=25)
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
dev.off()

### SVF multi-omics correlation/enrich
kg <- 3
kl <- 5
png("~/temp/cor/correlationPlot.png", 
    width = 720, height = 720)
list <- pheatmap(correlation,
                 cutree_col = kg, cutree_rows = kl)
dev.off()

a <- read.csv("/home/lifs/temp/cor3/genes_3.csv")
library(clusterProfiler)
library(cowplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
go_BP <- enrichGO(gene = a$x, OrgDb = 'org.Mm.eg.db', ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,keyType = "SYMBOL")
p <- barplot(go_BP,showCategory=10,drop=T, font.size = 15) +
  ggtitle("Biological Process") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
ggsave("~/temp/GOenrich_Biological_Process.pdf", plot = p, 
       device = "pdf", width = 8, height = 15/50*10, limitsize = FALSE)

### SVF multi-omics cumulative percent contain
a <- data_sub_classSum_stat2 %>%
  group_by(group) %>%
  mutate(perc = realmean / sum(realmean))
plot_cum_a <- ggplot(data = a, aes(x = group, y = perc, fill = Class)) + 
  geom_bar(stat="identity", width = 0.7)+
  theme_classic()+
  labs(title = "Cumulative lipid composition", 
       x = "group",
       y = "total concentration") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25), 
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
  scale_fill_manual(values = plottingPalettes(nClass, type = "discrete"))
ggsave(paste0(fileLoc, "headgroup_cumperc_", pname, ".pdf"), 
       plot = plot_cum_a, device = "pdf", 
       width = 6/7.7*(3.3/3*length(groupsLevel)+4.4), height = 10)



### SVF multi-omics differentiation-related genes heatmap
# find normalized data
a <- read.csv("./testData/SVFmultiomics_210118/output/normalized_counts_sorted.csv")
mygenes <- read.csv("./testData/SVFmultiomics_210118/input/RNAseq_difRel.csv")
b <- a[a$X %in% mygenes$gene_id, ]
rownames(b) <- b$X
b <- select(b, select = -X)
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
pdf(paste0("./testData/SVFmultiomics_210118/output/", "heatmap_allgroups.pdf"), 
    width=9, 
    height=8)
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
dev.off()
write.csv(b, "~/temp/normalized_counts_lipgenes.csv")



### SVF multi-omics find differentiation-related genes
genes <- dir("~/temp/cor3/", pattern = "^genes_")
alldata <- data.frame(
  x = character(),
  i = character()
)
for(i in 1:6){
  j <- paste0("genes_", i, ".csv")
  data <- read.csv(paste0("~/temp/cor3/", j))
  data <- data.frame(data, i)
  alldata <- rbind(alldata, data)
}
mygenes <- read.csv("./testData/SVFmultiomics_210118/input/RNAseq_difRel.csv")
mydata <- subset(alldata, subset = x %in% mygenes$gene_id)
write.csv(mydata, "~/temp/gene_loc.csv")

### SVF multi-omics differentiation-related genes heatmap
#b is above(subset of the large one) or use tidy format as seen below
b <- read.csv("~/temp/normalized_counts_lipgenes2.csv", row.names = 1)
correlation <- cor(lipid_data, t(b), method ="spearman")
correlation[is.na(correlation)] <- 0
list <- pheatmap::pheatmap(correlation, 
                 cutree_col = 4, cutree_rows = 12)
pdf("~/temp/cor/correlationPlot.pdf", 
    width = 5)
list <- pheatmap::pheatmap(correlation)
dev.off()

### See the lipid Class
lipids <- dir("~/temp/cor3/", pattern = "lipids_")
i <- lipids[1]
for(i in lipids){
  j <- read.csv(paste0("~/temp/cor3/", i))
  ind <- regexpr("[0-9]+", i)
  num <- substr(i, ind, attr(ind, "match.length")+ind-1)
  Class <- unique(gsub("(.*?)\\(.*", "\\1", j$x))
  myClass <- paste0(Class, collapse = ", ")
  write.csv(myClass, paste0("~/temp/cor3/lipidClass_", num, ".csv"))
}

### See the original data, in lipid class aspect
lipids <- dir("~/temp/cor3/", pattern = "lipids_")
i <- lipids[5]

j <- read.csv(paste0("~/temp/cor3/", i), stringsAsFactors = F)
ind <- regexpr("[0-9]+", i)
num <- substr(i, ind, attr(ind, "match.length")+ind-1)
subdata <- t(mSet$dataSet$proc[j$x])
dataSet$lipidName <- rownames(subdata)
rownames(subdata) <- NULL
subdata <- as.data.frame(subdata)
dataSet$data <- subdata
headgroupStat(dataSet = dataSet, mSet = NULL, 
              fileLoc = "~/temp/cor3/headgroup/")
lipSumClassHeatmapPlot(dataSet = dataSet, mSet = NULL, 
                       fileLoc = "~/temp/cor3/headgroup/")


### See the top correlation lipids(top9)
options(stringsAsFactors = F)
a <- read.csv("~/temp/cor3/correlation_7_3.csv", row.names = 1)
b <- read.csv("~/temp/cor3/correlation_2_4.csv", row.names = 1)
ab <- rbind(a, b)
corAvg <- apply(ab, 1, mean)
corAvg <- apply(a, 1, mean)
top9 <- names(corAvg[order(corAvg, decreasing = T)][1:9])
subdata <- as.data.frame(t(mSet$dataSet$proc[top9]))
subdata_plot <- subdata %>%
  rownames_to_column(var = "lipidName") %>%
  gather(key = "case", value = "conc", -lipidName) %>%
  mutate(group = dataSet$allgroups[case])
subdata_plot$group <- factor(subdata_plot$group, 
                             levels = c(dataSet$controlGrp, 
                                        unique(dataSet$allgroups[dataSet$allgroups != dataSet$controlGrp])))
plot_color <- ggplot(data = subdata_plot, aes(x = group, y = conc, color = group)) + 
  geom_boxplot(width = 0.8, size = 1.2)+ 
  facet_wrap(~lipidName, scales="free", ncol = 3) +
  scale_color_npg()+
  labs(x = "group",
       y = "Concentration",
       color = "group",
       title = "Lipids statistics"
       ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
    axis.title = element_text(size = 30), 
    axis.text = element_text(size = 30),
    legend.text = element_text(size = 30), 
    legend.title = element_text(size = 30), 
    strip.text = element_text(hjust = 0.5, size = 30, face = "bold")#,
    #axis.text.x = element_text(angle = 45)
  ) 
nclass <- 9
##used for color plot
w_plot <- 25/956*(82+(259/5*5+(288-259))*ifelse(nclass >= 3, 3, nclass)) +3
##used for color plot
h_plot <- 33/1264*(1224/4*ifelse(nclass%%3, nclass%/%3+1, nclass%/%3)+40) -3
ggsave("~/temp/top9.pdf", 
       plot = plot_color, device = "pdf", width = w_plot, height = h_plot, limitsize = FALSE)

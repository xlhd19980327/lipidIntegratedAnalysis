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

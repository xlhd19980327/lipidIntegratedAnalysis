outputLoc <- "~/temp/"
prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}
analOpt <- "E6"
dataSet <- readingLipidData(datafile = "/home/rivergroup/input/20211211_LP_mcx_Lipid_Data_forLINT_norm.csv",
                            sampleList = "/home/rivergroup/input/desc_forLINT.csv", 
                            controlGrp = "C6", dataType = "Lipids", delOddChainOpt = T)

dataset <- prepDataSet(analOpt)
mSet <- MARpreproc(dataSet = dataset, fileLoc = outputLoc, perc = 2/3*100, normOpt = "C")
data_type <- dataset$dataType
data_proc <- t(mSet[["dataSet"]][["proc"]])
dataset$lipidName <- rownames(data_proc)
rownames(data_proc) <- NULL
data_proc <- as.data.frame(data_proc)
dataset$data <- data_proc



### Heatmap with lipid class label(top xxx) ###
dataSet <- dataset
groupsLevel <- dataSet$groupsLevel
controlGrp <- dataSet$controlGrp
data_norm <- mSet$dataSet$norm
mSetObj <- mSet
#data_prenorm <- mSet$dataSet$proc
# Statistics 
topnum_in <- 75
if(length(groupsLevel) == 2){
  mSetObj <- Ttests.Anal(mSetObj)
  toplip <- names(sort(mSetObj$analSet$tt$p.value))[1:topnum_in]
  toplip_ind <- match(toplip, colnames(data_norm))
}else if(length(groupsLevel) > 2){
  mSetObj <- ANOVA.Anal(mSetObj)
  toplip <- names(sort(mSetObj$analSet$aov$p.value))[1:topnum_in]
  toplip_ind <- match(toplip, colnames(data_norm))
}
heatmapdata <- data_norm[toplip_ind]
tidy_ind <- match(names(sort(factor(dataSet[["allgroups"]], levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))), 
                  row.names(heatmapdata))
heatmapdata <- heatmapdata[tidy_ind, ]
datagroup <- factor(dataSet[["allgroups"]][match(row.names(heatmapdata), names(dataSet[["allgroups"]]))], 
                    levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
data <- t(heatmapdata)[apply(t(heatmapdata), 1, function(x) ifelse(is.na(sd(x, na.rm = T)), F, ifelse(sd(x, na.rm = T) == 0, F, T))), ]
# Get lipid class
source("./utilityFunc/getClassInfo.R")
classinfo <- sapply(rownames(data), getClassInfo, "LipidSearch")
annotation_row = data.frame(
  classinfo = factor(classinfo)
)
rownames(annotation_row) = rownames(data)
colorpars1 <- plottingPalettes(n = length(groupsLevel), type = "discrete")
colorpars2 <- plottingPalettes(n = length(levels(factor(classinfo))), type = "discrete")
data2 <- data[order(rownames(data)), ]
x <- pheatmap::pheatmap(mat = data2,
                       annotation = data.frame(group = factor(datagroup), row.names = rownames(heatmapdata)), 
                       annotation_row = annotation_row, 
                       fontsize_row = 5,
                       cluster_rows = F, 
                       cluster_cols = F,
                       scale = "row", 
                       #annotation_colors = list(group = colorpars1, classinfo = colorpars2),
                       clustering_distance_rows = "euclidean"
)
pdf(paste0(outputLoc, "heatmap_top", topnum_in, ".pdf"), width = 4)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()


### Heatmap about some specific lipids ###
data <- read.csv("~/temp/PC6.csv", row.names = 1)
subdata <- data[grepl("18:0", rownames(data)), ]
p.value <- apply(subdata, 1, function(x) t.test(x[1:3], x[4:6])$p.value)
p.value[p.value < 0.05]
x <- pheatmap::pheatmap(mat = subdata,
                        #fontsize_row = 5,
                        cluster_rows = F, 
                        cluster_cols = F,
                        scale = "row"
)
pdf(paste0(outputLoc, "heatmap_cus.pdf"), width = 3, height = 1.5)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()
#-----#
data <- read.csv("~/temp/DG6.csv", row.names = 1)
label <- c("16:0", "16:1", "18:1", "18:2", "18:3")
subdata <- data[grepl(paste(label, collapse = "|"), rownames(data)), ]
ind <- sapply(label, grep, rownames(subdata))
indv <- do.call(c, ind)
tidydata <- subdata[indv, ]
annotation_row = data.frame(
  info = rep(label, sapply(ind, length))
)
rownames(annotation_row) = rownames(tidydata)
x <- pheatmap::pheatmap(mat = tidydata,
                        annotation_row = annotation_row,
                        #fontsize_row = 5,
                        cluster_rows = F, 
                        cluster_cols = F,
                        scale = "row"
)
pdf(paste0(outputLoc, "heatmap_cus.pdf"), width = 4, height = 7)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()

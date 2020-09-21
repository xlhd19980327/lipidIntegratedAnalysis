### Get MS2 heatmap ###

## Source will offer the following contents:
## Function(s): getFAsInfo
if(dataType == "LipidSearch"){
  source("./utilityFunc/getFAsInfo.R")
}
if(dataType == "MS_DIAL"){
  source("./utilityFunc/getFAsInfo_msdial.R")
}

## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
lipids <- dataset$lipidName
# Get the FA info and MS1 info
fasInfo <- lapply(lipids, getFAsInfo)
ms1Info <- sapply(fasInfo, function(x) x$ms1)
data_tidy <- cbind(lipids, dataset$data, ms1 = ms1Info) %>%
  filter(ms1 == F) %>%
  select(-ms1)
dataset_tidy <- list(
  data = subset(data_tidy, select = -lipids), lipidName = data_tidy$lipids, 
  groupsLevel = dataset$groupsLevel, allgroups = dataset$allgroups, controlGrp = dataset$controlGrp, nsamples = dataset$nsamples, 
  dataType = dataset$dataType
)
mSet <- MARpreproc(dataSet = dataset_tidy)

dataSet = dataset_tidy
mSet = mSet
topnum = 75

groupsLevel <- dataSet$groupsLevel
controlGrp <- dataSet$controlGrp
data_norm <- mSet$dataSet$norm
data_prenorm <- mSet$dataSet$prenorm
colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])

## Do heatmap, top 75, norm == T
mSetObj = mSet
topnum_in = topnum
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
data <- t(heatmapdata)[apply(t(heatmapdata), 1, function(x) sd(x)!=0), ]
row.names(data) <- gsub(".*\\|(.*)", "\\1", row.names(data))
x <-pheatmap::pheatmap(mat = data, 
                       ##This color palette is ugly, use default
                       #color = plottingPalettes(100, type = "continuous"),
                       annotation = data.frame(group = datagroup, row.names = rownames(heatmapdata)), 
                       fontsize = 8, 
                       fontsize_row = 8, 
                       clustering_distance_rows = "euclidean", 
                       clustering_distance_cols = "euclidean", 
                       clustering_methods = "ward.D", 
                       cluster_rows = T, 
                       cluster_cols = F, 
                       scale = "row", 
                       annotation_colors = list(group = colorpars)
)








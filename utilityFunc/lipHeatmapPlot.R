lipHeatmapPlot <- function(dataSet, mSet, 
                           fileLoc, 
                           topnum = NA){
  groupsLevel <- dataSet$groupsLevel
  controlGrp <- dataSet$controlGrp
  data_norm <- mSet$dataSet$norm
  data_prenorm <- mSet$dataSet$prenorm
  colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
  names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
  
  
  ## Heatmap
  plotHeatmap <- function(#mSetObj = mSet, 
                          norm = T, topind = NA){
    if(norm == T){
      heatmapdata <- data_norm
    }else{
      heatmapdata <- data_prenorm
    }
    if(!all(is.na(topind))){
      heatmapdata <- heatmapdata[topind]
    }
    tidy_ind <- match(names(sort(factor(dataSet[["allgroups"]], levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))), 
                      row.names(heatmapdata))
    heatmapdata <- heatmapdata[tidy_ind, ]
    datagroup <- factor(dataSet[["allgroups"]][match(row.names(heatmapdata), names(dataSet[["allgroups"]]))], 
                        levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
    data <- t(heatmapdata)[apply(t(heatmapdata), 1, function(x) sd(x)!=0), ]
    if(nrow(data) < ncol(heatmapdata)){
      cat("Some lipids in low variation will be removed.\n")
    }
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
    pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                           "_vs_", controlGrp))
    if(all(is.na(topind))){
      pdf(paste0(fileLoc, "heatmap_", pname, ".pdf"), width=(nrow(heatmapdata)*25+300)/72, height=(ncol(heatmapdata)*18+150)/72)
    }else{
      pdf(paste0(fileLoc, "heatmap_top", length(topind), "_", pname, ".pdf"), width=(nrow(heatmapdata)*25+300)/72, height=(ncol(heatmapdata)*18+150)/72)
    }
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  plotSubHeatmap <- function(mSetObj = mSet, 
                             norm = T, topnum_in = topnum){
    if(length(groupsLevel) == 2){
      mSetObj <- Ttests.Anal(mSetObj)
      toplip <- names(sort(mSetObj$analSet$tt$p.value))[1:topnum_in]
      toplip_ind <- match(toplip, colnames(data_norm))
    }else if(length(groupsLevel) > 2){
      mSetObj <- ANOVA.Anal(mSetObj)
      toplip <- names(sort(mSetObj$analSet$aov$p.value))[1:topnum_in]
      toplip_ind <- match(toplip, colnames(data_norm))
    }
    plotHeatmap(norm = norm, topind = toplip_ind)
  }
  if(is.na(topnum)){
    #!!!Client options
    plotHeatmap(norm = T, topind = NA)
  }else{
    #!!!Client options
    plotHeatmap(norm = T, topind = NA)
    #!!!Client options
    plotSubHeatmap(mSet, 
                   norm = T)
  }
}

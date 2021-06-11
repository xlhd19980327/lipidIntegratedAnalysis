lipSumClassHeatmapPlot <- function(dataSet, mSet, 
                                   ignore = T, show_detail = F,
                                   fileLoc){
  dataType <- dataSet$dataType
  allgroups <- dataSet$allgroups
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  
  source("./utilityFunc/plottingPalettes.R")
  ## Source will offer the following contents:
  ## Function(s): getClassInfo
  source("./utilityFunc/getClassInfo.R")
  
  data_tidy <- dataSet[["data"]] %>%
    mutate(lipidName = dataSet[["lipidName"]], 
           Class = switch(dataSet$dataType,
                          LipidSearch = sapply(lipidName, getClassInfo, "LipidSearch", ignore = ignore), 
                          MS_DIAL = sapply(lipidName, getClassInfo, "MS_DIAL", ignore = ignore)))
  nclass <- length(unique(data_tidy$Class))
  
  datagroup <- factor(allgroups, 
                      levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
  colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
  names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
  
  if(show_detail == F){
    data_sub_classSum_stat <- data_tidy %>%
      ungroup() %>%
      dplyr::select(-lipidName) %>%
      gather(key = "case", value = "lipidsum", -Class) %>%
      group_by(Class, case) %>%
      summarise(lipidsum = sum(lipidsum, na.rm = T)) %>%
      mutate(group = allgroups[match(case, names(allgroups))]) 
    data_Class <- data_sub_classSum_stat %>% 
      select(-group) %>% 
      spread(key = "case", value = "lipidsum") %>%
      ungroup(Class) %>%
      column_to_rownames(var = "Class")
    data_Class <- data_Class[, match(names(allgroups), colnames(data_Class))]
    data_Class <- data_Class[apply(data_Class, 1, function(x) sd(x)!=0), ]
    # if(nrow(data_Class) < nrow(data_Class)){
    #   cat("Some lipids in low variation will not show in the plot.\n")
    # }
    
    x <- pheatmap::pheatmap(mat = data_Class,
                            annotation = data.frame(group = datagroup, row.names = colnames(data_Class)), 
                            fontsize_col = 20, 
                            fontsize_row = 20, 
                            fontsize = 12,
                            clustering_distance_rows = "euclidean", 
                            clustering_distance_cols = "euclidean", 
                            clustering_methods = "ward.D", 
                            cluster_rows = T, 
                            cluster_cols = F, 
                            scale = "row", 
                            annotation_colors = list(group = colorpars))
    pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                           "_vs_", controlGrp))
    lipnames <- rownames(data_Class)
    maxlen <- max(nchar(lipnames))
    pdf(paste0(fileLoc, "heatmap_lipClassSummary_", pname, ".pdf"), 
        width=8/15.3*(10.5/length(allgroups)*15+4.8), 
        #height=8/15.3*(12.6/nclass*10+2.7)
        height=20/15.3*(12.6/nclass*10+2.7)
        )
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }else{
    data_Class_det <- data_tidy %>%
      arrange(Class) %>%
      column_to_rownames(var = "lipidName")
    data_Class_det <- data_Class_det[apply(data_Class_det, 1, function(x) ifelse(is.na(sd(x[-length(x)], na.rm = T)), F, ifelse(sd(x[-length(x)], na.rm = T) == 0, F, T))), ]
    if(nrow(data_Class_det) < nrow(data_tidy)){
      cat("Some lipids in low variation will not show in the plot.\n")
    }
    ClassInfo <- data_Class_det$Class
    data_Class_det <- select(data_Class_det, -Class)
    data_Class_det <- data_Class_det[, match(colnames(data_Class_det), names(allgroups))]
    colorpars2 <- plottingPalettes(n = length(unique(ClassInfo))+length(groupsLevel), type = "discrete")[-1:-length(groupsLevel)]
    names(colorpars2) <- unique(ClassInfo)
    
    x <- pheatmap::pheatmap(mat = data_Class_det,
                            annotation_col = data.frame(group = datagroup, row.names = colnames(data_Class_det)), 
                            annotation_row = data.frame(lipidClass = ClassInfo, row.names = rownames(data_Class_det)),
                            # fontsize_col = 20, 
                            # fontsize_row = 20, 
                            #fontsize = 20,
                            clustering_distance_rows = "euclidean", 
                            clustering_distance_cols = "euclidean", 
                            clustering_methods = "ward.D", 
                            cluster_rows = F, 
                            cluster_cols = F, 
                            scale = "row", 
                            show_rownames = F, 
                            annotation_colors = list(group = colorpars, lipidClass = colorpars2))
    pdf(paste0(fileLoc, "heatmap_lipClassSummary_detail_", pname, ".pdf"), width=4, height=7)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
}

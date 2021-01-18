lipSumChainHeatmapPlot <- function(dataSet, mSet, 
                                   plotInfo = "FA_info", ignore = T, show_detail = F,
                                   fileLoc){
  
  allgroups <- dataSet$allgroups
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  dataType <- dataSet$dataType
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  
  datagroup <- factor(allgroups, 
                      levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
  colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
  names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
  
  source("./utilityFunc/FAchainStat.R")
  lipid_subclass_stat_output <- FAchainStat(dataSet=dataSet, mSet=mSet,  
                                            fileLoc = NULL, plotInfo = plotInfo, ignore = ignore,
                                            stat = F, stat2 = T)
  data_FAchain <- lipid_subclass_stat_output %>%
    mutate(chain = gsub(".*?([0-9]+):.*", "C\\1", subclass), 
           #unsaturate = as.numeric(gsub(".*?:([0-9]+).*", "\\1", subclass)), 
           Class = gsub("\\(.*\\)$", "", subclass)) %>%
    ungroup() %>%
    select(-subclass)
  data_FAchain_split <- split(data_FAchain, f = data_FAchain$Class)
  for(i in 1:length(data_FAchain_split)){
    data_i <- data_FAchain_split[[i]] %>%
      arrange(chain) %>%
      select(-Class) %>%
      group_by(chain) %>%
      mutate(id = paste0(chain, "_", 1:n())) %>%
      column_to_rownames(var = "id")
    
    if(show_detail == F){
      data_i <- data_i %>%
        gather(-chain, key = "group", value = "conc") %>% 
        group_by(group, chain) %>%
        summarise(lipsum = sum(conc, na.rm = T)) %>%
        spread(group, lipsum) %>%
        column_to_rownames(var = "chain")
      data_i <- data_i[, match(colnames(data_i), names(allgroups))]
      data_i <- data_i[apply(data_i, 1, function(x) sd(x)!=0), ]
      # if(nrow(data_Class) < nrow(data_Class)){
      #   cat("Some lipids in low variation will not show in the plot.\n")
      # }
      x <- pheatmap::pheatmap(mat = data_i,
                              annotation = data.frame(group = datagroup, row.names = colnames(data_i)), 
                              fontsize_col = 20, 
                              fontsize_row = 20, 
                              fontsize = 10,
                              main = names(data_FAchain_split)[i],
                              clustering_distance_rows = "euclidean", 
                              clustering_distance_cols = "euclidean", 
                              clustering_methods = "ward.D", 
                              cluster_rows = T, 
                              cluster_cols = F, 
                              scale = "row", 
                              annotation_colors = list(group = colorpars))
      lipnames <- rownames(data_i)
      maxlen <- max(nchar(lipnames))
      pdf(paste0(fileLoc, "heatmap_lipChainSummary_", names(data_FAchain_split)[i], "_", pname, ".pdf"), width=8, height=4)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
    
    
    
    
    else{
      data_i <- data_i[apply(data_i, 1, function(x) ifelse(is.na(sd(x[-length(x)], na.rm = T)), F, ifelse(sd(x[-length(x)], na.rm = T) == 0, F, T))), ]
      # if(nrow(data_Class_det) < nrow(data_tidy)){
      #   cat("Some lipids in low variation will not show in the plot.\n")
      # }
      chainInfo <- data_i$chain
      data_i <- select(data_i, -chain)
      data_i <- data_i[, match(colnames(data_i), names(allgroups))]
      colorpars2 <- plottingPalettes(n = length(unique(chainInfo))+length(groupsLevel), type = "discrete")[-1:-length(groupsLevel)]
      names(colorpars2) <- unique(chainInfo)
      x <- pheatmap::pheatmap(mat = data_i,
                              annotation_col = data.frame(group = datagroup, row.names = colnames(data_i)), 
                              annotation_row = data.frame(lipidChain = factor(chainInfo), row.names = rownames(data_i)),
                              # fontsize_col = 20, 
                              # fontsize_row = 20, 
                              fontsize = 20,
                              clustering_distance_rows = "euclidean", 
                              clustering_distance_cols = "euclidean", 
                              clustering_methods = "ward.D", 
                              cluster_rows = F, 
                              cluster_cols = F, 
                              scale = "row", 
                              show_rownames = F, 
                              annotation_colors = list(group = colorpars, lipidChain = colorpars2), 
                              main = names(data_FAchain_split)[i])
      pdf(paste0(fileLoc, "heatmap_lipChainSummary_", names(data_FAchain_split)[i], "_detail_", pname, ".pdf"), width=7, height=7)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
  }
}

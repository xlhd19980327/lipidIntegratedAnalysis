rnaMiArPCAPlot <- function(dataProc = dataProc_RNA, 
                           fileLoc){
  eset <- dataProc$eset
  controlGrp <- dataProc$dataSet$controlGrp
  groupsLevel <- dataProc$dataSet$groupsLevel
  
  ## PCA & plot
  pcaAnal <- prcomp(t(eset), center = TRUE, scale = F)
  pclevel <- c(1, 2)
  pcaScore <- as.data.frame(pcaAnal$x)[, pclevel] %>%
    rownames_to_column(var = "sampleLabel") 
  pcaScore <- pcaScore %>%
    mutate(group = factor(dataProc$dataSet$sampleInfo$conditions[match(pcaScore$sampleLabel, dataProc$dataSet$sampleInfo$samples)], 
                          levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))
  pcaVar <- summary(pcaAnal)$importance[2, ]
  ellipse_coor_all <- data.frame()
  for(i in groupsLevel){
    groupVar <- var(subset(pcaScore, subset = group == i)[, 2:3], na.rm = T)
    groupMean <- cbind(mean(subset(pcaScore, subset = group == i)[, 2], na.rm = T), 
                       mean(subset(pcaScore, subset = group == i)[, 3], na.rm = T))
    ellipse_coor <- as.data.frame(ellipse::ellipse(groupVar, centre = groupMean, 
                                                   level = 0.95, npoints = 100))
    colnames(ellipse_coor) <- c("x", "y")
    ellipse_coor_all <- rbind(ellipse_coor_all, 
                              cbind(ellipse_coor, group = i))
  }
  colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
  colorpars_fill <- alpha(colorpars, 0.3)
  names(colorpars_fill) <- names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
  pcaScorePlot <- ggplot() +
    geom_polygon(ellipse_coor_all, mapping = aes(x = x, y = y, 
                                                 color = factor(group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])), 
                                                 fill = factor(group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))) +
    geom_point(pcaScore, mapping =aes(x = PC1, y = PC2, color = group)) +
    scale_color_manual(values = colorpars) +
    scale_fill_manual(values = colorpars_fill) +
    theme_bw() + 
    labs(x = paste0("PC", pclevel[1], "(", round(100 * pcaVar[pclevel[1]], 1), "%)"), 
         y = paste0("PC", pclevel[2], "(", round(100 * pcaVar[pclevel[2]], 1), "%)"), 
         color = "group", 
         fill = "group",
         title = "PCA Score Plot") +
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 15), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  ggsave(paste0(fileLoc, "PCA_score_plot_", pname, ".pdf"), plot = pcaScorePlot, 
         device = "pdf", width = 9, height = 9)
}
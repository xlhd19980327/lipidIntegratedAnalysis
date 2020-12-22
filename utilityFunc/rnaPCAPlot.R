rnaPCAPlot <- function(dataProc, fileLoc){
  rld <- dataProc$rld
  groupsLevel <- dataProc$dataSet$groupsLevel
  controlGrp <- dataProc$dataSet$controlGrp
  batremvOpt <- ifelse("batch" %in% colnames(dataProc$dataSet$sampleInfo), T, F)
  
  ##use VST or rlog data to do PCA(DESeq2 recommanded, here use rlog)
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
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 15), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))
  ggsave(paste0(fileLoc, "PCA_score_plot_all.pdf"), plot = pcaScorePlot, 
         device = "pdf", width = 9, height = 9)
}

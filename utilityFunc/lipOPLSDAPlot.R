lipOPLSDAPlot <- function(dataSet, mSet,  
                          fileLoc){
  groupsLevel <- dataSet$groupsLevel
  controlGrp <- dataSet$controlGrp
  data_norm <- mSet$dataSet$norm
  
  cls <- scale(as.numeric(mSet$dataSet$cls))[, 1]
  datmat <- as.matrix(mSet$dataSet$norm)
  cv.num <- min(7, dim(mSet$dataSet$norm)[1] - 1)
  res <- MetaboAnalystR:::perform_opls(datmat, cls, predI = 1, 
                                       permI = 0, orthoI = NA, crossvalI = cv.num)
  
  lv1 <- res[["scoreMN"]][, 1]
  lv2 <- res[["orthoScoreMN"]][, 1]
  xlabel <- paste0("T score[1]", "(", round(100 * res$modelDF["p1", "R2X"], 1), "%)")
  ylabel <- paste0("Ortho-T score[1]", "(", round(100 * res$modelDF["o1", "R2X"], 1), "%)")
  text.lbls <- substr(rownames(data_norm), 1, 12)
  oplsdaScoreScore <- cbind(lv1, lv2) %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleLabel") %>%
    mutate(group = factor(mSet$dataSet$cls, 
                          levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))
  
  ellipse_coor_all <- data.frame()
  for (i in groupsLevel) {
    groupVar <- var(subset(oplsdaScoreScore, subset = group == i)[, 2:3], na.rm = T)
    groupMean <- cbind(mean(subset(oplsdaScoreScore, subset = group == i)[, 2], na.rm = T), 
                       mean(subset(oplsdaScoreScore, subset = group == i)[, 3], na.rm = T))
    ellipse_coor <- as.data.frame(ellipse::ellipse(groupVar, centre = groupMean, 
                                                   level = 0.95, npoints = 100))
    colnames(ellipse_coor) <- c("x", "y")
    ellipse_coor_all <- rbind(ellipse_coor_all, 
                              cbind(ellipse_coor, group = i))
  }
  
  colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
  colorpars_fill <- alpha(colorpars, 0.3)
  names(colorpars_fill) <- names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
  
  oplsdaScorePlot <- ggplot() +
    geom_polygon(ellipse_coor_all, mapping = aes(x = x, y = y, 
                                                 color = factor(group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])), 
                                                 fill = factor(group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))) +
    geom_point(oplsdaScoreScore, mapping =aes(x = lv1, y = lv2, color = group)) +
    scale_color_manual(values = colorpars) +
    scale_fill_manual(values = colorpars_fill) +
    theme_bw() + 
    labs(x = xlabel, 
         y = ylabel, 
         color = "group", 
         fill = "group",
         title = "OPLS-DA Score Plot") +
    theme(legend.position = "right", 
          plot.title = element_text(hjust = 0.5, size = 25, face = "bold"), 
          axis.title = element_text(size = 25), 
          axis.text = element_text(size = 25),
          legend.text = element_text(size = 25), 
          legend.title = element_text(size = 25)) 
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  ggsave(paste0(fileLoc, "OPLSDA_score_plot_", pname, ".pdf"), plot = oplsdaScorePlot, 
         device = "pdf", width = 9, height = 9)
}
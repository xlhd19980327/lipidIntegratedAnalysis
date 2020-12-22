lipVolcanoPlot <- function(dataSet, mSet, showLipClass = F, 
                           fileLoc,
                           paired = F, pval.type = "raw", fcthresh = 2.0, pthresh = 0.1, 
                           showtop = 10, stat = F, ignore = T){
  source("./utilityFunc/plottingPalettes.R")
  ## Source will offer the following contents:
  ## Function(s): getClassInfo
  source("./utilityFunc/getClassInfo.R")
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  data_norm <- mSet$dataSet$norm
  #MetaboAnalystR no longer use mSet$dataSet$row.norm
  #data_norm2 <- mSet$dataSet$row.norm
  data_norm2 <- qs::qread("row_norm.qs")
  
  ## Volcano analysis & plot
  #paired <- FALSE
  equal.var <- TRUE
  #!!!Client options: "raw" or "fdr"
  #pval.type <- "raw"
  cprind <- mSet[["dataSet"]][["cls"]]  == controlGrp
  exper <- groupsLevel[groupsLevel != controlGrp]
  expind <- mSet[["dataSet"]][["cls"]]  == exper
  p.value <- apply(as.matrix(data_norm), 2, 
                   function(x) {
                     tmp <- tryCatch(t.test(x[cprind], x[expind], paired = paired, 
                                            var.equal = equal.var), 
                                     error = function(e){cat("Note: A NA generate!\n")}
                     )
                     if (is.null(tmp)) {
                       return(NA)
                     }
                     else {
                       return(tmp$p.value)
                     }
                   })
  if(pval.type == "fdr"){
    p.value <- p.adjust(p.value, "fdr")
  }
  cprmean <- colMeans(data_norm2[cprind, ])
  expmean <- colMeans(data_norm2[expind, ])
  ratio <- expmean/cprmean
  fc.log <- signif(log2(ratio), 5)
  if(stat){
    result <- as.data.frame(cbind(p.value, fc.log)) %>%
      rownames_to_column(var = "lipid") 
    return(result)
  }
  p.log <- -log10(p.value)
  result <- cbind(p.log, fc.log)
  write.csv(result, paste0(fileLoc, "volcano_data.csv"))
  
  
  #!!!Client options: Fold change threshold
  #fcthresh <- 2.0
  #!!!Client options: Fold change p.value
  #pthresh <- 0.1
  volcano.data <- result %>%
    as.data.frame() %>%
    rownames_to_column(var = "lipid") %>%
    mutate(
      up = ifelse(fc.log > log(fcthresh, 2) & p.log > -log10(pthresh), 1, 0),
      down = ifelse(fc.log < -log(fcthresh, 2) & p.log > -log10(pthresh), 2, 0),
      regState = sapply(up + down, function(x) switch(x, "upreg", "downreg")), 
      regState = sapply(regState, function(x) ifelse(is.null(x), 0, x))
    )
  
  volcano.data_reg <- subset(volcano.data, subset = regState != 0)
  volcano.data_unreg <- subset(volcano.data, subset = regState == 0)
  if(showtop != 0){
    if(showtop > nrow(volcano.data_reg)){
      cat("No enough significant features! Show all the feature labels.\n")
    }
    volcano.data_top <- volcano.data_reg[order(-abs(volcano.data_reg$fc.log)), ][1:showtop, ]
  }else{
    volcano.data_top <- volcano.data_reg[order(-abs(volcano.data_reg$fc.log)), ][0, ]
  }
  volcano.plot <- ggplot() +
    geom_point(data = subset(volcano.data, subset = regState != 0), 
               aes(x = fc.log, y = p.log, color = factor(regState)))+ 
    geom_point(data = subset(volcano.data, subset = regState == 0), 
               aes(x = fc.log, y = p.log), color = "gray") +
    geom_text_repel(data = volcano.data_top, mapping = aes(x = fc.log, y = p.log, label = lipid), 
                    hjust = 0, nudge_x = 0.05) +
    theme_bw() +
    geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
    labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "regulation", 
         title = paste0(exper, " vs ", controlGrp)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 15), 
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))+ 
    scale_color_aaas()
  ggsave(paste0(fileLoc, "volcano_reg_", exper, ".pdf"), plot = volcano.plot, 
         device = "pdf", width = 9, height = 9)
  
  if(showLipClass){
    volcano.data <- volcano.data %>%
      mutate(
        Class = switch(dataSet$dataType,
                       LipidSearch = sapply(lipid, getClassInfo, "LipidSearch", ignore = ignore), 
                       MS_DIAL = sapply(lipid, getClassInfo, "MS_DIAL", ignore = ignore))
      )
    nClass <- length(unique(volcano.data$Class))
    volcano.plot_class <- ggplot() +
      geom_point(data = subset(volcano.data, subset = regState != 0), 
                 aes(x = fc.log, y = p.log, color = reorder(Class, -abs(fc.log), mean)))+ 
      geom_point(data = subset(volcano.data, subset = regState == 0), 
                 aes(x = fc.log, y = p.log), color = "gray") +
      theme_bw() +
      geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
      geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
      labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "Class", 
           title = paste0(exper, " vs ", controlGrp)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 20), 
            axis.title = element_text(size = 15), 
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12)) +
      scale_color_manual(values = plottingPalettes(nClass, type = "discrete"))
    ggsave(paste0(fileLoc, "volcano_regClass_", exper, ".pdf"), plot = volcano.plot_class, 
           device = "pdf", width = 9, height = 9)
  }
  
}

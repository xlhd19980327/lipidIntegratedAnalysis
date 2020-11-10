rnaVolcanoPlot <- function(DEAresult = DEAresult, type = "RNAseq",
                           fileLoc, fcthresh = 2.0, pthresh = 0.1,
                           #Show top gene labels in volcano or set "F" to do not show any
                           showtop = 20){
  resLFC <- DEAresult$resLFC
  experGrp <- DEAresult$experGrp
  controlGrp <- DEAresult$dataProc$dataSet$controlGrp
  
  if(type == "MiAr"){
    resLFC <- resLFC %>%
      dplyr::rename(log2FoldChange = logFC, 
                    padj = adj.P.Val) 
  }
  
  volcano.data <- resLFC %>%
    mutate(
      log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
      padj = ifelse(is.na(padj), 1, padj),
      up = ifelse(log2FoldChange > log(fcthresh, 2) & padj < pthresh, 1, 0),
      down = ifelse(log2FoldChange < -log(fcthresh, 2) & padj < pthresh, 2, 0),
      regState = sapply(up + down, function(x) switch(x, "upreg", "downreg")), 
      regState = sapply(regState, function(x) ifelse(is.null(x), 0, x))
    )
  volcano.data_reg <- subset(volcano.data, subset = regState != 0)
  volcano.data_unreg <- subset(volcano.data, subset = regState == 0)
  #!!!!!WARNING: use p_value as rank standard, may use others(eg. mad?)
  if(showtop != F){
    if(showtop > nrow(volcano.data_reg)){
      cat("Not enough significant feature! Show all the feature label.\n")
    }
    volcano.data_top <- volcano.data_reg[order(volcano.data_reg$padj), ][1:showtop, ]
  }else{
    volcano.data_top <- volcano.data_reg[order(volcano.data_reg$padj), ][0, ]
  }
  volcano.plot <- ggplot(mapping = aes(x = log2FoldChange, y = -log(padj, 10))) +
    geom_point(data = volcano.data_reg, mapping = aes(color = factor(regState))) + 
    geom_point(data = volcano.data_unreg, color = "gray") +
    geom_text_repel(data = volcano.data_top, mapping = aes(label = gene), 
                    hjust = 0, nudge_x = 0.05) +
    theme_bw() +
    geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
    labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "regulation", 
         title = paste0(experGrp, " vs ", controlGrp)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))+ 
    scale_color_aaas()
  ggsave(paste0(fileLoc, "volcano_", experGrp, "_vs_", controlGrp, ".pdf"), plot = volcano.plot, 
         device = "pdf", width = 9, height = 9)
}

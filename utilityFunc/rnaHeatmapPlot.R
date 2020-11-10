rnaHeatmapPlot <- function(DEAresult, showtop = 75, showallgroups = F, type = "RNAseq", 
                           fileLoc){
  if(showallgroups == F){
    controlGrp <- DEAresult$dataProc$dataSet$controlGrp
    sampleInfo <- DEAresult$dataProc$dataSet$sampleInfo
    groupsLevel <- DEAresult$dataProc$dataSet$groupsLevel
    experGrp <- DEAresult$experGrp
    
    colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
    names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
    legend_annotation <- data.frame(conditions = sampleInfo$conditions)
    rownames(legend_annotation) <- sampleInfo$samples
    
    if(type == "RNAseq"){
      resLFC <- DEAresult$resLFC
      normalized_counts <- DEAresult$dataProc$normalized_counts
      dds <- DEAresult$dataProc$dds
      
      volcano.data <- resLFC %>%
        mutate(padj = ifelse(is.na(padj), 1, padj))
      if(showtop > nrow(volcano.data)){
        cat("Not enough significant feature! Show all the feature label.\n")
      }
      #!!!!!WARNING: use p_value as rank standard, may use others(eg. mad?)
      heatmap.data_topgene <- volcano.data[order(volcano.data$padj), ][1:showtop, ]$gene
      heatmap.data_top <- normalized_counts[match(heatmap.data_topgene, rownames(normalized_counts)), 
                                            colData(dds)$conditions %in% c(controlGrp, experGrp)]
    }
    if(type == "MiAr"){
      fit <- DEAresult$fit
      resLFC <- topTable(fit, adjust="BH", number = showtop)
      normalized_counts <- DEAresult[["dataProc"]][["eset"]]
      
      heatmap.data_topgene <- rownames(resLFC)
      heatmap.data_top <- normalized_counts[match(heatmap.data_topgene, rownames(normalized_counts)), 
                                            sampleInfo$conditions %in% c(controlGrp, experGrp)]
      
    }
    collegend_annotation <- list(conditions = colorpars[names(colorpars) %in% c(experGrp, controlGrp)])
  } else{
    #Use dataProc to do overall heatmap
    controlGrp <- DEAresult$dataSet$controlGrp
    sampleInfo <- DEAresult$dataSet$sampleInfo
    groupsLevel <- DEAresult$dataSet$groupsLevel
    
    colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
    names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
    legend_annotation <- data.frame(conditions = sampleInfo$conditions)
    rownames(legend_annotation) <- sampleInfo$samples
    
    if(type == "RNAseq"){
      normalized_counts <- DEAresult$normalized_counts
    }
    if(type == "MiAr"){
      normalized_counts <- DEAresult$eset
    }
    
    ##We really donot recommand client to use showtop
    #heatmap.data_top <- normalized_counts[1:showtop, ]
    heatmap.data_top <- normalized_counts
    collegend_annotation <- list(conditions = colorpars)
  }

  tidy_ind <- match(sampleInfo$samples[order(factor(sampleInfo$conditions, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp])))], 
                    colnames(heatmap.data_top))
  tidy_ind <- tidy_ind[!is.na(tidy_ind)]
  heatmap.data_top <- heatmap.data_top[, tidy_ind]
  
  ntop <- nrow(heatmap.data_top)
  heatmap.data_top <- heatmap.data_top[apply(heatmap.data_top, 1, function(x) sd(x)!=0), ]
  if(nrow(heatmap.data_top) < ntop){
    cat("Some genes in low variation will be removed.\n")
  }
  if(showallgroups){
    heatmap <- pheatmap::pheatmap(mat = heatmap.data_top, 
                                  annotation = legend_annotation, 
                                  color = plottingPalettes(100, type = "continuous"),
                                  #one-by-one usage
                                  annotation_colors = collegend_annotation,
                                  cluster_rows = T, 
                                  cluster_cols = T, 
                                  scale = "row"
    )
    pdf(paste0(fileLoc, "heatmap_allgroups.pdf"), 
        width=7, 
        height=20)
  }else{
    heatmap <- pheatmap::pheatmap(mat = heatmap.data_top, 
                                  annotation = legend_annotation, 
                                  color = plottingPalettes(100, type = "continuous"),
                                  #one-by-one usage
                                  annotation_colors = collegend_annotation,
                                  fontsize = 8, 
                                  fontsize_row = 8, 
                                  cluster_rows = T, 
                                  cluster_cols = F, 
                                  scale = "row"
    )
    pdf(paste0(fileLoc, "heatmap_top", showtop, 
               ifelse(showallgroups, "_allgroups", 
                      paste0("_", experGrp, "_vs_", controlGrp)), ".pdf"), 
        width=(ncol(heatmap.data_top)*25+300)/72, height=(nrow(heatmap.data_top)*18+150)/72)
  }
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  dev.off()
}

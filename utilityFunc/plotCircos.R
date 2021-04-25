plotCircos <- function(gofileLoc, corfileLoc, k, j, 
                       thresh1 = 0.6, topnum = 20, outfileLoc){
  #k: gene, j: met
  filenames <- paste0(corfileLoc, "correlation_", j, "_", k, ".csv")
  cordata <- read.csv(filenames, row.names = 1)
  lipids <- rownames(cordata)
  genes <- colnames(cordata)
  load(paste0(gofileLoc, "data_circos.RData"))
  godata <- go@result
  #thresh1 <- 0.6
  subdatas <- lapply(cordata, function(x) lipids[abs(x) > thresh1])
  subcors <- lapply(cordata, function(x) x[abs(x) > thresh1])
  thresh2 <- 0.05
  #topnum <- 20
  subgodata <- subset(godata, p.adjust < thresh2)
  if(nrow(subgodata) < topnum){
    topnum <- nrow(subgodata)
    cat("Enrichment item number may not fit your demand. Show all the enrichment items in the plot.")
  }
  subgodata <- subgodata[order(subgodata$p.adjust), ][1:topnum, ]
  subgogenes <- strsplit(subgodata$geneID, "/")
  
  library(paletteer)
  library(circlize)
  gene_corfilt <- unlist(mapply(rep, names(subcors), ifelse(sapply(subcors, length) == 0, 0, 1)))
  gene_gofilt <- unlist(subgogenes)
  ind1 <- gene_corfilt %in% gene_gofilt
  if(all(!ind1)){
    cat("Please check your correlation threshold. It may be too strict to filter.\n")
    return(1)
  }
  gene_filt <- gene_corfilt[gene_corfilt %in% gene_gofilt]
  lipid_filt <- subdatas[gene_filt]
  gl_link <- mapply(cbind, names(lipid_filt), lipid_filt)
  gl_link <- do.call(rbind, gl_link)
  gl_link_ind <- data.frame(
    geneind = match(gl_link[, 1], gene_filt), 
    lipind = match(gl_link[, 2], unique(unlist(lipid_filt))),
    value = mapply(function(x, y){cordata[x, y]}, gl_link[, 2], gl_link[, 1])
  )
  colind1 <- match(gl_link_ind$value, sort(unique(gl_link_ind$value)))
  col1 <- as.character(paletteer_c("ggthemes::Classic Area Green", n = length(unique(gl_link_ind$value))))
  #col1 <- gsub("FF", "", col1)
  gl_link_ind <- cbind(gl_link_ind, 
                       col = col1[colind1])
  go_link <- lapply(subgogenes, function(x){
    ind <- match(x, gene_filt)
    ind <- ind[!is.na(ind)]
    return(ind)
  })
  subgodata_item <- subgodata$Description
  subgo_filt <- subgodata_item[unlist(mapply(rep, 1:length(subgodata_item), sapply(go_link, function(x) ifelse(length(x) == 0, 0, 1))))]
  indtemp <- sapply(go_link, length)
  indtemp <- indtemp[indtemp != 0]
  go_link_ind <- data.frame(
    goind = unlist(mapply(rep, 1:length(subgo_filt), indtemp)), 
    geneind = unlist(go_link)
  )
  go_link_ind <- cbind(go_link_ind, 
                       value = subgodata$p.adjust[go_link_ind$goind])
  colind2 <- match(go_link_ind$value, sort(unique(go_link_ind$value), decreasing = T))
  col2 <- as.character(paletteer_c("ggthemes::Classic Area Red", n = length(unique(go_link_ind$value))))
  #col2 <- gsub("FF", "", col2)
  go_link_ind <- cbind(go_link_ind, 
                       col = col2[colind2])
  
  sectors <- c("Genes", "Lipids", "Gene_Ontology")
  n1 <- length(gene_filt)
  n2 <- length(unique(unlist(lipid_filt)))
  n3 <- length(subgo_filt)
  df <- data.frame(
    sectors = sectors, 
    xlim = c(n1, n2, n3)
  )
  circos.clear()
  pdf(paste0(outfileLoc, "circosPlot.pdf"), width = 8)
  circos.initialize(factor(df$sectors, levels = df$sectors), 
                    xlim = matrix(c(rep(0,3),df$xlim),ncol=2))
  circos.track(ylim=c(0,1),panel.fun=function(x,y) {
    item=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
    circos.text(mean(xlim),mean(ylim),item,cex=0.9,col="grey40",
                facing="bending.inside",niceFacing=TRUE, font = 2)
  },bg.col="grey90",bg.border=F,track.height=0.08)
  circos.par("track.margin" = c(0, 0.1))
  circos.par("cell.padding" = c(0, 1, 0, 1))
  circos.track(ylim = c(0,10), panel.fun=function(x,y) {
    features <- CELL_META$xrange
    data <- switch(CELL_META$sector.index, 
                   Genes = gene_filt, 
                   Lipids = unique(unlist(lipid_filt)), 
                   Gene_Ontology = subgodata_item)
    for(i in 1:features){
      circos.text(i, 0, data[i], facing = "clockwise", adj = c(0, 0.5),
                  cex = 0.5, niceFacing = T, font = 2)
    }
  }, bg.border=F)
  circos.par("track.margin" = c(0, 0))
  for(i in 1:nrow(gl_link_ind)){
    circos.link("Genes", gl_link_ind$geneind[i], "Lipids", gl_link_ind$lipind[i], 
                col = gl_link_ind$col[i])
  }
  for(i in 1:nrow(go_link_ind)){
    circos.link("Gene_Ontology", go_link_ind$goind[i], "Genes", go_link_ind$geneind[i], 
                col = go_link_ind$col[i])
  }
  dev.off()
  
}

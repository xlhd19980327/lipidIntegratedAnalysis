## Do correlation
source("./branch/correlation/readingLipidData_cor2.R")
source("./branch/correlation/readingRNAData_cor.R")

library(pheatmap)
library(tidyverse)
library(MetaboAnalystR)
library(DESeq2)
library(limma)

library(dbscan)
library(MCL)
library(scales)
library(reshape)
options(stringsAsFactors = F)

#Shell interaction
#library(getopt)
library(optparse)

option_list <- list( 
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  make_option(c("-t", "--data_type"), action="store"), 
  #make_option(c("-f", "--feature_field"), action="store", default = NA, type = "character"), 
  #make_option(c("-n", "--NA_string"), action="store", default = NULL, type = "character"), 
  make_option(c("-l", "--lipid_odd_chain_deletion"), action="store", default = F),
  make_option(c("-m", "--na_percent"), action="store", default = 67),
  
  make_option(c("-j", "--input_file_rna"), action="store"),
  make_option(c("-e", "--description_file_rna"), action="store"), 
  make_option(c("-u", "--data_type_rna"), action="store"), 
  
  make_option(c("-g", "--gene_split"), action="store", default = 4), 
  make_option(c("-k", "--lipid_split"), action="store", default = 4),
  make_option(c("-n", "--use_quantile"), action="store", default = T), 
  make_option(c("-s", "--max_thresh"), action="store", default = 0.8), 
  make_option(c("-o", "--output"), action="store"), 
  
  make_option(c("-b", "--clustering_type"), action="store", default = "hierarchical"), 
  make_option(c("-c", "--dbscan_l_pts"), action="store", default = 4), 
  make_option(c("-f", "--dbscan_g_pts"), action="store", default = 3), 
  make_option(c("-p", "--markov_quantile"), action="store", default = 0.55)
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Data input
lipid_data <- readingLipidData_cor(datafile = opt$input_file, sampleList = opt$description_file,
                               dataType = opt$data_type, delOddChainOpt = opt$lipid_odd_chain_deletion, 
                               perc = opt$na_percent)
if(opt$data_type_rna == "Proteins"){
  gene_data <- readingLipidData_cor(datafile = opt$input_file_rna, sampleList = opt$description_file_rna,
                                dataType = "Proteins", delOddChainOpt = F, 
                                perc = 100)
  gene_data <- t(gene_data)
}else{
  gene_data <- readingRNAData(datafile = opt$input_file_rna, 
                            sampleList = opt$description_file_rna,
                            type = opt$data_type_rna)
}


##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

### Waiting time(~2+ min, larger data may require more)
correlation <- cor(lipid_data, t(gene_data), method ="spearman")
correlation[is.na(correlation)] <- 0
#correlation <- correlation[,which(apply(correlation, 2, FUN = sd) != 0)]
#correlation <- correlation[which(apply(correlation, 1, FUN = sd) != 0),]

cutoff <- function(x){
  x<-abs(x)
  x<-max(x)
}
if(opt$use_quantile){
  #70% quantile
  max_list <- apply(correlation,2,max)
  value <- quantile(max_list,0.7)
}else{
  #cor thresh
  value <- opt$max_thresh
}

correlation <- correlation[,which(apply(correlation, 2, cutoff)>value)]## set a parameter for cutoff, default =0.8
correlation <- correlation[which(apply(correlation, 1, cutoff)>value),]
write.csv(correlation, paste0(opt$output, "correlation.csv"))

if(opt$clustering_type == "hierarchical"){
  pdf(paste0(opt$output, "correlationPlot.pdf"))
  ### Waiting time(~2+ min, larger data may require more)
  list <- pheatmap(correlation, 
                   cutree_col = opt$gene_split, cutree_rows = opt$lipid_split)
  dev.off()
  
  cutgene <- cutree(list$tree_col, k = opt$gene_split)[list[["tree_col"]][["order"]]]
  cutlipid <- cutree(list$tree_row, k = opt$lipid_split)[list[["tree_row"]][["order"]]]
  genegaps <- c(names(cutgene[1]), 
                names(which((cutgene[-1] - cutgene[-length(cutgene)]) != 0)))
  lipidgaps <- c(names(cutlipid[1]), 
                 names(which((cutlipid[-1] - cutlipid[-length(cutlipid)]) != 0)))
  for(i in 1:opt$lipid_split){
    for(j in 1:opt$gene_split){
      subdata <- correlation[names(cutlipid)[cutlipid == i], names(cutgene)[cutgene == j], drop = F]
      ith <- which(lipidgaps %in% rownames(subdata))
      jth <- which(genegaps %in% colnames(subdata))
      sublipids <- rownames(subdata)
      subgenes <- colnames(subdata)
      write.csv(subdata, paste0(opt$output, "correlation_", ith, "_", jth, ".csv"))
      write.csv(sublipids, paste0(opt$output, "lipids_", ith, ".csv"), row.names = F)
      write.csv(subgenes, paste0(opt$output, "genes_", jth, ".csv"), row.names = F)
    }
  }
}
if(opt$clustering_type == "k_means"){
  kl <- opt$lipid_split
  kg <- opt$gene_split
  kmeans.lipid.cluster <- kmeans(correlation, kl)
  kmeans.gene.cluster <- kmeans(t(correlation), kg)
  
  cluster.l <- data.frame(kmeans.lipid.cluster[["cluster"]])
  cluster.l[,'lipid'] <- rownames(cluster.l)
  cluster.l <- cluster.l[order(cluster.l$kmeans.lipid.cluster...cluster...),]
  colnames(cluster.l) <- c('cluster.l', 'lipid')
  
  cluster.g <- data.frame(kmeans.gene.cluster[["cluster"]])
  cluster.g[,'gene'] <- rownames(cluster.g)
  cluster.g <- cluster.g[order(cluster.g$kmeans.gene.cluster...cluster...),]
  colnames(cluster.g) <- c('cluster.g', 'gene')
  
  mat <- correlation[order(kmeans.lipid.cluster[["cluster"]]),order(kmeans.gene.cluster[["cluster"]])]
  mat.m <- melt(mat)
  colnames(mat.m) <- c('lipid', 'gene', 'value')
  mat.m <- merge(mat.m, cluster.g, by = 'gene')
  mat.m <- merge(mat.m, cluster.l, by = 'lipid')
}
if(opt$clustering_type == "DBSCAN"){
  MinPts.l <- opt$dbscan_l_pts
  MinPts.g <- opt$dbscan_g_pts
  Min <- min(correlation)
  Max <- max(correlation)
  Rank <- Max - Min
  correlation_rank <- (correlation-Min)/Rank
  
  dbscan.lipid.cluster <- hdbscan(correlation_rank, minPts = MinPts.l)
  dbscan.gene.cluster <- hdbscan(t(correlation_rank), minPts = MinPts.g)
  
  cluster.l <- data.frame(dbscan.lipid.cluster[["cluster"]])
  cluster.l[,'lipid'] <- rownames(correlation_rank)
  cluster.l <- cluster.l[order(cluster.l$dbscan.lipid.cluster...cluster...),]
  colnames(cluster.l) <- c('cluster.l', 'lipid')
  
  cluster.g <- data.frame(dbscan.gene.cluster[["cluster"]])
  cluster.g[,'gene'] <- colnames(correlation_rank)
  cluster.g <- cluster.g[order(cluster.g$dbscan.gene.cluster...cluster...),]
  colnames(cluster.g) <- c('cluster.g', 'gene')
  
  mat <- correlation_rank[order(cluster.l$cluster.l),order(cluster.g$cluster.g)]
  mat.m <- melt(mat)
  colnames(mat.m) <- c('lipid', 'gene', 'value')
  mat.m <- merge(mat.m, cluster.g, by = 'gene')
  mat.m <- merge(mat.m, cluster.l, by = 'lipid')
}
if(opt$clustering_type == "MCL"){
  lipid_cor <- cor(t(correlation))
  lipid_cor <- (lipid_cor - min(lipid_cor))/(max(lipid_cor) - min(lipid_cor))
  lipid_cut <- quantile(lipid_cor,opt$markov_quantile)
  lipid_cor[lipid_cor<=lipid_cut] <- 0
  
  gene_cor <- cor(correlation)
  gene_cor <- (gene_cor - min(gene_cor))/(max(gene_cor) - min(gene_cor))
  gene_cut <- quantile(gene_cor,opt$markov_quantile)
  gene_cor[gene_cor<=gene_cut] <- 0
  mcl.lipid.cluster <- mcl(lipid_cor, addLoops= T, allow1 = T)
  mcl.gene.cluster <- mcl(gene_cor, addLoops = T, allow1 = T)
  
  cluster.l <- data.frame(mcl.lipid.cluster[["Cluster"]])
  cluster.l[,'lipid'] <- rownames(correlation)
  colnames(cluster.l) <- c('cluster.l', 'lipid')
  cluster.l[,'raw.order'] <- as.numeric(rownames(cluster.l))
  cluster.l <- arrange(cluster.l, cluster.l)
  cluster.l[,'new.order'] <- order(cluster.l$cluster.l)
  cluster.l <- arrange(cluster.l, raw.order)
  
  cluster.g <- data.frame(mcl.gene.cluster[["Cluster"]])
  cluster.g[,'gene'] <- colnames(correlation)
  colnames(cluster.g) <- c('cluster.g', 'gene')
  cluster.g[,'raw.order'] <- as.numeric(rownames(cluster.g))
  cluster.g <- arrange(cluster.g, cluster.g)
  cluster.g[,'new.order'] <- order(cluster.g$cluster.g)
  cluster.g <- arrange(cluster.g, raw.order)
  
  mat <- correlation[order(cluster.l$new.order),order(cluster.g$cluster.g)]
  
  mat.m <- melt(mat)
  colnames(mat.m) <- c('lipid', 'gene', 'value')
  mat.m <- merge(mat.m, cluster.g, by = 'gene')
  mat.m <- merge(mat.m, cluster.l, by = 'lipid')
}
if(opt$clustering_type != "hierarchical"){
  p <- ggplot(mat.m,aes(gene, lipid))+geom_tile(aes(fill=value)) +
    scale_fill_gradientn(colours = c('#03A9F4','#81D4FA','#FFF9C4','#FFAB91','#FF5722')) +
    theme(text=element_text(size=4),
          axis.text.x=element_text(angle=90,vjust=0), 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          axis.title = element_blank()) +
    facet_grid(cluster.l~cluster.g, scales='free', space='free')
  ggsave(paste0(opt$output, "correlationPlot.pdf"), p)
  
  lipid_group <- split(cluster.l$lipid, cluster.l$cluster.l)
  gene_group <- split(cluster.g$gene, cluster.g$cluster.g)
  if(opt$clustering_type != "k_means"){
    kl <- length(lipid_group)
    kg <- length(gene_group)
  }
  for(i in 1:kl){
    ith <- as.character(i)
    sublipids <- lipid_group[[ith]]
    write.csv(sublipids, paste0(opt$output, "lipids_", ith, ".csv"), row.names = F)
    for(j in 1:kg){
      jth <- as.character(j)
      subgenes <- gene_group[[jth]]
      subdata <- correlation[sublipids, subgenes]
      write.csv(subdata, paste0(opt$output, "correlation_", ith, "_", jth, ".csv"))
      write.csv(subgenes, paste0(opt$output, "genes_", jth, ".csv"), row.names = F)
    }
  }
}

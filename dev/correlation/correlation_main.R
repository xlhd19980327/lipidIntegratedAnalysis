## Do correlation
source("./branch/correlation/readingLipidData_cor2.R")
source("./branch/correlation/readingRNAData_cor.R")

library(pheatmap)
library(tidyverse)
library(MetaboAnalystR)
library(DESeq2)
library(limma)
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
  make_option(c("-o", "--output"), action="store")
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


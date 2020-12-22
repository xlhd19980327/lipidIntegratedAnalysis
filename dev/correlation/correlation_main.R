## Do correlation
source("./branch/correlation/readingLipidData_cor.R")
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
  make_option(c("-f", "--feature_field"), action="store", default = NA, type = "character"), 
  make_option(c("-n", "--NA_string"), action="store", default = NULL, type = "character"), 
  make_option(c("-j", "--input_file_rna"), action="store"),
  make_option(c("-e", "--description_file_rna"), action="store"), 
  make_option(c("-u", "--data_type_rna"), action="store"), 
  make_option(c("-g", "--gene_split"), action="store", default = 4), 
  make_option(c("-l", "--lipid_split"), action="store", default = 4),
  make_option(c("-o", "--output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Data input
lipid_data <- readingLipidData(datafile = opt$input_file, sampleList = opt$description_file,
                               controlGrp = "", dataType = opt$data_type, delOddChainOpt = F,
                               lipField = opt$feature_field, na.char = opt$NA_string)
gene_data <- readingRNAData(datafile = opt$input_file_rna, 
                            sampleList = opt$description_file_rna,
                            type = opt$data_type_rna)

##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

### Waiting time(~2+ min, larger data may require more)
correlation <- cor(lipid_data, t(gene_data), use='pairwise.complete.obs', method ="spearman")
correlation[is.na(correlation)] <- 0
correlation <- correlation[,which(apply(correlation, 2, FUN = sd) != 0)]
correlation <- correlation[which(apply(correlation, 1, FUN = sd) != 0),]
png(paste0(opt$output, "correlationPlot.png"), 
    width = 720, height = 720)
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


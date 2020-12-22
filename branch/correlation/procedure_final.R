## Do correlation
library(pheatmap)
library(tidyverse)
library(MetaboAnalystR)
library(DESeq2)
library(limma)
options(stringsAsFactors = F)
source("./branch/correlation/readingLipidData_cor.R")
source("./branch/correlation/readingRNAData_cor.R")
kg <- 4
kl <- 4

## Data input
lipid_data <- readingLipidData(datafile = "./branch/benchmark/input/HANLipidMediator_imm_forcor.CSV", 
                               controlGrp = "", dataType = "Metabolites", delOddChainOpt = F,
                               lipField = "LipidMediator", sampleList = "./branch/benchmark/input/HANsampleList_lipmid.csv")
gene_data <- readingRNAData(datafile = "./branch/benchmark/input/HANgene_tidy.CSV", 
                              sampleList = "./branch/benchmark/input/HANsampleList.CSV",
                              type = "RNAseq")

##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

### Waiting time(~2+ min, larger data may require more)
correlation <- cor(lipid_data, t(gene_data), use='pairwise.complete.obs', method ="spearman")
correlation[is.na(correlation)] <- 0
correlation <- correlation[,which(apply(correlation, 2, FUN = sd) != 0)]
correlation <- correlation[which(apply(correlation, 1, FUN = sd) != 0),]
### Waiting time(~2+ min, larger data may require more)
list <- pheatmap(correlation,
                 cutree_col = kg, cutree_rows = kl)
png("~/temp/cor/correlationPlot.png", 
    width = 720, height = 720)
list
dev.off()

cutgene <- cutree(list$tree_col, k = kg)[list[["tree_col"]][["order"]]]
cutlipid <- cutree(list$tree_row, k = kl)[list[["tree_row"]][["order"]]]
genegaps <- c(names(cutgene[1]), 
               names(which((cutgene[-1] - cutgene[-length(cutgene)]) != 0)))
lipidgaps <- c(names(cutlipid[1]), 
               names(which((cutlipid[-1] - cutlipid[-length(cutlipid)]) != 0)))
for(i in 1:kl){
  for(j in 1:kg){
    subdata <- correlation[names(cutlipid)[cutlipid == i], names(cutgene)[cutgene == j], drop = F]
    ith <- which(lipidgaps %in% rownames(subdata))
    jth <- which(genegaps %in% colnames(subdata))
    sublipids <- rownames(subdata)
    subgenes <- colnames(subdata)
    write.csv(subdata, paste0("~/temp/cor/", "correlation_", ith, "_", jth, ".csv"))
    write.csv(sublipids, paste0("~/temp/cor/", "lipids_", ith, ".csv"), row.names = F)
    write.csv(subgenes, paste0("~/temp/cor/", "genes_", jth, ".csv"), row.names = F)
  }
}


## Do correlation
library(tidyverse)
library(dplyr)
library(pheatmap)
options(stringsAsFactors = F)
source("./branch/correlation/readingRNAData_cor.R")
source("./branch/correlation/readingLipidData_cor.R")

h_gene=1
h_lipid=1
cutoff_correlation=0.7

## Data input
lipiddataSet <- readingLipidData(datafile = "./branch/benchmark/input/HANLipidMediator_imm_forcor.CSV", 
                                            controlGrp = "", dataType = "Metabolites", delOddChainOpt = F,
                                            lipField = "LipidMediator", sampleList = "./branch/benchmark/input/HANsampleList_lipmid.csv", 
                                            fileLoc = "./branch/benchmark/output/")
genedataSet <- readingRNAData(datafile = "./branch/benchmark/input/HANgene_tidy.CSV", 
                              sampleList = "./branch/benchmark/input/HANsampleList.CSV")

## Tidy the input files
alllipidnames <- lipiddataSet$lipidName
lipid_data <- as.data.frame(lipiddataSet$data)
rownames(lipid_data) <- paste0("L", 1:nrow(lipid_data)) 
lipidnamelink <- data.frame(alllipidnames)
rownames(lipidnamelink) <- paste0("L", 1:nrow(lipidnamelink))
allgenenames <- rownames(genedataSet$data)
gene_data <- as.data.frame(genedataSet$data)
rownames(gene_data) <- paste0("G", 1:nrow(gene_data))
genenamelink <- data.frame(allgenenames)
rownames(genenamelink) <- paste0("G", 1:nrow(genenamelink))

##!!!!!WARNING: This code requires same amount of samples between lipid data and gene data
##!!!!!WARNING: May handle this here later

lipid_data[is.na(lipid_data)] <- 0.5*min(lipid_data, na.rm = T)
### Waiting time(~15+ min)
correlation <- cor(t(lipid_data), t(gene_data), use='pairwise.complete.obs', method ="spearman")
correlation[is.na(correlation)] <- 0

## prepare function
mean_function <- function(x){
  d <- dim(x)
  s <- sum(abs(x))
  s/(d[1]*d[2])
}

merge_list_function <- function(x){
  n <- length(x)
  if(n == 0){
    y <- NULL
  }
  else{
    if(n>1){
      y <- x[[1]]
      for(z in 2:n){
        merge(y,x[[z]])
      }
    }
    else{
      y <- x[[1]]
    }
  }
  return(y)
}
extract_name <- function(listname, direction){
  name <- c()
  if(direction == 'R'){
    for(i in 1:length(listname)){
      name <- c(name, rownames(listname[[i]]))
    }
  }
  if(direction == 'C'){
    for(i in 1:length(listname)){
      name <- c(name, colnames(listname[[i]]))
    }
  }
  return(name)
}
link <- function(namelist){
  n <- dim(namelist)[1]
  namelist[,2]<-NA
  for(i in 1:n){
    namelist[i,2] <- paste0(rownames(namelist)[i], '-',namelist[i,1])
  }
  return(namelist)
}

#extract target lipid and gene ~5min
list <- pheatmap(correlation)
gene_cluster_num <- cutree(list$tree_col, h=h_gene)
gene_cluster <- split(data.frame(t(correlation)), gene_cluster_num)
gene_group_number <- length(gene_cluster)
y=1
result_list <- list()
for (i in 1:gene_group_number){
  gene_number <- dim(gene_cluster[[i]])[1]
  if(gene_number > 1){
    sub_list <- pheatmap(gene_cluster[[i]])
    lipid_cluster_num <- cutree(sub_list$tree_col, h=h_lipid)
    lipid_cluster <- split(data.frame(t(gene_cluster[[i]])), lipid_cluster_num)
    lipid_group_number <- length(lipid_cluster)
    sub_result <- list()
    x=1
    for(j in 1:lipid_group_number){
      sub_mean <- mean_function(lipid_cluster[[j]])
      if(sub_mean > cutoff_correlation){
        sub_result[[x]] <- lipid_cluster[[j]]
        x <- x+1
      }
    }
    result_list[[y]] <- merge_list_function(sub_result)
    y <- y+1
  }
  else{
    df <- gene_cluster[[i]]
    df[abs(df)<cutoff_correlation] <- NA
    df <- t(na.omit(t(df)))
    result_list[[y]] <- t(data.frame(df))
    y <- y+1
  }
  
}

genename<-extract_name(result_list, 'C')
lipidname <- extract_name(result_list, 'R')


### Output
correlation_filt <- correlation[rownames(correlation) %in% unique(lipidname), 
                                colnames(correlation) %in% unique(genename)]
rownames(correlation) <- alllipidnames[as.integer(gsub("^L([0-9]+)", "\\1", 
                                                      rownames(correlation)))]
colnames(correlation) <- allgenenames[as.integer(gsub("^G([0-9]+)", "\\1", 
                                                       colnames(correlation)))]
list_addlipn <- pheatmap(correlation)

rownames(correlation_filt) <- alllipidnames[as.integer(gsub("^L([0-9]+)", "\\1", 
                                                       rownames(correlation_filt)))]
colnames(correlation_filt) <- allgenenames[as.integer(gsub("^G([0-9]+)", "\\1", 
                                                            colnames(correlation_filt)))]
list_filt <- pheatmap(correlation_filt)
# genenamelink <- link(genenamelink)
# lipidnamelink <- link(lipidnamelink)
# 
# genename_test <- as.integer(gsub('G[0-9]+','G'))

### See the names 
up <- c("Maresin 1", "5,6-diHETE ", "5-HEPE ", "5,15-diHETE ", "18-HEPE ", 
        "5-HETE ", "17-HDHA ", "4-HDHA ", "EPA", "12-HETE ", "14-HDHA ", 
        "17R-RvD1 ", "7S,14S-diHDHA ", "RvE2 ", "5,15-diHETE ", "RvE3 ", 
        "13-HDHA ", "Maresin 2", "7-HDHA ", "15-HEPE ", "DHA", "15R-LXA4 ", 
        "AA", "10S,17S-diHDHA", "15-HETE ", "12-HEPE ")
down <- c("AT-LXB4 ", "RvD2 ", "RvD3 ", "17R-RvD1 ", "RvD6 ", "7S,12-trans Mar1")

write.csv(up, "~/temp/up_lipmed.csv")
write.csv(down, "~/temp/down_lipmed.csv")

correlation_filt_up <- correlation_filt[rownames(correlation_filt) %in% up, ]
up_summ <- apply(correlation_filt_up, 2, median)
# up_gene <- rbind(data.frame(up_gene = names(up_summ[up_summ < 0]), median = up_summ[up_summ < 0]), 
#                  data.frame(up_gene = names(up_summ[up_summ > 0]), median = up_summ[up_summ > 0]))
up_regdown_gene <- data.frame(up_gene = names(up_summ[up_summ < 0]), median = up_summ[up_summ < 0])
up_regup_gene <- data.frame(up_gene = names(up_summ[up_summ > 0]), median = up_summ[up_summ > 0])
write.csv(up_regdown_gene, "~/temp/up_regdown_gene.csv")
write.csv(up_regup_gene, "~/temp/up_regup_gene.csv")

correlation_filt_down <- correlation_filt[rownames(correlation_filt) %in% down, ]
down_summ <- apply(correlation_filt_down, 2, median)
down_regdown_gene <- data.frame(down_gene = names(down_summ[down_summ < 0]), median = down_summ[down_summ < 0])
down_regup_gene <- data.frame(down_gene = names(down_summ[down_summ > 0]), median = down_summ[down_summ > 0])
write.csv(down_regdown_gene, "~/temp/down_regdown_gene.csv")
write.csv(down_regup_gene, "~/temp/down_regup_gene.csv")

RvD2 <- c("RvD2 ")
correlation_rvd2 <- correlation[rownames(correlation) %in% RvD2, ]
#data <- rbind(correlation_rvd2, correlation_rvd2)
#pheatmap(data)
correlation_rvd2_regup_gene <- correlation_rvd2[correlation_rvd2 > 0.8]
correlation_rvd2_regdown_gene <- correlation_rvd2[correlation_rvd2 < -0.8]
write.csv(correlation_rvd2_regup_gene, "~/temp/correlation_rvd2_regup_gene.csv")
write.csv(correlation_rvd2_regdown_gene, "~/temp/correlation_rvd2_regdown_gene.csv")




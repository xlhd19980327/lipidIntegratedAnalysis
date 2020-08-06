### For MS-DIAL data ###
### NOTE1: Clients should add group info to the first line of the file. ###
### NOTE1: If using MS Excel, should follow the steps to get csv file: ###
### NOTE1: Click File->Export->Change File Type ->Choose "CSV" options ###
### NOTE1: Or may have some characters garbled, see NOTE1-ref in the code ###
### NOTE2: Clients should put control group should be on the first column ###
### NOTE3: Data should have the following columns: ###
### NOTE3: Metabolite.name ###

### Loading data and environment settings ###
library(tidyverse)
options(stringsAsFactors = F)
setwd("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/LMPDhsa/testData/zsy_DGATinhibitors/ms_dial_DATA/")
data <- read.csv("integ_posnegdata_cos7_addgroupinfo.csv", na.strings = "N/A", skip = 1)
#NOTE1-ref: may have the first character garbled
allgroups <- scan("integ_posnegdata_cos7_addgroupinfo.csv", what = "character", nlines = 1, sep = ",", quote = "\"")
notdataColsLen <- sum(allgroups == '')
allgroups <- allgroups[allgroups != '']
groupsLevel <- unique(allgroups)
nsamples <- table(factor(allgroups, levels = groupsLevel))
getSmapleName <- function(sam, n){
  seqs <- 1:n
  return(paste(sam, seqs, sep = "_"))
}
samplesName <- unlist(mapply(getSmapleName, groupsLevel, nsamples, SIMPLIFY = F))
colnames(data)[(notdataColsLen+1):(notdataColsLen+length(allgroups))] <- samplesName
#NOTE2-ref: should indicate the control group or set on the first column
controlGrp <- groupsLevel[1]

### Clean and Tidy the data ###
## Delete NA (or 0) value
# The lipid is not NA (or 0) at least one sample in a group remains
# ref: 3 samples in a group -- select <= 0.67(2/3) missing value
# ref: 5 samples in a group -- select <= 0.8(4/5) missing value
# Missing value will be replaced by half of the minimum positive value in the original data
percent <- 0.67
data[is.na(data)] <- 0
data_intensity <- data[, -1:-notdataColsLen]
ind_ok <- apply(data_intensity == 0, 1, sum)/ncol(data_intensity) <= percent
data <- data[ind_ok, ]
# Missing value will be replaced by half of the minimum positive value in the original data
minInten <- min(data_intensity[data_intensity > 0], na.rm = T)/2
data[data == 0] <- minInten
## Delete odd FA chain lipids
lipids <- data$Metabolite.name
m <- gregexpr("[0-9]*:", lipids)
fachain <- regmatches(lipids, m)
my.even <-  function(x){
  x <- gsub(":", "", x, perl = T)
  x <- as.numeric(x)
  result <- ifelse(sum(x %% 2) == 0, T, F)
  return(result)
}
ind <- lapply(fachain, my.even)
ind <- unlist(ind)
data <- data[ind, ]
## Only use the lipid class that the library exsists
lipClasses <- c("Cer", "ChE", "DG", "LPC", "LPE", "LPG", "LPI", 
                "MG", "PA", "PC", "PE", "PG", "PI", "PS", 
                "SM", "SPH", "TG", 
                "FA", "LPS", "CL")
#NOTE: class "SPH", "CL", "ChE" will not be considered later due to no occurence in the data
#NOTE: some Glycerophospholipids(eg, "PA", though they donnot occur in the data) will still be considered later due to their similar form to other Glycerophospholipids
Class <- gsub("(.*?) .*", "\\1", data$Metabolite.name, perl = T)
data_sub <- subset(data, subset = Class %in% lipClasses)
## Duplication handle
#Delete some "0:0" info
data_sub_all <- cbind(data_sub, 
                      lipidName = gsub("(.*)\\/0:0$", "\\1", data_sub$Metabolite.name))
data_sub_dup <- data_sub_all %>%
  group_by(lipidName) %>%
  filter(n() > 1)
data_sub_dup_list <- split(data_sub_dup, data_sub_dup$lipidName)
# Rule2: duplicated lipids in different modes- chooese stronger intensity
#        duplicated lipids in the same modes but have different adduct - choose stronger intensity
#        i.e. duplicated lipids all chooese stronger intensity
getSoloInten <- function(x){
  data <- x[, (notdataColsLen+1):(notdataColsLen+length(allgroups))]
  result <- rep(0, length(allgroups))
  for(i in 1:nrow(data)){
    one <- unlist(data[i, ])
    result <- ifelse(one > result, one, result)
  }
  return(result)
}
# End of rule2
data_sub_dup_soloInten <- sapply(data_sub_dup_list, getSoloInten)
# Merge duplication and single ones
data_sub_sing_allhandle <- data_sub_all %>%
  select(-c(1:all_of(notdataColsLen))) %>%
  group_by(lipidName) %>%
  filter(n() == 1) 
data_sub_sing_allhandle <- data_sub_sing_allhandle[, c((length(allgroups)+1), 1:length(allgroups))]
data_sub_dup_allhandle <- rownames_to_column(as.data.frame(t(data_sub_dup_soloInten)))
colnames(data_sub_dup_allhandle) <- colnames(data_sub_sing_allhandle)
data_sub_allhandle <- bind_rows(data_sub_sing_allhandle, data_sub_dup_allhandle)
data_sub_allhandle <- cbind(data_sub_allhandle, 
                            Class = gsub("(.*?) .*", "\\1", data_sub_allhandle$lipidName))
write.csv(data_sub_allhandle, "data_sub_allhandle_cos7.csv", row.names = F)
## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
lipids <- data_sub_allhandle$lipidName
# Get the FA info and MS1 info
getFAsInfo <- function(i){
  Class <- gsub("(.*?) .*", "\\1", i, perl = T)
  if(Class %in% c("Cer", "SM")){ #may use Class %in% c("Cer", "SM", "SPH")
    if(grepl("\\|", i)){
      spbase <- gsub(".*\\|.*?([0-9]+:[0-9]+;[0-9]O).*", "\\1", i, perl = T)
      spbase_OH <- gsub("[0-9]+:[0-9]+;([0-9])O", "\\1", spbase, perl = T)
      spbase_OH_label <- switch (spbase_OH,
                                 '1' = "m",
                                 '2' = "d", 
                                 '3' = "t"
      )
      spbase_carbon <- gsub("([0-9]+:[0-9]+);[0-9]O", "\\1", spbase, perl = T)
      spbase_tidy <- paste0(spbase_OH_label, spbase_carbon)
      fa <- gsub(".*/([0-9]+:[0-9]+).*", "\\1", i, perl = T)
      fas <- c(spbase_tidy, fa)
    } else {
      spall <- gsub(".*?([0-9]+:[0-9]+;[0-9]O).*", "\\1", i, perl = T)
      spall_OH <- gsub("[0-9]+:[0-9]+;([0-9])O", "\\1", spall, perl = T)
      spall_OH_label <- switch (spall_OH,
                                '1' = "m",
                                '2' = "d", 
                                '3' = "t"
      )
      spall_carbon <- gsub("([0-9]+:[0-9]+);[0-9]O", "\\1", spall, perl = T)
      spall_tidy <- paste0(spall_OH_label, spall_carbon)
      fas <- spall_tidy
    }
  } else if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                         "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                         "LPS")){#will not contain "CL" or "ChE"
    # All ignore e/p(O-/P-) connection
    if(grepl("\\|", i)){
      fainfo <- gsub(".*\\|(.*)", "\\1", i)
      fas <- strsplit(fainfo, "[_/]")
      m <- gregexpr("[0-9]+:[0-9]+", fas)
      fas <- unlist(regmatches(as.character(fas), m))
    } else {
      m <- gregexpr("[0-9]+:[0-9]+", i)
      fas <- unlist(regmatches(i, m))
    }
  }
  n <- length(fas)
  if(n > 1){
    ms1 <- F
  } else{
    if(Class %in% c("Cer", "DG", "PA", "PC", "PE", "PG", "PI", "PS", "SM", "TG")){
      ms1 <- T
    } else {
      ms1 <- F
    }
  }
  result <- list(fas, ms1)
  names(result) <- c(Class, "ms1")
  return(result)
}
fasInfo <- lapply(lipids, getFAsInfo)
ms1Info <- sapply(fasInfo, function(x) x$ms1)
data_sub_allhandle <- cbind(data_sub_allhandle, ms1 = ms1Info)
data_sub_ms1 <- subset(data_sub_allhandle, subset = ms1 == T)
data_sub_ms2 <- subset(data_sub_allhandle, subset = ms1 == F)
# Tidy and integrate itensity of lipid class containing FA chain info
# Use MS2 to do the later statistics only
lipid_subclass_handle <- data.frame()
for(i in 1:nrow(data_sub_allhandle)){
  # e,p connection info will be ignored
  fa <- fasInfo[[i]][[1]]
  Class <- names(fasInfo[[i]][1])
  subclass <- paste0(Class, "(", fa, ")")
  lipid_subclass_handle <- rbind(lipid_subclass_handle, 
                                 cbind(subclass = subclass, 
                                       bind_rows(replicate(length(subclass), data_sub_allhandle[i, ], simplify = FALSE))))
}
lipid_subclass_handle <- subset(lipid_subclass_handle, subset = ms1 == F)
write.csv(lipid_subclass_handle, "lipid_subclass_handle_cos7.csv", row.names = F)

### Tidy for and do Visualization ###
## Use data_sub_allhandle(Delete duplication) to calculate itensity of each lipid class
data_sub_classSum_stat <- data_sub_allhandle %>%
  ungroup() %>%
  select(-ms1, -lipidName) %>%
  gather(key = "case", value = "lipidsum", -Class) %>%
  group_by(Class, case) %>%
  summarise(lipidsum = sum(lipidsum)) %>%
  mutate(group = gsub("(.*?)_[0-9]+$", "\\1", case)) 
data_sub_classSum_stat2 <- data_sub_classSum_stat %>%
  group_by(Class, group) %>%
  summarise(mean = mean(lipidsum) / n(), 
            realmean = mean(lipidsum),
            sd = sd(lipidsum))
data_sub_classSum_p <- split(data_sub_classSum_stat, data_sub_classSum_stat$Class)
getPValue <- function(x){
  Class <- unique(x$Class)
  #Totally equal to t.test() one by one
  p <- pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value
  p_tidy <- gather(rownames_to_column(as.data.frame(p), var = "group"), 
                   -group, key = "cprGroup", value = "p")
  p_tidy <- subset(p_tidy, cprGroup == controlGrp, select = -cprGroup) ## the control group
  result <- cbind(p_tidy, Class = Class)
  return(result)
}
data_sub_classSum_stat3 <- lapply(data_sub_classSum_p, getPValue)
data_sub_classSum_stat3 <- do.call(rbind, data_sub_classSum_stat3)
data_sub_classSum_integStat <- left_join(data_sub_classSum_stat, data_sub_classSum_stat2)
data_sub_classSum_integStat <- left_join(data_sub_classSum_integStat, data_sub_classSum_stat3)
addSigLabel <- function(x){
  x[is.na(x)] <- 1
  label <- c('', "*", "**", "***")
  allLable <- ifelse(x >= 0.05, label[1], NA)
  allLable <- ifelse(x < 0.05 & x >= 0.01, label[2], allLable)
  allLable <- ifelse(x < 0.01 & x >= 0.001, label[3], allLable)
  allLable <- ifelse(x < 0.001, label[4], allLable)
}
sigLabel <- addSigLabel(data_sub_classSum_integStat$p)
data_sub_classSum_integStat <- cbind(data_sub_classSum_integStat, sigLabel = sigLabel)
ggplot(data = data_sub_classSum_integStat, aes(x = group)) +
  geom_bar(aes(y = mean, fill = group), stat = "identity") +
  geom_errorbar(aes(ymin = realmean - sd, ymax = realmean + sd, color = group), width = 0.2) +
  geom_dotplot(aes(y = lipidsum), binaxis='y', stackdir='center') +
  geom_text(aes(y = realmean+1.5*sd, label = sigLabel),
            size = 3, fontface = "bold", color = "red") +
  theme_classic() +
  facet_wrap(~Class, scales="free") 
## Use lipid_subclass_handle to calculate itensity of lipid class containing FA chain info(subclass)
dupCalcHandle <- function(x){
  divisor <- c(1, 2, 3)
  #contain every lipids in the db except "CL"
  allDivisor <- ifelse(x %in% c("ChE", "LPC", "LPE", "LPG", "LPI", "MG", "SPH", "FA", "LPA", "LPS"), divisor[1], NA)
  allDivisor <- ifelse(x %in% c("DG", "PA", "PC", "PE", "PG", "PI", "PS", "Cer", "SM"), divisor[2], allDivisor)
  allDivisor <- ifelse(x %in% c("TG"), divisor[3], allDivisor)
  return(allDivisor)
}
lipid_subclass_handle <- cbind(lipid_subclass_handle, 
                               divisor = dupCalcHandle(lipid_subclass_handle$Class))
lipid_subclass_stat <- lipid_subclass_handle %>%
  ungroup() %>%
  select(-ms1, -lipidName) %>%
  gather(key = "case", value = "lipidsum", -subclass, -Class, -divisor) %>%
  group_by(subclass, case, Class) %>%
  summarise(lipidsum = sum(lipidsum / divisor)) %>%
  mutate(group = gsub("(.*?)_[0-9]+$", "\\1", case))
lipid_subclass_stat2 <- lipid_subclass_stat %>%
  group_by(subclass, group) %>%
  summarise(mean = mean(lipidsum) / n(), 
            realmean = mean(lipidsum) ,
            sd = sd(lipidsum))
lipid_subclass_stat_p <- split(lipid_subclass_stat, lipid_subclass_stat$subclass)
getPValue <- function(x){
  subclass <- unique(x$subclass)
  #Totally equal to t.test() one by one
  p <- pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value
  p_tidy <- gather(rownames_to_column(as.data.frame(p), var = "group"), 
                   -group, key = "cprGroup", value = "p")
  p_tidy <- subset(p_tidy, cprGroup == controlGrp, select = -cprGroup) ## the control group
  result <- cbind(p_tidy, subclass = subclass)
  return(result)
}
lipid_subclass_stat3 <- lapply(lipid_subclass_stat_p, getPValue)
lipid_subclass_stat3 <- do.call(rbind, lipid_subclass_stat3)
lipid_subclass_integStat <- left_join(lipid_subclass_stat, lipid_subclass_stat2)
lipid_subclass_integStat <- left_join(lipid_subclass_integStat, lipid_subclass_stat3)
sigLabel2 <- addSigLabel(lipid_subclass_integStat$p)
lipid_subclass_integStat <- cbind(lipid_subclass_integStat, sigLabel = sigLabel2)
for(i in lipClasses){
  oneLipClassData <- subset(lipid_subclass_integStat, 
                            subset = Class == i)
  if(nrow(oneLipClassData) != 0){
    ggplot(data = oneLipClassData, aes(x = group)) +
      geom_bar(aes(y = mean, fill = group), stat = "identity") +
      geom_errorbar(aes(ymin = realmean - sd, ymax = realmean + sd, color = group), width = 0.2) +
      geom_dotplot(aes(y = lipidsum), binaxis='y', stackdir='center') +
      geom_text(aes(y = realmean+1.5*sd, label = sigLabel),
                size = 3, fontface = "bold", color = "red") +
      theme_classic() +
      facet_wrap(~subclass, scales="free") 
    ggsave(paste0("./lipsubclassFigures_facet_cos7/", i, ".pdf"), dpi = 300, width = 15, height = 10)
  }
}

### Network construction && Statistics use our methods###
## Data preparation
# Tidy the lipid_subclass_integStat and retain the useful info
lipid_subclass_tidyStat <- lipid_subclass_integStat %>%
  ungroup() %>%
  select(subclass, Class, group, realmean) %>%
  distinct() %>%
  group_by(subclass) %>%
  mutate(ind = ifelse(group == controlGrp, "A", "B")) %>% 
  arrange(ind, .by_group = T) %>% ## put controlGrp into ind[1] so that we can recognize the control group
  mutate(log2FC = log2(realmean / realmean[1])) %>% 
  filter(group != controlGrp) %>%
  select(-ind)
write.csv(lipid_subclass_tidyStat, "lipid_subclass_tidyStat_cos7.csv", row.names = F)
# Lipid reaction library file input
lipReact <- read.csv("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/LMPDhsa/lip_gene_intergReac/hsa_lipidreact.csv")
allReact <- read.csv("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/LMPDhsa/hsa_all_integData.csv")  
integReact <- inner_join(lipReact, allReact, by = "reaction_id")
integReact <- subset(integReact, 
                     select = c("reaction_id", "gene_symbol", "R_lipid.x", "P_lipid.x"))
colnames(integReact) <- c("reaction_id", "gene_symbol", "R_lipid", "P_lipid")
#use unique/distinct
integReact <- unique(integReact)
# After add libnames to the lipid_subclass_stat manually...
# ...
# !!!!!WARNING: Some libnames will be wrong, cause carbon base info add
lipidsubclass <- read.csv("lipid_subclass_tidyStat_cos7_addlibnames.csv")
lipsubname <- levels(factor(lipidsubclass$libnames))
lipsubname <- unlist(strsplit(lipsubname, ";"))
lipsubname <- levels(factor(lipsubname))
# Only use the reactions that lipid Class both corresponding to Reactant and Product
ind <- integReact$R_lipid %in% lipsubname & integReact$P_lipid %in% lipsubname
integReact_filter <- integReact[ind, ]
# Calculate library lipid names ind from data 
cal_libInd <- function(x, y){
  libname <- unlist(strsplit(x, ";"))
  getInd <- function(i, lib = y){
    return(which(lib == i))
  }
  ind <- lapply(libname, getInd)
  ind <- unlist(ind)
  return(ind)
}
R_ind <- lapply(lipidsubclass$libnames, cal_libInd, integReact_filter$R_lipid)
P_ind <- lapply(lipidsubclass$libnames, cal_libInd, integReact_filter$P_lipid)
# Transfer and integrate library lipid names ind to data ind
integ_ind <- function(libInd, type = c("R_dataInd", "P_dataInd")){
  n <- length(libInd)
  result <- data.frame()
  for(i in 1:n){
    dataIndi <- i
    libIndi <- libInd[[i]]
    result <- rbind(result, 
                    data.frame(libInd = libIndi, rep(dataIndi, length(libIndi))))
  }
  colnames(result)[2] <- type
  return(result)
}
R_ind <- integ_ind(R_ind, "R_dataInd")
P_ind <- integ_ind(P_ind, "P_dataInd")
reaction_ind <- full_join(R_ind, P_ind)
reaction_ind_list <- split(reaction_ind, reaction_ind$libInd)
getReactionInfo <- function(ind, data = lipidsubclass){
  libind <- unique(ind$libind)
  getInd <- function(x, y = ind){
    res <- y[[x]][!is.na(y[[x]])]
    res <- unique(res)
    return(res)
  }
  R_lipid_ind <- getInd("R_dataInd")
  R_lipid <- data[R_lipid_ind, ]
  P_lipid_ind <- getInd("P_dataInd")
  P_lipid <- data[P_lipid_ind, ]
  return(list(
    R_lipid = R_lipid, 
    P_lipid = P_lipid
  ))
}
reaction_info <- lapply(reaction_ind_list, getReactionInfo)
names(reaction_info) <- integReact_filter[as.numeric(names(reaction_info)), ]$gene_symbol
# Use two sample comparation once
lipsample <- "OADGAT1i"
## Network construction(for cytoscape visualization)
network <- data.frame()
for(i in 1:length(reaction_info)){
  data <- reaction_info[[i]]
  gene <- names(reaction_info)[i]
  R_info <- subset(data[["R_lipid"]], subset = group == lipsample)$subclass
  P_info <- subset(data[["P_lipid"]], subset = group == lipsample)$subclass
  n <- max(length(R_info), length(P_info))
  length(R_info) <- n
  length(P_info) <- n
  integ_info <- data.frame(R = R_info, P = P_info)
  network_i <- apply(integ_info, 1, 
                     function(x) data.frame(from = c(x[1], gene), to = c(gene, x[2])))
  network_i <- do.call(rbind, network_i)
  network_i <- network_i[complete.cases(network_i), ]
  network <- rbind(network, network_i)
}
nodeAttr <- rbind(
  subset(lipidsubclass, select = c(subclass, log2FC, Class)), 
  data.frame(subclass = names(reaction_info), log2FC = NA, Class = "gene")
)
colnames(nodeAttr) <- c("nodes", "log2FC", "Class")
write.csv(network, paste0("network_", lipsample, ".csv"), row.names = F)
write.csv(nodeAttr, paste0("nodeAttr_", lipsample, ".csv"), row.names = F)
## Statistics use our methods
getRegState <- function(data){
  R_info <- subset(data[["R_lipid"]], subset = group == lipsample)
  P_info <- subset(data[["P_lipid"]], subset = group == lipsample)
  getInfo <- function(info, type = c("R", "P")){
    n <- length(info$log2FC)
    up <- sum(info$log2FC > log2(1.5))
    down <- sum(info$log2FC < -log2(1.5))
    zero <- n - up - down 
    log2FC_avg <- sum(info$log2FC) / n
    up_pct <- up / n
    down_pct <- down / n
    if(up > down){
      result <- "up"
    } else if(up == down){
      result <- "zero"
    } else{
      result <- "down"
    }
    data_result <- data.frame(up, down, zero, log2FC_avg, n, result, up_pct, down_pct)
    colnames(data_result) <- paste0(type, "_", colnames(data_result))
    return(data_result)
  }
  R_result <- getInfo(R_info, type = "R")
  P_result <- getInfo(P_info, type = "P")
  P_R <- P_result$P_log2FC_avg - R_result$R_log2FC_avg
  # Calculate p-value
  # Rmarginal <- ifelse(P_R > 0, R_result$R_down, R_result$R_up)
  # Pmarginal <- ifelse(P_R > 0, P_result$P_down, P_result$P_up)
  Rmarginal <- ifelse(R_result$R_log2FC_avg > 0, R_result$R_down, R_result$R_up)
  Pmarginal <- ifelse(P_result$P_log2FC_avg > 0, P_result$P_down, P_result$P_up)
  PallValue <- Pmarginal:P_result$P_n
  P_pvalue <- min(sum(choose(P_result$P_n, PallValue)*(1/3)^PallValue), 1)
  RallValue <- Rmarginal:R_result$R_n
  R_pvalue <- min(sum(choose(R_result$R_n, RallValue)*(1/3)^RallValue), 1)
  p_value_strict <- 1 - (1-P_pvalue) * (1-R_pvalue)
  regState <- cbind(R_result, P_result, P_R, Rmarginal, Pmarginal, p_value = p_value_strict)
  return(regState)
}
regState <- lapply(reaction_info, getRegState)
regState <- cbind(gene = names(regState), 
                  do.call(rbind, regState))
write.csv(regState, "regState_2.csv", row.names = F)

### For MetaboAnalyst analysis ###
#Results will retain the ':' or '_' character(i.e. the original characters)
library(MetaboAnalystR)
mSet<-InitDataObjects("conc", "stat", FALSE)
##!!!NOTE: Replace "Replacing_with_your_file_path" to the file name
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "colu", "disc");
mSet<-SanityCheckData(mSet)
#Retain the ':' or '_' character(i.e. the original characters)
mynames <- mSet[["dataSet"]][["orig.var.nms"]]
mSet[["dataSet"]][["cmpd"]] <- mynames
names(mSet[["dataSet"]][["orig.var.nms"]]) <- mynames
#colnames(mSet[["dataSet"]][["orig"]]) <- mynames
colnames(mSet[["dataSet"]][["preproc"]]) <- mynames
mSet<-ReplaceMin(mSet)
#Use Interquantile range (IQR) for Data Filtering
mSet<-FilterVariable(mSet, "iqr", "F", 25)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
#!!!NOTE: volcano plot should concern which the control group is
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 1, 0.75,F, 0.1, TRUE, "raw")
mSet<-PlotVolcano(mSet, "volcano_2_",1, "pdf", 72, width=NA)
#t.test, for subheatmap plotting use
mSet <- Ttests.Anal(mSet)
mSet<-PCA.Anal(mSet)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotHeatMap(mSet, "heatmap_3_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F)
mSet<-PlotSubHeatMap(mSet, "heatmap_2_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 75, "detail", F, T, T, F)

## Temp procedure ##
## Customize heatmap plotting(names will not be restricted to 18 length)
my.plot <- function (mSetObj = NA, imgName, format = "png", dpi = 72, width = NA, 
                     dataOpt, scaleOpt, smplDist, clstDist, palette, viewOpt = "detail", 
                     rowV = T, colV = T, var.inx = NA, border = T, grp.ave = F) 
{
  filenm = paste0(imgName, ".json")
  mSetObj$analSet$htmap <- list(dist.par = smplDist, clust.par = clstDist)
  if (dataOpt == "norm") {
    my.data <- mSetObj$dataSet$norm
  }
  else {
    my.data <- mSetObj$dataSet$prenorm
  }
  if (is.na(var.inx)) {
    hc.dat <- as.matrix(my.data)
  }
  else {
    hc.dat <- as.matrix(my.data[, var.inx])
  }
  #colnames(hc.dat) <- substr(colnames(hc.dat), 1, 18)
  if (mSetObj$dataSet$type.cls.lbl == "integer") {
    hc.cls <- as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls])
  }
  else {
    hc.cls <- mSetObj$dataSet$cls
  }
  if (grp.ave) {
    lvs <- levels(hc.cls)
    my.mns <- matrix(ncol = ncol(hc.dat), nrow = length(lvs))
    for (i in 1:length(lvs)) {
      inx <- hc.cls == lvs[i]
      my.mns[i, ] <- apply(hc.dat[inx, ], 2, mean)
    }
    rownames(my.mns) <- lvs
    colnames(my.mns) <- colnames(hc.dat)
    hc.dat <- my.mns
    hc.cls <- as.factor(lvs)
  }
  if (palette == "gbr") {
    colors <- colorRampPalette(c("green", "black", "red"), 
                               space = "rgb")(256)
  }
  else if (palette == "heat") {
    colors <- heat.colors(256)
  }
  else if (palette == "topo") {
    colors <- topo.colors(256)
  }
  else if (palette == "gray") {
    colors <- colorRampPalette(c("grey90", "grey10"), space = "rgb")(256)
  }
  else {
    colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, 
                                                            "RdBu"))(256))
  }
  imgName = paste(imgName, "dpi", dpi, ".", format, sep = "")
  if (is.na(width)) {
    minW <- 630
    myW <- nrow(hc.dat) * 18 + 150
    if (myW < minW) {
      myW <- minW
    }
    w <- round(myW/72, 2)
  }
  else if (width == 0) {
    w <- 7.2
  }
  else {
    w <- 7.2
  }
  mSetObj$imgSet$heatmap <- imgName
  myH <- ncol(hc.dat) * 18 + 150
  h <- round(myH/72, 2)
  if (viewOpt == "overview") {
    if (is.na(width)) {
      if (w > 9) {
        w <- 9
      }
    }
    else if (width == 0) {
      if (w > 7.2) {
        w <- 7.2
      }
    }
    else {
      w <- 7.2
    }
    if (h > w) {
      h <- w
    }
    mSetObj$imgSet$heatmap <- imgName
  }
  if (grp.ave) {
    w <- nrow(hc.dat) * 25 + 300
    w <- round(w/72, 2)
  }
  if (border) {
    border.col <- "grey60"
  }
  else {
    border.col <- NA
  }
  if (format == "pdf") {
    pdf(file = imgName, width = w, height = h, bg = "white", 
        onefile = FALSE)
  }
  else {
    Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, 
                 width = w, height = h, type = format, bg = "white")
  }
  if (mSetObj$dataSet$cls.type == "disc") {
    annotation <- data.frame(class = hc.cls)
    rownames(annotation) <- rownames(hc.dat)
    if (palette == "gray") {
      cols <- MetaboAnalystR:::GetColorSchema(mSetObj, T)
      uniq.cols <- unique(cols)
    }
    else {
      cols <- MetaboAnalystR:::GetColorSchema(mSetObj)
      uniq.cols <- unique(cols)
    }
    if (mSetObj$dataSet$type.cls.lbl == "integer") {
      cls <- as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls])
    }
    else {
      cls <- mSetObj$dataSet$cls
    }
    names(uniq.cols) <- unique(as.character(sort(cls)))
    ann_colors <- list(class = uniq.cols)
    pheatmap::pheatmap(t(hc.dat), annotation = annotation, 
                       fontsize = 8, fontsize_row = 8, clustering_distance_rows = smplDist, 
                       clustering_distance_cols = smplDist, clustering_method = clstDist, 
                       border_color = border.col, cluster_rows = colV, 
                       cluster_cols = rowV, scale = scaleOpt, color = colors, 
                       annotation_colors = ann_colors)
    dat = t(hc.dat)
    if (scaleOpt == "row") {
      res <- t(apply(dat, 1, function(x) {
        as.numeric(cut(x, breaks = 30))
      }))
    }
    else {
      res <- t(apply(dat, 2, function(x) {
        as.numeric(cut(x, breaks = 30))
      }))
    }
    colnames(dat) = NULL
    netData <- list(data = res, annotation = annotation, 
                    smp.nms = colnames(t(hc.dat)), met.nms = rownames(t(hc.dat)), 
                    colors = colors)
    sink(filenm)
    cat(RJSONIO::toJSON(netData))
    sink()
  }
  else {
    heatmap(hc.dat, Rowv = rowTree, Colv = colTree, col = colors, 
            scale = "column")
    dat = t(hc.dat)
    if (scaleOpt == "row") {
      res <- t(apply(dat, 1, function(x) {
        as.numeric(cut(x, breaks = 30))
      }))
    }
    else {
      res <- t(apply(dat, 2, function(x) {
        as.numeric(cut(x, breaks = 30))
      }))
    }
    colnames(dat) = NULL
    netData <- list(data = res, annotation = "NA", smp.nms = colnames(t(hc.dat)), 
                    met.nms = rownames(t(hc.dat)), colors = colors)
    sink(filenm)
    cat(RJSONIO::toJSON(netData))
    sink()
  }
  dev.off()
  return(mSetObj)
}
my.subplot <- function (mSetObj = NA, imgName, format = "png", dpi = 72, width = NA, 
                        dataOpt, scaleOpt, smplDist, clstDist, palette, method.nm, 
                        top.num, viewOpt, rowV = T, colV = T, border = T, grp.ave = F) 
{
  var.nms = colnames(mSetObj$dataSet$norm)
  if (top.num < length(var.nms)) {
    if (method.nm == "tanova") {
      if (MetaboAnalystR:::GetGroupNumber(mSetObj) == 2) {
        if (is.null(mSetObj$analSet$tt)) {
          Ttests.Anal(mSetObj)
        }
        var.nms <- names(sort(mSetObj$analSet$tt$p.value))[1:top.num]
      }
      else {
        if (is.null(mSetObj$analSet$aov)) {
          ANOVA.Anal(mSetObj)
        }
        var.nms <- names(sort(mSetObj$analSet$aov$p.value))[1:top.num]
      }
    }
    else if (method.nm == "cor") {
      if (is.null(mSetObj$analSet$cor.res)) {
        Match.Pattern(mSetObj)
      }
      cor.res <- mSetObj$analSet$cor.res
      ord.inx <- order(cor.res[, 3])
      cor.res <- cor.res[ord.inx, ]
      ord.inx <- order(cor.res[, 1])
      cor.res <- cor.res[ord.inx, ]
      var.nms <- rownames(cor.res)[1:top.num]
    }
    else if (method.nm == "vip") {
      if (is.null(mSetObj$analSet$plsda)) {
        PLSR.Anal(mSetObj)
        PLSDA.CV(mSetObj)
      }
      vip.vars <- mSetObj$analSet$plsda$vip.mat[, 1]
      var.nms <- names(rev(sort(vip.vars)))[1:top.num]
    }
    else if (method.nm == "rf") {
      if (is.null(analSet$rf)) {
        RF.Anal(mSetObj)
      }
      var.nms <- GetRFSigRowNames()[1:top.num]
    }
  }
  var.inx <- match(var.nms, colnames(mSetObj$dataSet$norm))
  my.plot(mSetObj, imgName, format, dpi, width, dataOpt, 
          scaleOpt, smplDist, clstDist, palette, viewOpt, rowV, 
          colV, var.inx, border, grp.ave)
}

mSet<-my.subplot(mSet, "heatmap_2_", "pdf", 300, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 75, "detail", F, T, T, F)


### Lipidr procedure
##!! Not complete yet!
options(stringsAsFactors = F)
library(lipidr)
data <- read.csv("data_sub_allhandle_cos7.csv")
data <- subset(data, select = -Class)
lipid <- data$lipidName
## other modify
getTidyInfo <- function(i){
  Class <- gsub("(.*?) .*", "\\1", i, perl = T)
  if(Class %in% c("Cer", "SM")){ #may use Class %in% c("Cer", "SM", "SPH")
    if(grepl("\\|", i)){
      spbase <- gsub(".*\\|.*?([0-9]+:[0-9]+;[0-9]O).*", "\\1", i, perl = T)
      spbase_OH <- gsub("[0-9]+:[0-9]+;([0-9])O", "\\1", spbase, perl = T)
      spbase_OH_label <- switch (spbase_OH,
                                 '1' = "m",
                                 '2' = "d", 
                                 '3' = "t"
      )
      spbase_carbon <- gsub("([0-9]+:[0-9]+);[0-9]O", "\\1", spbase, perl = T)
      spbase_tidy <- paste0(spbase_OH_label, spbase_carbon)
      fa <- gsub(".*/([0-9]+:[0-9]+).*", "\\1", i, perl = T)
      fas <- c(spbase_tidy, fa)
    } else {
      spall <- gsub(".*?([0-9]+:[0-9]+;[0-9]O).*", "\\1", i, perl = T)
      spall_OH <- gsub("[0-9]+:[0-9]+;([0-9])O", "\\1", spall, perl = T)
      spall_OH_label <- switch (spall_OH,
                                '1' = "m",
                                '2' = "d", 
                                '3' = "t"
      )
      spall_carbon <- gsub("([0-9]+:[0-9]+);[0-9]O", "\\1", spall, perl = T)
      spall_tidy <- paste0(spall_OH_label, spall_carbon)
      fas <- spall_tidy
    }
  } else if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                         "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                         "LPS")){#will not contain "CL" or "ChE"
    # Should handle e/p(O-/P-) connection, e/p will write in the front of the FA chain info
    i <- gsub("O\\-", "e", i)
    i <- gsub("P\\-", "p", i)
    if(grepl("\\|", i)){
      fainfo <- gsub(".*\\|(.*)", "\\1", i)
      fas <- unlist(strsplit(fainfo, "[_/]"))
      m <- gregexpr("(e|p)?[0-9]+:[0-9]+", fas)
      fas <- unlist(regmatches(as.character(fas), m))
    } else {
      m <- gregexpr("(e|p)?[0-9]+:[0-9]+", i)
      fas <- unlist(regmatches(i, m))
    }
  }
  result <- paste0(Class, " ", paste0(fas, collapse = "/"))
  return(result)
}
lipid <- sapply(lipid, getTidyInfo)
data$lipidName <- lipid
## Annotation files
Sample <- colnames(data)[-1]
SampleGroup <- gsub("(.*)_[0-9]+", "\\1", Sample)
SampleNum <- gsub(".*_([0-9]+)", "\\1", Sample)
anno <- data.frame(Sample, SampleGroup, SampleNum)
d <- as_lipidomics_experiment(data)
d <- add_sample_annotation(d, anno)
plot_samples(d, type = "tic", log = TRUE)
plot_molecules(d, "sd", measure = "Area", log = FALSE)




### For LipidSearch data ###
### NOTE: First merge pos+neg data(if they are separate)
### NOTE1: Clients should add group info to the first line of the file. ###
### NOTE1: If using MS Excel, should follow the steps to get csv file: ###
### NOTE1: Click File->Export->Change File Type ->Choose "CSV" options ###
### NOTE1: Or may have some characters garbled, see NOTE1-ref in the code ###
### NOTE2: Clients should put control group should be on the first column ###
### NOTE3: Data should have the following columns: ###
### NOTE3: Class, FattyAcid ###

### Loading data and environment settings ###
library(tidyverse)
options(stringsAsFactors = F)
setwd("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/procedure/github/testData/zsy_DGATinhibitors/cos7Data/input")
data <- read.csv("Cos7_integ.csv", skip = 1)
#NOTE1-ref: may have the first character garbled
allgroups <- scan("Cos7_integ.csv", what = "character", nlines = 1, sep = ",", quote = "\"")
notdataColsLen <- sum(allgroups == '')
allgroups <- allgroups[allgroups != '']
groupsLevel <- unique(allgroups)
nsamples <- table(factor(allgroups, levels = groupsLevel))
names(allgroups) <- colnames(data)[(notdataColsLen+1):(notdataColsLen+length(allgroups))]
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
fas <- data$FattyAcid
m <- gregexpr("[0-9]*:", fas)
fachain <- regmatches(fas, m)
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
#NOTE: class "CL" will not be considered later due to no occurence in the data
#NOTE: some Glycerophospholipids(eg, "PA", though they donnot occur in the data) will still be considered later due to their similar form to other Glycerophospholipids
Class <- data$Class
data_sub <- subset(data, subset = Class %in% lipClasses)
## Duplication handle
data_sub_all <- cbind(data_sub, 
                      lipidName = paste0(data_sub$Class, data_sub$FattyAcid))
data_sub_dup <- data_sub_all %>%
  group_by(lipidName) %>%
  filter(n() > 1)
data_sub_dup_list <- split(data_sub_dup, data_sub_dup$lipidName)
# Rule2: duplicated lipids in different modes- chooese stronger intensity
#        duplicated lipids in the same modes but have different adduct - choose stronger intensity
#        i.e. duplicated lipids all choose stronger intensity
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
                            Class = gsub("(.*?)\\(.*", "\\1", data_sub_allhandle$lipidName))
write.csv(data_sub_allhandle, "data_sub_allhandle_cos7.csv", row.names = F)
## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
lipids <- data_sub_allhandle$lipidName
# Get the FA info and MS1 info
getFAsInfo <- function(i){
  Class <- gsub("(.*?)\\(.*", "\\1", i, perl = T)
  if(Class %in% c("Cer", "SM", "SPH")){ 
    # Ignore "+O" info
    if(grepl("_", i)){
      spbase <- gsub(".*\\(([mdt][0-9]+:[0-9]+)_.*", "\\1", i, perl = T)
      fa <- gsub(".*_([0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- c(spbase, fa)
    } else {
      spall <- gsub(".*\\(([mdt][0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- spall
    }
  } else if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                         "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                         "LPS", "ChE")){#will not contain "CL"
    # All ignore e/p(O-/P-) connection
    if(grepl("_", i)){
      fainfo <- gsub(".*\\((.*)\\)", "\\1", i)
      fas <- strsplit(fainfo, "[_/]")
      m <- gregexpr("[0-9]+:[0-9]+", fas)
      fas <- unlist(regmatches(as.character(fas), m))
    } else {
      m <- gregexpr("[0-9]+:[0-9]+", i)
      fas <- unlist(regmatches(i, m))
      ## ChE may have some 0:0 info 
      if(length(fas) == 0){
        fas <- "0:0"
      }
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
    ggsave(paste0("./lipsubclassFigures_facet/", i, ".pdf"), dpi = 300, width = 15, height = 10)
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
lipReact <- read.csv("hsa_lipidreact.csv")
allReact <- read.csv("hsa_all_integData.csv")  
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
lipsample <- "COSOADGAT1i"
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
# ## Statistics use our methods
# getRegState <- function(data){
#   R_info <- subset(data[["R_lipid"]], subset = group == lipsample)
#   P_info <- subset(data[["P_lipid"]], subset = group == lipsample)
#   getInfo <- function(info, type = c("R", "P")){
#     n <- length(info$log2FC)
#     up <- sum(info$log2FC > log2(1.5))
#     down <- sum(info$log2FC < -log2(1.5))
#     zero <- n - up - down 
#     log2FC_avg <- sum(info$log2FC) / n
#     up_pct <- up / n
#     down_pct <- down / n
#     if(up > down){
#       result <- "up"
#     } else if(up == down){
#       result <- "zero"
#     } else{
#       result <- "down"
#     }
#     data_result <- data.frame(up, down, zero, log2FC_avg, n, result, up_pct, down_pct)
#     colnames(data_result) <- paste0(type, "_", colnames(data_result))
#     return(data_result)
#   }
#   R_result <- getInfo(R_info, type = "R")
#   P_result <- getInfo(P_info, type = "P")
#   P_R <- P_result$P_log2FC_avg - R_result$R_log2FC_avg
#   # Calculate p-value
#   # Rmarginal <- ifelse(P_R > 0, R_result$R_down, R_result$R_up)
#   # Pmarginal <- ifelse(P_R > 0, P_result$P_down, P_result$P_up)
#   Rmarginal <- ifelse(R_result$R_log2FC_avg > 0, R_result$R_down, R_result$R_up)
#   Pmarginal <- ifelse(P_result$P_log2FC_avg > 0, P_result$P_down, P_result$P_up)
#   PallValue <- Pmarginal:P_result$P_n
#   P_pvalue <- min(sum(choose(P_result$P_n, PallValue)*(1/3)^PallValue), 1)
#   RallValue <- Rmarginal:R_result$R_n
#   R_pvalue <- min(sum(choose(R_result$R_n, RallValue)*(1/3)^RallValue), 1)
#   p_value_strict <- 1 - (1-P_pvalue) * (1-R_pvalue)
#   regState <- cbind(R_result, P_result, P_R, Rmarginal, Pmarginal, p_value = p_value_strict)
#   return(regState)
# }
# regState <- lapply(reaction_info, getRegState)
# regState <- cbind(gene = names(regState), 
#                   do.call(rbind, regState))
# write.csv(regState, "regState_2.csv", row.names = F)

## Statistics use our methods, use t-test to test the statistics(log2FC difference)
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
    data_result <- data.frame(log2FC_avg, n, up_pct, down_pct)
    colnames(data_result) <- paste0(type, "_", colnames(data_result))
    return(data_result)
  }
  R_result <- getInfo(R_info, type = "R")
  P_result <- getInfo(P_info, type = "P")
  P_R <- P_result$P_log2FC_avg - R_result$R_log2FC_avg
  p_value <- ifelse(P_result$P_n >= 3 & R_result$R_n >= 3, 
                    t.test(R_info$log2FC, P_info$log2FC)$p.value, 
                    NA)
  regState <- cbind(R_result, P_result, P_R, p_value = p_value)
  return(regState)
}
regState <- lapply(reaction_info, getRegState)
regState <- cbind(gene = names(regState), 
                  do.call(rbind, regState))
write.csv(regState, "regState_3.csv", row.names = F)

### Lipidr procedure
options(stringsAsFactors = F)
library(lipidr)
use_interactive_graphics(interactive=T)
data <- read.csv("data_sub_allhandle_cos7.csv")
data <- subset(data, select = -Class)
lipid <- data$lipidName
## Modify the format
getTidyInfo <- function(i){
  Class <- gsub("(.*?)\\(.*", "\\1", i, perl = T)
  if(Class %in% c("Cer", "SM", "SPH")){ 
    # Ignore "+O" info
    if(grepl("_", i)){
      spbase <- gsub(".*\\(([mdt][0-9]+:[0-9]+)_.*", "\\1", i, perl = T)
      fa <- gsub(".*_([0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- c(spbase, fa)
    } else {
      spall <- gsub(".*\\(([mdt][0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- spall
    }
  } else if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                         "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                         "LPS", "ChE")){#will not contain "CL"
    # Not ignore e/p(O-/P-) connection, but lipidr may not use the e/p info
    if(grepl("_", i)){
      fainfo <- gsub(".*\\((.*)\\)", "\\1", i)
      fas <- strsplit(fainfo, "[_/]")
      m <- gregexpr("[0-9]+:[0-9]+(e|p)?", fas)
      fas <- unlist(regmatches(as.character(fas), m))
    } else {
      m <- gregexpr("[0-9]+:[0-9]+(e|p)?", i)
      fas <- unlist(regmatches(i, m))
    }
    fas <- gsub("(.*?)(e|p)", "\\2\\1", fas)
    if(length(fas) == 0){#the ChE may has no fas info
      fas <- "0:0"
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
##NOTE: Due to neglect of the "+O" symbol, lipidr may send out a "duplicate" warning
##Here, we ignore it
d <- add_sample_annotation(d, anno)
plot_samples(d, type = "tic", log = T)
#plot_molecules(d, "sd", measure = "Area", log = FALSE)
plot_lipidclass(d, "boxplot")
lipid_classes = rowData(d)$Class %in% c("Cer", "DG", "PE", "PG",
                                        "PC", "LPC", "PI", "PS", 
                                        "SM", "SPH", "TG")
d = d[lipid_classes, ]
## Normalization
d_summarized = summarize_transitions(d, method = "average")
d_normalized = normalize_pqn(d_summarized, measure = "Area", exclude = "blank", log = TRUE)
plot_samples(d_normalized, "boxplot")
## Analysis
mvaresults = mva(d_normalized, measure="Area", method="PCA")
plot_mva(mvaresults, color_by="SampleGroup", components = c(1,2))
de_results = de_analysis(
  data=d_normalized, 
  COS7OADGAT1.2i - COS7OADGAT,
  measure="Area"
)
significant_molecules(de_results)
plot_results_volcano(de_results, show.labels = F)
enrich_results = lsea(de_results, rank.by = "logFC")
significant_lipidsets(enrich_results)
plot_class_enrichment(de_results, significant_lipidsets(enrich_results))
plot_chain_distribution(de_results)
library(ggplot2)
ggsave("chains.pdf", width = 10, height = 20)



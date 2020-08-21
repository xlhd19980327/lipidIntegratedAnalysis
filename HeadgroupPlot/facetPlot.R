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

### Tidy for and do Visualization ###
## Use data_sub_allhandle(Delete duplication) to calculate itensity of each lipid class
data_sub_classSum_stat <- data_sub_allhandle %>%
  ungroup() %>%
  select(-ms1, -lipidName) %>%
  gather(key = "case", value = "lipidsum", -Class) %>%
  group_by(Class, case) %>%
  summarise(lipidsum = sum(lipidsum)) %>%
  mutate(group = allgroups[match(case, names(allgroups))]) 
data_sub_classSum_stat2 <- data_sub_classSum_stat %>%
  group_by(Class, group) %>%
  summarise(realmean = mean(lipidsum),
            sd = sd(lipidsum))
data_sub_classSum_p <- split(data_sub_classSum_stat, data_sub_classSum_stat$Class)
getPValue <- function(x){
  Class <- unique(x$Class)
  #Totally equal to t.test() one by one
  p <- pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value
  p1 <- p[rownames(p) == controlGrp, ]
  if(length(p1) != 0){
    names(p1) <- colnames(p)
  }
  p2 <- p[, colnames(p) == controlGrp]
  if(length(p2) != 0){
    names(p2) <- rownames(p)
  }
  p_tidy <- c(p1, p2)
  p_tidy <- p_tidy[!is.na(p_tidy)]
  if(length(p_tidy) == 0){
    result <- NULL
  }else {
    result <- data.frame(
      group = c(names(p_tidy), controlGrp),
      p = c(p_tidy, NA), 
      Class = Class
    )
  }
  return(result)
}
data_sub_classSum_stat3 <- lapply(data_sub_classSum_p, getPValue)
data_sub_classSum_stat3 <- do.call(rbind, data_sub_classSum_stat3)
addSigLabel <- function(x){
  x[is.na(x)] <- 1
  label <- c('', "*", "**", "***")
  allLable <- ifelse(x >= 0.05, label[1], NA)
  allLable <- ifelse(x < 0.05 & x >= 0.01, label[2], allLable)
  allLable <- ifelse(x < 0.01 & x >= 0.001, label[3], allLable)
  allLable <- ifelse(x < 0.001, label[4], allLable)
}
sigLabel <- addSigLabel(data_sub_classSum_stat3$p)
data_sub_classSum_stat3 <- cbind(data_sub_classSum_stat3, sigLabel = sigLabel)
data_sub_classSum_integStat <- left_join(data_sub_classSum_stat2, data_sub_classSum_stat3)
ggplot() +
  geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean, fill = group), stat = "identity") +
  geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd, color = group), width = 0.2) +
  geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
  geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+2*sd, label = sigLabel),
            size = 3, fontface = "bold", color = "red") +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        line = element_line(colour = "black", size = 1, 
                            linetype = 1, lineend = "butt")) +
  facet_wrap(~Class, scales="free") 

ggplot() +
  geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean), fill = "white", color = "black", stat = "identity", 
           width = 0.5, position=position_dodge(10)) +
  geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd), color = "black", width = 0.2) +
  geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center', 
               binwidth = 5) +
  geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+2*sd, label = sigLabel),
            size = 3, fontface = "bold", color = "red") +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        line = element_line(colour = "black", size = 1, 
                            linetype = 1, lineend = "butt"), 
        axis.title = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~Class, scales="free") 


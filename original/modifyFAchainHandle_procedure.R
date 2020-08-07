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
setwd("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/LMPDhsa/testData/zsy_DGATinhibitors/statistics")
data <- read.csv("Cos7_integ.csv", skip = 1)
#NOTE1-ref: may have the first character garbled
allgroups <- scan("Cos7_integ.csv", what = "character", nlines = 1, sep = ",", quote = "\"")
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
# Missing value will be replaced by 0
percent <- 0.67
data[is.na(data)] <- 0
data_intensity <- data[, -1:-notdataColsLen]
ind_ok <- apply(data_intensity == 0, 1, sum)/ncol(data_intensity) <= percent
data <- data[ind_ok, ]
# # Missing value will be replaced by half of the minimum positive value in the original data
# minInten <- min(data_intensity[data_intensity > 0], na.rm = T)/2
# data[data == 0] <- minInten
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
#Delete some "0:0" info
data_sub_all <- cbind(data_sub, 
                      lipidName = paste0(data_sub$Class, data_sub$FattyAcid))
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
write.csv(data_sub_allhandle, "data_sub_allhandle_cos7.csv", row.names = F)
## Normalize the data(PQN) (or not)
lipids <- data_sub_allhandle$lipidName
dataMatx <- subset(data_sub_allhandle, select = -lipidName)
dataMatx[dataMatx == 0] <- NA
## Uncomment the line below can be normalized
dataMatx <- dataMatx/apply(dataMatx/rowMeans(dataMatx, na.rm = TRUE), 2, median, na.rm = TRUE)
data_sub_allhandle <- cbind(lipidName = lipids, dataMatx, 
                            Class = gsub("(.*?)\\(.*", "\\1", data_sub_allhandle$lipidName))
## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
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
# # Tidy and integrate itensity of lipid class containing FA chain info
# # Use MS2 to do the later statistics only
# lipid_subclass_handle <- data.frame()
# for(i in 1:nrow(data_sub_allhandle)){
#   # e,p connection info will be ignored
#   fa <- fasInfo[[i]][[1]]
#   Class <- names(fasInfo[[i]][1])
#   subclass <- paste0(Class, "(", fa, ")")
#   lipid_subclass_handle <- rbind(lipid_subclass_handle, 
#                                  cbind(subclass = subclass, 
#                                        bind_rows(replicate(length(subclass), data_sub_allhandle[i, ], simplify = FALSE))))
# }
# lipid_subclass_handle <- subset(lipid_subclass_handle, subset = ms1 == F)
# write.csv(lipid_subclass_handle, "lipid_subclass_handle_cos7.csv", row.names = F)

# Tidy and integrate itensity of lipid class containing FA chain info
# Use MS1 & MS2 to do the later statistics, code above should be commented manually
lipid_subclass_handle <- data.frame()
for(i in 1:nrow(data_sub_allhandle)){
  # e,p connection info will be ignored
  fa <- fasInfo[[i]][[1]]
  chains <- sum(as.numeric(gsub(".*?([0-9]+):.*", "\\1", fa)))
  unsaturate <- sum(as.numeric(gsub(".*?:([0-9]+).*", "\\1", fa)))
  Class <- names(fasInfo[[i]][1])
  subclass <- paste0(Class, "(", chains, ":", unsaturate, ")")
  lipid_subclass_handle <- rbind(lipid_subclass_handle, 
                                 cbind(subclass = subclass, 
                                       bind_rows(replicate(length(subclass), data_sub_allhandle[i, ], simplify = FALSE))))
}
write.csv(lipid_subclass_handle, "lipid_subclass_handle_all_cos7.csv", row.names = F)

# ### Tidy for and do Visualization ###
# ## Use data_sub_allhandle(Delete duplication) to calculate itensity of each lipid class
# data_sub_classSum_stat <- data_sub_allhandle %>%
#   ungroup() %>%
#   select(-ms1, -lipidName) %>%
#   gather(key = "case", value = "lipidsum", -Class) %>%
#   group_by(Class, case) %>%
#   summarise(lipidsum = sum(lipidsum)) %>%
#   mutate(group = gsub("(.*?)_[0-9]+$", "\\1", case)) 
# data_sub_classSum_stat2 <- data_sub_classSum_stat %>%
#   group_by(Class, group) %>%
#   summarise(mean = mean(lipidsum) / n(), 
#             realmean = mean(lipidsum),
#             sd = sd(lipidsum))
# data_sub_classSum_p <- split(data_sub_classSum_stat, data_sub_classSum_stat$Class)
# getPValue <- function(x){
#   Class <- unique(x$Class)
#   #Totally equal to t.test() one by one
#   p <- pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value
#   p_tidy <- gather(rownames_to_column(as.data.frame(p), var = "group"), 
#                    -group, key = "cprGroup", value = "p")
#   p_tidy <- subset(p_tidy, cprGroup == controlGrp, select = -cprGroup) ## the control group
#   result <- cbind(p_tidy, Class = Class)
#   return(result)
# }
# data_sub_classSum_stat3 <- lapply(data_sub_classSum_p, getPValue)
# data_sub_classSum_stat3 <- do.call(rbind, data_sub_classSum_stat3)
# data_sub_classSum_integStat <- left_join(data_sub_classSum_stat, data_sub_classSum_stat2)
# data_sub_classSum_integStat <- left_join(data_sub_classSum_integStat, data_sub_classSum_stat3)
addSigLabel <- function(x){
  x[is.na(x)] <- 1
  label <- c('', "*", "**", "***")
  allLable <- ifelse(x >= 0.05, label[1], NA)
  allLable <- ifelse(x < 0.05 & x >= 0.01, label[2], allLable)
  allLable <- ifelse(x < 0.01 & x >= 0.001, label[3], allLable)
  allLable <- ifelse(x < 0.001, label[4], allLable)
}
# sigLabel <- addSigLabel(data_sub_classSum_integStat$p)
# data_sub_classSum_integStat <- cbind(data_sub_classSum_integStat, sigLabel = sigLabel)
# ggplot(data = data_sub_classSum_integStat, aes(x = group)) +
#   geom_bar(aes(y = mean, fill = group), stat = "identity") +
#   geom_errorbar(aes(ymin = realmean - sd, ymax = realmean + sd, color = group), width = 0.2) +
#   geom_dotplot(aes(y = lipidsum), binaxis='y', stackdir='center') +
#   geom_text(aes(y = realmean+1.5*sd, label = sigLabel),
#             size = 3, fontface = "bold", color = "red") +
#   theme_classic() +
#   facet_wrap(~Class, scales="free") 

## Use lipid_subclass_handle to calculate itensity of lipid class containing FA chain info(subclass)
lipid_subclass_stat <- lipid_subclass_handle %>%
  ungroup() %>%
  select(-ms1, -lipidName) %>%
  gather(key = "case", value = "lipidsum", -subclass, -Class) %>%
  group_by(subclass, case, Class) %>% 
  filter(!is.na(lipidsum)) %>%    #delete low abundance of the lipid signal
  summarise(lipidsum = sum(lipidsum) / n()) %>% #calc abundance of the lipid in a sample
  mutate(group = gsub("(.*?)_[0-9]+$", "\\1", case)) %>%
  ungroup() %>%
  group_by(group, subclass) %>%
  filter(n() >= 3) #if the lipid signal of one group are detected in no more than 2 samples, the lipid will be dropped
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
  if(nrow(p_tidy) == 0){
    result <- NULL
  }else {
    p_tidy <- subset(p_tidy, cprGroup == controlGrp, select = -cprGroup) ## the control group
    if(nrow(p_tidy) == 0){
      result <- NULL
    }else {
      result <- cbind(p_tidy, subclass = subclass)
    }
  }
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
    ggsave(paste0("./lipsubclassFigures_facet9/", i, ".pdf"), dpi = 300, width = 15, height = 10)
  }
}

### Visualize with the tile plot
realmean <- lipid_subclass_stat2$realmean
names(realmean) <- lipid_subclass_stat2$group
getlog2FC <- function(x){
  cpr <- x[names(x) == controlGrp]
  if(length(cpr) != 0){
    res <- log2(x/cpr)
  }else{
    res <- rep(NA, length(x))
  }
}
log2FC <- unlist(tapply(realmean, lipid_subclass_stat2$subclass, getlog2FC))
lipid_subclass_stat_tile <- lipid_subclass_stat2 %>%
  add_column(log2FC = log2FC) %>%
  #ignore other info(eg. e/p, d)
  mutate(chain = as.numeric(gsub(".*?([0-9]+):.*", "\\1", subclass)), 
         unsaturate = as.numeric(gsub(".*?:([0-9]+).*", "\\1", subclass)), 
         Class = gsub("(.*?)\\(.*", "\\1", subclass))
for(i in groupsLevel[groupsLevel != controlGrp]){
  oneGrpdata <- subset(lipid_subclass_stat_tile, 
                       subset = group == i & !is.na(log2FC))
  ggplot(oneGrpdata, aes(chain, unsaturate, fill = log2FC)) + 
    geom_tile() + 
    facet_wrap(~Class) + 
    scale_x_continuous(breaks = seq(0, max(oneGrpdata$chain), 2))+
    xlab("Chain length") + 
    ylab("Chain unsaturation") + 
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = 0) #+ 
    #geom_text(aes(label = ifelse(nmolecules > 1, nmolecules, 
    #                             "")))
  ggsave(paste0("./lipsubclassFigures_facet10/", i, ".pdf"), dpi = 300, width = 25, height = 10)
}

### Get different chain level and unsaturated level visualization
## Only use MS2 data to do statistics
lipid_subclass_stat_cinfo <- lipid_subclass_handle %>%
  ungroup() %>%
  select(-ms1, -lipidName) %>%
  gather(key = "case", value = "lipidsum", -subclass, -Class) %>%
  #ignore other info(eg. e/p, d)
  mutate(chain = gsub(".*?([0-9]+):.*", "\\1", subclass), 
         #unsaturate = gsub(".*?:([0-9]+).*", "\\1", subclass),
         #HG_unsaturate = paste0(Class, "_", unsaturate), 
         HG_chain = paste0(Class, "_", chain)) %>%
  group_by(HG_chain, case, Class) %>% 
  filter(!is.na(lipidsum)) %>%    #delete low abundance of the lipid signal
  summarise(lipidsum = sum(lipidsum) / n()) %>% #calc abundance of the lipid in a sample
  mutate(group = gsub("(.*?)_[0-9]+$", "\\1", case)) %>%
  ungroup() %>%
  group_by(group, HG_chain) %>%
  filter(n() >= 3) #if the lipid signal of one group are detected in no more than 2 samples, the lipid will be dropped
lipid_subclass_stat2_cinfo <- lipid_subclass_stat_cinfo %>%
  group_by(HG_chain, group) %>%
  summarise(mean = mean(lipidsum) / n(), 
            realmean = mean(lipidsum) ,
            sd = sd(lipidsum))
lipid_subclass_stat_cinfo_p <- split(lipid_subclass_stat_cinfo, lipid_subclass_stat_cinfo$HG_chain)
getPValue_cs <- function(x, type = c("chain", "unsaturate")){
  type <- paste0("HG_", type)
  subclass <- unique(x[[type]])
  #Totally equal to t.test() one by one
  p <- pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value
  p_tidy <- gather(rownames_to_column(as.data.frame(p), var = "group"), 
                   -group, key = "cprGroup", value = "p")
  if(nrow(p_tidy) == 0){
    result <- NULL
  }else {
    p_tidy <- subset(p_tidy, cprGroup == controlGrp, select = -cprGroup) ## the control group
    if(nrow(p_tidy) == 0){
      result <- NULL
    }else {
      result <- cbind(p_tidy, subclass)
      colnames(result)[ncol(result)] <- type
    }
  }
  return(result)
}
lipid_subclass_stat3_c <- lapply(lipid_subclass_stat_cinfo_p, getPValue_cs, "chain")
lipid_subclass_stat3_c <- do.call(rbind, lipid_subclass_stat3_c)
lipid_subclass_integStat_c <- left_join(lipid_subclass_stat_cinfo, lipid_subclass_stat2_cinfo)
lipid_subclass_integStat_c <- left_join(lipid_subclass_integStat_c, lipid_subclass_stat3_c)
sigLabel3 <- addSigLabel(lipid_subclass_integStat_c$p)
lipid_subclass_integStat_c <- cbind(lipid_subclass_integStat_c, sigLabel = sigLabel3)
for(i in lipClasses){
  oneLipClassData <- subset(lipid_subclass_integStat_c, 
                            subset = Class == i)
  if(nrow(oneLipClassData) != 0){
    ggplot(data = oneLipClassData, aes(x = group)) +
      geom_bar(aes(y = mean, fill = group), stat = "identity") +
      geom_errorbar(aes(ymin = realmean - sd, ymax = realmean + sd, color = group), width = 0.2) +
      geom_dotplot(aes(y = lipidsum), binaxis='y', stackdir='center') +
      geom_text(aes(y = realmean+1.5*sd, label = sigLabel),
                size = 3, fontface = "bold", color = "red") +
      theme_classic() +
      facet_wrap(~HG_chain, scales="free") 
    ggsave(paste0("./lipsubclassFigures_facet4/", i, ".pdf"), dpi = 300, width = 15, height = 10)
  }
}










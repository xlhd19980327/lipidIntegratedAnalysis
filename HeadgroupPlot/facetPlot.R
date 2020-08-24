### For LipidSearch data ###

## Source will offer the following variables:
## ??allgroups, groupsLevel, mSet
source("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/procedure/github/lipidpreproc/LipidSearchData_preprocess.R")

data_tidy <- as.data.frame(t(mSet[["dataSet"]][["preproc"]])) %>%
  rownames_to_column(var = "lipidName") %>%
  mutate(Class = gsub("(.*?)\\(.*", "\\1", lipidName))

##!!!!!WARNING: Class statistics will only contain the following lipid class:
## Cer, SM, SPH, FA, MG, DG, TG, PA, PC, PE, PG, PI, PS, LPA, LPC, LPE, LPG, LPI, LPS, ChE

## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
lipids <- data_tidy$lipidName
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
  }else{
    result <- list(otherClass = Class, ms1 = NA)
    return(result)
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
data_tidy <- cbind(data_tidy, ms1 = ms1Info)
data_sub_ms1 <- subset(data_tidy, subset = ms1 == T)
data_sub_ms2 <- subset(data_tidy, subset = ms1 == F)
# Tidy and integrate itensity of lipid class containing FA chain info
# Use MS2 to do the later statistics only
lipid_subclass_handle <- data.frame()
for(i in 1:nrow(data_tidy)){
  # e,p connection info will be ignored
  fa <- fasInfo[[i]][[1]]
  Class <- names(fasInfo[[i]][1])
  subclass <- paste0(Class, "(", fa, ")")
  lipid_subclass_handle <- rbind(lipid_subclass_handle, 
                                 cbind(subclass = subclass, 
                                       bind_rows(replicate(length(subclass), data_tidy[i, ], simplify = FALSE))))
}
lipid_subclass_handle <- subset(lipid_subclass_handle, subset = ms1 == F)

### Tidy for and do Visualization ###
## Use data_tidy(Delete duplication) to calculate itensity of each lipid class
data_sub_classSum_stat <- data_tidy %>%
  ungroup() %>%
  select(-ms1, -lipidName) %>%
  gather(key = "case", value = "lipidsum", -Class) %>%
  group_by(Class, case) %>%
  summarise(lipidsum = sum(lipidsum, na.rm = T)) %>%
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


# Use two sample comparation once
lipsample <- "DGAT1i+2i"
fileLoc = "./testData/zsy_DGATinhibitors/HeLaData/output/"
source("./utilityFunc/addLibnames.R")
#source("./utilityFunc/FAchainStat.R")

## Data preparation
# Tidy the lipid_subclass_integStat and retain the useful info
#lipid_subclass_integStat <- FAchainStat(..., stat = F)
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
# Lipid reaction library file input
lipReact <- read.csv("./testData/zsy_DGATinhibitors/cos7Data/input/hsa_lipidreact.csv")
#Just need contain these colums: "reaction_id", "gene_symbol"
allReact <- read.csv("./testData/zsy_DGATinhibitors/cos7Data/input/hsa_all_integData.csv")  
#use unique/distinct, other info will not consider here
allReact <- unique(subset(allReact, select = c("reaction_id", "gene_symbol")))
integReact <- inner_join(lipReact, allReact, by = "reaction_id")
integReact <- subset(integReact, 
                     select = c("reaction_id", "gene_symbol", "R_lipid", "P_lipid"))

# !!!!!WARNING: Some libnames will be wrong, cause carbon base info add
lipidsubclass <- cbind(lipid_subclass_tidyStat, 
                       libnames = apply(lipid_subclass_tidyStat, 1, addLibnames))
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
write.csv(regState, 
          paste0(fileLoc, "regStat_", lipsample, ".csv"), row.names = F)

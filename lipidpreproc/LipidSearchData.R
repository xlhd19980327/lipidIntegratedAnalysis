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
nmin <- min(nsamples)
percent <- 1-1/nmin
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
## Duplication handle
# Do not do any subset
data_sub <- data
data_sub_all <- cbind(data_sub, 
                      lipidName = paste0(data_sub$Class, data_sub$FattyAcid))
data_sub_dup <- data_sub_all %>%
  group_by(lipidName) %>%
  filter(n() > 1)
data_sub_dup_list <- split(data_sub_dup, data_sub_dup$lipidName)
# Rule: duplicated lipids in different modes- chooese stronger intensity
#       duplicated lipids in the same modes but have different adduct - choose stronger intensity
#       i.e. duplicated lipids all choose stronger intensity
getSoloInten <- function(x){
  data <- x[, (notdataColsLen+1):(notdataColsLen+length(allgroups))]
  result <- rep(0, length(allgroups))
  for(i in 1:nrow(data)){
    one <- unlist(data[i, ])
    result <- ifelse(one > result, one, result)
  }
  return(result)
}
data_sub_dup_soloInten <- sapply(data_sub_dup_list, getSoloInten)
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




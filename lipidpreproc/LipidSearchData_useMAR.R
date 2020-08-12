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
library(MetaboAnalystR)
options(stringsAsFactors = F)
setwd("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/procedure/github/testData/zsy_DGATinhibitors/cos7Data/input")
data <- read.csv("Cos7_integ.csv", skip = 1)
#NOTE1-ref: may have the first character garbled
allgroups <- scan("Cos7_integ.csv", what = "character", nlines = 1, sep = ",", quote = "\"")
notdataColsLen <- sum(allgroups == '')
allgroups <- allgroups[allgroups != '']
groupsLevel <- unique(allgroups)
nsamples <- table(factor(allgroups, levels = groupsLevel))
#NOTE2-ref: should indicate the control group or set on the first column
controlGrp <- groupsLevel[1]

### Check Data Integrity ###
if(min(nsamples) < 3){
  stop("At least one group have no more than 2 replicates, PROGRAM EXIT!")
}

### Formatting data ###
## Change 0 to NA
data[data == 0] <- NA

### Clean and Tidy the data ###
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

### Input data into MetaboAnalystR ###
data <- data_sub_allhandle[, -ncol(data_sub_allhandle)]
mSet<-InitDataObjects("conc", "stat", FALSE)
## The code following is adapted from MetaboAnalystR::Read.TextData()
mSet$dataSet$cls.type <- "disc"
mSet$dataSet$format <- "colu"
var.nms <- data$lipidName
data <- data[, -1]
smpl.nms <- colnames(data)
cls.lbl <- allgroups
conc <- t(data)
mSet$dataSet$type.cls.lbl <- class(cls.lbl)
orig.var.nms <- var.nms
names(orig.var.nms) <- var.nms
rownames(conc) <- smpl.nms
colnames(conc) <- var.nms
mSet$dataSet$orig.cls <- mSet$dataSet$cls <- as.factor(as.character(cls.lbl))
mSet$dataSet$cmpd <- var.nms
mSet$dataSet$mumType <- "table"
mSet$dataSet$orig.var.nms <- orig.var.nms
mSet$dataSet$orig <- conc
mSet$msgSet$read.msg <- c(paste("The uploaded data file contains ", 
                                     nrow(conc), " (samples) by ", ncol(conc), " (", tolower(GetVariableLabel(mSet$dataSet$type)), 
                                     ") data matrix.", sep = ""))
## End of rewriting from MetaboAnalystR::Read.TextData()
mSet<-SanityCheckData(mSet)
#Retain the ':' or '_' character(i.e. the original characters)
mynames <- mSet[["dataSet"]][["orig.var.nms"]]
mSet[["dataSet"]][["cmpd"]] <- mynames
names(mSet[["dataSet"]][["orig.var.nms"]]) <- mynames
## Delete NA (or 0) value related argument 
# Default methods: 
# 1. Remove missing value
# The lipid is not NA (or 0) at least one sample in a group remains
# eg. 3 samples in a group -- select <= 0.67(2/3) missing value
# eg. 5 samples in a group -- select <= 0.8(4/5) missing value
# Default will use this method to remove NA data
# 2. Imputation of missing value
# MetaboAnalystR offer many methods to do imputation
# Default will use a small value(min(data)/2) to replace NA
nmin <- min(nsamples)
percent <- 1-1/nmin
handleMissingData <- function(data, remove = T, imput = "min"){
  if(remove == T){
    data<-RemoveMissingPercent(data, percent=percent)
  }
  data<-ImputeVar(data, method=imput)
  return(data)
}
#!!!Client options
mSet <- handleMissingData(mSet, 
                          remove = T, imput = "min")
## Normalization / Scaling
# Default methods: 
# 1. Normalization method: Probabilistic Quotient Normalization(PQN) without using a reference sample
# 2. Transformation method: None
# 3. Scaling method: Auto scaling(mean-centered and divided by the standard deviation of each variable)
mSet<-PreparePrenormData(mSet)
#!!!Client options
mSet<-Normalization(mSet, 
                    rowNorm = "ProbNormT", transNorm = "NULL", scaleNorm = "AutoNorm", ref = NULL, 
                    ratio=FALSE, ratioNum=20)




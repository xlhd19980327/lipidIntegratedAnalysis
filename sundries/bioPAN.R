source("./utilityFunc/readingLipidData.R")
source("./utilityFunc/MARpreproc.R")

outputLoc <- "./testData/052921sffData_CCl4liver/output/"
prepDataSet <- function(x, dataset = dataSet){
  ind <- dataset$allgroups %in% c(x, dataset$controlGrp)
  ind2 <- dataset$groupsLevel %in% c(x, dataset$controlGrp)
  dataset$data <- dataset$data[, ind]
  dataset$groupsLevel <- dataset$groupsLevel[ind2]
  dataset$allgroups <- dataset$allgroups[ind]
  return(dataset)
}

analOpt <- "48hrs"
dataSet <- readingLipidData(datafile = "./testData/052921sffData_CCl4liver/input/data.csv",
                            sampleList = "./testData/052921sffData_CCl4liver/input/des.csv", 
                            controlGrp = "24hrs", dataType = "Lipids", delOddChainOpt = T)
cat(paste0(analOpt, " will be analyzed with ", dataSet$controlGrp, "\n"))
dataset <- prepDataSet(analOpt)
mSet <- MARpreproc(dataSet = dataset, fileLoc = outputLoc, perc = 2/3*100)
data_type <- dataset$dataType
data_proc <- t(mSet[["dataSet"]][["proc"]])
dataset$lipidName <- rownames(data_proc)
rownames(data_proc) <- NULL
data_proc <- as.data.frame(data_proc)
dataset$data <- data_proc

dataSet <- dataset
#Always consider e/p bond in BioPAN will make results more precise
ignore <- F
allgroups <- dataSet$allgroups
controlGrp <- dataSet$controlGrp
groupsLevel <- dataSet$groupsLevel
dataType <- dataSet$dataType

## Source will offer the following contents:
## Function(s): getFAsInfo
if(dataType == "LipidSearch"){
  source("./utilityFunc/getFAsInfo.R")
}
if(dataType == "MS_DIAL"){
  source("./utilityFunc/getFAsInfo_msdial.R")
}
## Source will offer the following contents:
## Function(s): getClassInfo
source("./utilityFunc/getClassInfo.R")

data_tidy <- dataSet[["data"]] %>%
  mutate(lipidName = dataSet[["lipidName"]], 
         Class = switch(dataSet$dataType,
                        LipidSearch = sapply(lipidName, getClassInfo, "LipidSearch", ignore = ignore), 
                        MS_DIAL = sapply(lipidName, getClassInfo, "MS_DIAL", ignore = ignore)))
lipids <- data_tidy$lipidName
# Get the FA info and MS1 info(all_info)
fasInfo <- lapply(lipids, getFAsInfo, ignore)
## Calculate the sum of the lipid subchain and put the info into 'subclass' colum
lipid_subclass_handle <- data.frame()
for(i in 1:nrow(data_tidy)){
  # e,p connection info will be ignored
  fa <- fasInfo[[i]][[1]]
  Class <- names(fasInfo[[i]][1])
  if(any(fa %in% "UnknownPattern")){
    subclass <- paste0(Class, "(", fa, ")")
  }else{
    chains <- sum(as.numeric(gsub(".*?([0-9]+):.*", "\\1", fa)))
    unsaturate <- sum(as.numeric(gsub(".*?:([0-9]+).*", "\\1", fa)))
    subclass <- paste0(Class, "(", chains, ":", unsaturate, ")")
  }
  
  lipid_subclass_handle <- rbind(lipid_subclass_handle, 
                                 cbind(subclass = subclass, 
                                       bind_rows(replicate(length(subclass), data_tidy[i, ], simplify = FALSE))))
}
write.csv(lipid_subclass_handle, paste0(outputLoc, "out_subclass_all_info.csv"), 
          row.names = F)

## Duplication handle(sum of all the duplications)
## O/P names handle
## Colnames handle to fulfill BioPAN file requirement
lipsubclass <- lipid_subclass_handle$subclass
lipsubclass <- gsub("(.*?)\\((O|P)\\)(.*)", "\\2\\-\\1\\3", lipsubclass)
lipid_subclass_handle$subclass <- lipsubclass
lipid_subclass_all <- lipid_subclass_handle %>%
  ungroup() %>%
  select(-Class, -lipidName) %>%
  filter(!grepl("UnknownPattern", subclass)) %>%
  gather(key = "case", value = "lipidsum", -subclass) %>%
  group_by(subclass, case) %>% 
  filter(!is.na(lipidsum)) %>%    #delete low abundance of the lipid signal
  summarise(lipidsum = sum(lipidsum)) %>%
  spread(key = case, value = lipidsum)
sortedgroups <- sapply(colnames(lipid_subclass_all)[2:ncol(lipid_subclass_all)], 
                       function(x) allgroups[[x]])
nsamples <- table(factor(allgroups, levels = groupsLevel))
labedgroups <- lapply(sort(groupsLevel), 
                      function(x) paste0(x, "_", 1:nsamples[[x]]))
labedgroups <- do.call(c, labedgroups)
colnames(lipid_subclass_all) <- c("subclass", 
                                  labedgroups[order(order(sortedgroups))])
write.csv(lipid_subclass_all, paste0(outputLoc, "out_forBioPAN.csv"), 
          row.names = F)

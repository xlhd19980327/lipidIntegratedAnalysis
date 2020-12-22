MARpreproc <- function(dataSet, fileLoc){
  data <- dataSet$data
  lipidName <- dataSet$lipidName
  allgroups <- dataSet$allgroups
  nsamples <- dataSet$nsamples
  
  
  ### Input data into MetaboAnalystR ###
  mSet<-InitDataObjects("conc", "stat", FALSE)
  ## The code following is adapted from MetaboAnalystR::Read.TextData()
  mSet$dataSet$cls.type <- "disc"
  mSet$dataSet$format <- "colu"
  var.nms <- lipidName
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
  qs::qsave(conc, file = "data_orig.qs")
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
    data<-ImputeMissingVar(data, method=imput)
    return(data)
  }
  #!!!Client options
  mSet <- handleMissingData(mSet, 
                            remove = T, imput = "min")
  ## Data filtering
  # The purpose of the data filtering is to identify and remove variables that are unlikely to be of use when modeling the data. 
  # Default methods:
  # None for targeted lipidomics and most other lipidomics data
  #!!!Client options
  mSet<-FilterVariable(mSet, 
                       filter = "none", qcFilter = "F", rsd = 25)
  if(!is.na(fileLoc)){
    write.csv(t(mSet[["dataSet"]][["proc"]]), paste0(fileLoc, "data_tidy.csv"))
  }
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
  return(mSet)
}

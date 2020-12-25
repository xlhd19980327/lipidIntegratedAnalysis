readingLipidData <- function(datafile, sampleList, controlGrp = "", dataType, 
                             lipField = NA, delOddChainOpt = F, na.char = NULL){
  if(dataType == "LipidSearch"){
    if(is.na(lipField)){
      lipField <- "LipidIon"
    }
  }else if(dataType == "MS_DIAL"){
    if(is.na(lipField)){
      lipField <- "Metabolite name"
    }
  }else{
    if(is.na(lipField)){
      stop("You should indicate the feature colum instead of NA, PROGRAM EXIT!")
    }
  }
  #!!!!!WARNING: The NA string may be others, may add the char clients customized 
  data <- read.csv(datafile, na.strings = c("N/A", "NA", na.char))
  firstline <- scan(datafile, what = "character", nlines = 1, sep = ",", quote = "\"", 
                    na.strings = c("N/A", "NA"))
  colnames(data) <- firstline
  if(!lipField %in% colnames(data)){
    stop("You should indicate the right feature colum, PROGRAM EXIT!")
  }
  sampleInfo <- read.csv(sampleList)
  allgroups <- sampleInfo$conditions
  sampleInd <- colnames(data) %in% sampleInfo$samples 
  groupsLevel <- unique(allgroups)
  nsamples <- table(factor(allgroups, levels = groupsLevel))
  names(allgroups) <- sampleInfo$samples
  #NOTE2-ref: should indicate the control group or set on the first row in sampleList
  if(controlGrp == ""){
    controlGrp <- groupsLevel[1]
  }
  
  ### Check Data Integrity ###
  # if(min(nsamples) < 3){
  #   stop("At least one group have no more than 2 replicates, PROGRAM EXIT!")
  # }
  if(any(!(apply(data[, sampleInd], 2, is.numeric)))){
    stop("Data contain non-numeric variable. You may check the NA string and offer na.char parameter. PROGRAM EXIT!")
  }
  
  ### Formatting data ###
  ## Change 0 to NA
  data[data == 0] <- NA
  
  ### Clean and Tidy the data ###
  ## Delete odd FA chain lipids
  delOddChain <- function(x, 
                          delOddChainOpt_in = delOddChainOpt){
    if(delOddChainOpt_in){
      fas <- switch(dataType, 
                    LipidSearch = gsub(".*?(\\(.*\\))[\\+\\-].*", "\\1", data[[lipField]]), 
                    MS_DIAL = data[[lipField]])
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
      res <- x[ind, ]
    }else{
      res <- x
    }
    return(res)
  }
  data <- delOddChain(data)
  ## Duplication handle
  data_sub_all <- cbind(data[, sampleInd], 
                        lipidName = switch(dataType, 
                                           LipidSearch = gsub("(.*\\))[\\+\\-].*", "\\1", data[[lipField]]), 
                                           #Delete some "0:0" info
                                           MS_DIAL = gsub("(.*)\\/0:0$", "\\1", data[[lipField]]), 
                                           Metabolites = data[[lipField]], 
                                           Proteins = data[[lipField]]))
  data_sub_dup <- data_sub_all %>%
    group_by(lipidName) %>%
    filter(n() > 1)
  data_sub_sing_allhandle <- data_sub_all %>%
    group_by(lipidName) %>%
    filter(n() == 1) 
  if(nrow(data_sub_dup) == 0){
    data_sub_allhandle <- data_sub_sing_allhandle
  }else{
    data_sub_dup_list <- split(data_sub_dup, data_sub_dup$lipidName)
    # Rule: duplicated lipids in different modes- chooese stronger intensity
    #       duplicated lipids in the same modes but have different adduct - choose stronger intensity
    #       i.e. duplicated lipids all choose stronger intensity
    #       The signal have the most sum intensity
    getSoloInten <- function(x){
      data <- x[, 1:(length(allgroups))]
      sigSums <- rowSums(data)
      ind <- match(max(sigSums), sigSums)
      result <- x[ind, ]
      return(result)
    }
    data_sub_dup_soloInten <- lapply(data_sub_dup_list, getSoloInten)
    data_sub_dup_allhandle <- do.call(rbind, data_sub_dup_soloInten)
    data_sub_allhandle <- bind_rows(data_sub_sing_allhandle, data_sub_dup_allhandle)
  }
  
  data <- data_sub_allhandle
  lipidName <- data$lipidName
  data <- subset(data, select = -lipidName)
  data_output <- cbind(lipidName = lipidName, data)
  
  
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
  
  return(mSet[["dataSet"]][["norm"]])
}

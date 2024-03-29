readingLipidData <- function(datafile, controlGrp = "", dataType, 
                             lipField = NA, delOddChainOpt = T, fileLoc, na.char = NULL){
  source("./branch/benchmark/code_extra/handleLipidName_HCC.R")
  
  if(dataType == "LipidSearch"){
    if(is.na(lipField)){
      lipField <- "LipidIon"
    }
  }
  if(dataType == "MS_DIAL"){
    if(is.na(lipField)){
      lipField <- "Metabolite.name"
    }
  }
  #!!!!!WARNING: The NA string may be others, may add the char clients customized 
  data <- read.csv(datafile, skip = 1, na.strings = c("N/A", "NA", na.char))
  #NOTE1-ref: may have the first character garbled
  firstline <- scan(datafile, what = "character", nlines = 1, sep = ",", quote = "\"", 
                    na.strings = c("N/A", "NA"))
  #notdataColsLen <- sum(firstline == '')
  allgroups <- firstline[firstline != '']
  sampleInd <- which(firstline != '')
  groupsLevel <- unique(allgroups)
  nsamples <- table(factor(allgroups, levels = groupsLevel))
  names(allgroups) <- colnames(data)[sampleInd]
  #NOTE2-ref: should indicate the control group or set on the first column
  if(controlGrp == ""){
    controlGrp <- groupsLevel[1]
  }
  
  
  ### Check Data Integrity ###
  if(min(nsamples) < 3){
    stop("At least one group have no more than 2 replicates, PROGRAM EXIT!")
  }
  if(any(!(apply(data[, which(firstline != "")], 2, is.numeric)))){
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
                                           HCC = sapply(data[[lipField]], handleLipidName_HCC)))
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
  if(!is.na(fileLoc)){
    write.csv(data_output, paste0(fileLoc, "data_tidy.csv"), row.names = F)
  }
  return(list(
    data = data, lipidName = lipidName, 
    groupsLevel = groupsLevel, allgroups = allgroups, controlGrp = controlGrp, nsamples = nsamples, 
    dataType = dataType
  ))
}

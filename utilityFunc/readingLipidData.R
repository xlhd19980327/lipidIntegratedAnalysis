readingLipidData <- function(datafile, controlGrp){
  data <- read.csv(datafile, skip = 1)
  #NOTE1-ref: may have the first character garbled
  allgroups <- scan(datafile, what = "character", nlines = 1, sep = ",", quote = "\"")
  notdataColsLen <- sum(allgroups == '')
  allgroups <- allgroups[allgroups != '']
  groupsLevel <- unique(allgroups)
  nsamples <- table(factor(allgroups, levels = groupsLevel))
  names(allgroups) <- colnames(data)[(notdataColsLen+1):(notdataColsLen+length(allgroups))]
  #NOTE2-ref: should indicate the control group or set on the first column
  if(controlGrp == ""){
    controlGrp <- groupsLevel[1]
  }
  
  
  ### Check Data Integrity ###
  if(min(nsamples) < 3){
    stop("At least one group have no more than 2 replicates, PROGRAM EXIT!")
  }
  
  ### Formatting data ###
  ## Change 0 to NA
  data[data == 0] <- NA
  
  ### Clean and Tidy the data ###
  ## Delete odd FA chain lipids
  #!!!Client options: Client should choose whether to delete odd FA chain lipid signal, default TRUE
  delOddChainOpt <- T
  delOddChain <- function(x, 
                          delOddChainOpt = T){
    if(delOddChainOpt){
      fas <- x$FattyAcid
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
  data_sub_all <- cbind(data, 
                        lipidName = paste0(data$Class, data$FattyAcid))
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
  data_sub_sing_allhandle <- data_sub_all[, -1:-notdataColsLen] %>%
    group_by(lipidName) %>%
    filter(n() == 1) 
  data_sub_sing_allhandle <- data_sub_sing_allhandle[, c((length(allgroups)+1), 1:length(allgroups))]
  data_sub_dup_allhandle <- rownames_to_column(as.data.frame(t(data_sub_dup_soloInten)))
  colnames(data_sub_dup_allhandle) <- colnames(data_sub_sing_allhandle)
  data_sub_allhandle <- bind_rows(data_sub_sing_allhandle, data_sub_dup_allhandle)
  
  data <- data_sub_allhandle
  lipidName <- data$lipidName
  data <- subset(data, select = -lipidName)
  return(list(
    data = data, lipidName = lipidName, 
    groupsLevel = groupsLevel, allgroups = allgroups, controlGrp = controlGrp, nsamples = nsamples
  ))
}

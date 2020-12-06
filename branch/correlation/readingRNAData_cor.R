readingRNAData <- function(datafile, sampleList, controlGrp = ""){
  if(grepl(".*\\.tsv$", datafile, ignore.case = T)){
    firstline <- scan(datafile, what = "character", nlines = 1, sep = "", quote = "", 
                      na.strings = c("N/A", "NA"))
    data <- read.table(datafile,
                     header=T, com='', 
                     quote='',check.names=F, sep="")
  }
  if(grepl(".*\\.csv$", datafile, ignore.case = T)){
    firstline <- scan(datafile, what = "character", nlines = 1, sep = ",", quote = "\"", 
                      na.strings = c("N/A", "NA"))
    data <- read.csv(datafile)
  }
  sampleInfo <- read.csv(sampleList)
  if(controlGrp == ""){
    controlGrp <- sampleInfo$conditions[1]
  }
  if(any(duplicated(data[, 1]))){
    cat("Duplicate gene identifier detected! Choose stronger intensity for analysis.\n")
    data_list <- split(data, data[,1])
    # Rule: duplicated gene symbols will choose stronger intensity
    #       The signal have the most sum intensity
    getSoloInten <- function(x){
      if(ncol(x) > 1){
        data <- x[, match(sampleInfo$samples, firstline)]
        sigSums <- rowSums(data)
        ind <- match(max(sigSums), sigSums)
        result <- x[ind, ]
      }else{
        result <- x
      }
      return(result)
    }
    data_soloInten <- lapply(data_list, getSoloInten)
    data <- do.call(rbind, data_soloInten)
  }
  rownames(data) <- data[, 1]
  data <- data[, -1]

  colnames(data) <- firstline[-1]
  
  data <- select(data, sampleInfo$samples)
  #data <- data[, match(sampleInfo$samples, firstline)-1]
  groupsLevel <- unique(sampleInfo$conditions)
  sampleInfo$conditions <- factor(sampleInfo$conditions, 
                                  levels = c(groupsLevel[groupsLevel == controlGrp], groupsLevel[groupsLevel != controlGrp]))
  
  return(list(
    data = data, sampleInfo = sampleInfo, 
    groupsLevel = groupsLevel, controlGrp = controlGrp
  ))
}

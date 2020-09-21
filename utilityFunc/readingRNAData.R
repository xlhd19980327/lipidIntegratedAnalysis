readingRNAData <- function(datafile, sampleList, controlGrp = ""){
  data <- read.table(datafile,
                     header=T, row.names=1, com='', 
                     quote='',check.names=F, sep="")
  sampleInfo <- read.csv(sampleList)
  if(controlGrp == ""){
    controlGrp <- sampleInfo$conditions[1]
  }
  data <- select(data, sampleInfo$samples)
  groupsLevel <- unique(sampleInfo$conditions)
  sampleInfo$conditions <- factor(sampleInfo$conditions, 
                                  levels = c(groupsLevel[groupsLevel == controlGrp], groupsLevel[groupsLevel != controlGrp]))
  
  ## Delete low gene abundance feature
  data <- data[rowSums(data)>2, ]
  return(list(
    data = data, sampleInfo = sampleInfo, 
    groupsLevel = groupsLevel, controlGrp = controlGrp
  ))
}

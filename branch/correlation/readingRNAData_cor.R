readingRNAData <- function(datafile, sampleList, controlGrp = "", type){
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
  
  if(type == "RNAseq"){
    ## Delete low gene abundance feature
    data <- data[rowSums(data)>2, ]
    batremvOpt <- ifelse("batch" %in% colnames(sampleInfo), T, F)
    designFormula <- eval(parse(text = ifelse(batremvOpt, "~ conditions + batch", "~ conditions")))
    ddsFullCountTable <- tryCatch(DESeqDataSetFromMatrix(countData = data,
                                                         colData = sampleInfo, 
                                                         design= designFormula), 
                                  error = function(e){if(conditionMessage(e) == "some values in assay are not integers"){
                                    cat("Convert some float to integer!\n")
                                    data <- round(data)
                                    DESeqDataSetFromMatrix(countData = data,
                                                           colData = sampleInfo, 
                                                           design= designFormula)
                                  }else{
                                    stop(paste0("Error in read counts input: ", conditionMessage(e), 
                                                ". PROGRAM EXIT!"))
                                  }})
    # dds <- tryCatch(DESeq(ddsFullCountTable), 
    #                 error = function(e){if(grepl("all gene-wise dispersion estimates are within 2 orders of magnitude", 
    #                                              conditionMessage(e))){
    #                   cat("Little variance among groups detect! Using proper data to do test is recommanded.\n")
    #                   dds <- estimateDispersionsGeneEst(ddsFullCountTable)
    #                   dispersions(dds) <- mcols(dds)$dispGeneEst
    #                   DESeq(dds)
    #                 } else{
    #                   stop(paste0("Error in DESeq: ", conditionMessage(e), 
    #                               ". PROGRAM EXIT!"))
    #                 }})
    dds <- DESeq(ddsFullCountTable)
    #Normalization(Note: only providing counts scaled by size or normalization factors, 
    #experiment design info only used in the other normalization methods.
    #In our code, this only scale data by size factors, which account for differences in sequencing depth)
    normalized_counts <- counts(dds, normalized=TRUE)
    
    return(normalized_counts)
  }
  if(type == "MiAr"){
    eset <- data
    
    eset <- normalizeBetweenArrays(eset)
    
    return(eset)
  }

}

DESeq2preproc <- function(dataSet = dataSet_RNA, 
                          fileLoc){
  data <- dataSet$data
  sampleInfo <- dataSet$sampleInfo
  
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
  dds <- DESeq(ddsFullCountTable)
  #Normalization(Note: only providing counts scaled by size or normalization factors, 
  #experiment design info only used in the other normalization methods.
  #In our code, this only scale data by size factors, which account for differences in sequencing depth)
  normalized_counts <- counts(dds, normalized=TRUE)
  #Sorting by mad(Median Absolute Deviation) value(Genes with greater 'overall' differences are ranked higher)
  normalized_counts_mad <- apply(normalized_counts, 1, mad)
  normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  write.csv(normalized_counts, paste0(fileLoc, "normalized_counts_sorted.csv"))
  #log2 transformation and normalized(may use "vst" instead)
  rld <- rlog(dds, blind=FALSE)
  if(batremvOpt){
    ## Using limma::removeBatchEffect to remove batch effect, it will change the original count matrix
    assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)
  }
  return(list(
    rld = rld, dds = dds, normalized_counts = normalized_counts, 
    dataSet = dataSet
  ))
}

DEAnalysis <- function(dataProc, experGrp, fileLoc){
  dds <- dataProc$dds
  controlGrp <- dataProc$dataSet$controlGrp
  
  res <- results(dds, contrast = c("conditions", experGrp, controlGrp))
  coefVar <- paste0("conditions_", experGrp, "_vs_", controlGrp)
  resLFC <- tryCatch(lfcShrink(dds, coef=coefVar, type="apeglm"), 
                     error = function(e){
                       cat("Using DESeq2 default method to calculate LFC\n")
                       return(lfcShrink(dds, coef=coefVar, type="normal"))
                     })
  ## DEA output
  deaResult <- do.call(cbind, resLFC@listData) %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene")
  write.csv(deaResult, paste0(fileLoc, "DEgeneStatistics_", 
                              experGrp, "_vs_", controlGrp, ".csv"), 
            row.names = F)
  return(list(
    dataProc = dataProc, resLFC = deaResult, experGrp = experGrp
  ))
}

DEAnalysis_limma <- function(dataProc = dataProc_RNA, experGrp,
                             fileLoc){
  sampleList <- dataProc$dataSet$sampleInfo
  eset <- dataProc$eset
  controlGrp <- dataProc$dataSet$controlGrp
  
  batremvOpt <- ifelse("batch" %in% colnames(sampleList), T, F)
  f <- factor(sampleList$conditions)
  designFormula <- eval(parse(text = ifelse(batremvOpt, "~ batch + f", "~ 0+f")))
  design <- model.matrix(designFormula)
  colnames(design) <- levels(f)
  
  fit <- lmFit(eset, design)
  constrast.matrix <- makeContrasts(contrasts=paste0(experGrp, '-', controlGrp),
                                    levels = design)
  fit2 <- contrasts.fit(fit, constrast.matrix)
  fit2 <- eBayes(fit2)
  resLFC <- topTable(fit2, adjust="BH", number = Inf)
  
  ## DEA output
  deaResult <- resLFC[order(resLFC$adj.P.Val), ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene")
  write.csv(deaResult, paste0(fileLoc, "DEgeneStatistics_", 
                              experGrp, "_vs_", controlGrp, ".csv"), 
            row.names = F)
  return(list(
    dataProc = dataProc, resLFC = deaResult, experGrp = experGrp, fit = fit2
  ))
}
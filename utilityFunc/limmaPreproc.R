limmaPreproc <- function(dataSet = dataSet_RNA, norm = T, 
                         fileLoc){
  eset <- dataSet$data
  sampleList <- dataSet$sampleInfo
  
  if(norm){
    par(mfrow=c(1,2))
    pdf(file = paste0(fileLoc, "/normalization_boxplot.pdf"))
    p <- boxplot(eset,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
    title(main = list('Before normalization',cex = 2 ,font = 2),
          xlab = list('Sample list',cex = 1.5,font = 2),
          ylab = '',line = 0.7)
    mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
    eset <- normalizeBetweenArrays(eset) 
    p1 <- boxplot(eset,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F)
    title(main = list('Normalization',cex = 2 ,font = 2),
          xlab = list('Sample list',cex = 1.5,font = 2),
          ylab = '',line = 0.7)
    mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
    dev.off()
  }else{
    pdf(file = paste0(fileLoc, "/expression_boxplot.pdf"))
    p <- boxplot(eset,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
    title(main = list('Expression Value',cex = 2 ,font = 2),
          xlab = list('Sample list',cex = 1.5,font = 2),
          ylab = '',line = 0.7)
    mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
    dev.off()
  }
  batremvOpt <- ifelse("batch" %in% colnames(sampleList), T, F)
  if(batremvOpt){
    batch <- sampleList$batch
    eset <- removeBatchEffect(eset, batch)
  }
  return(list(
    eset = eset, dataSet = dataSet
  ))
}
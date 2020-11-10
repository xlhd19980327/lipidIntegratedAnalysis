### This procedure only suitable for: 
### Single-Channel Microarray data, 
### Input file:
### unnormalized/suitably pre-processed normalized Expression List file
### (may use upstream anlysis package to get this file, not support R object input yet)
### 

library(limma)
###!!!Client options: control group, default will use group of the first sample in sample-list file
controlGrp <- "HCC_FASNpos"
#!!!Client options: Client should choose which group to do the comparison(DEA) or do group by group comparison("group_by_group")
experGrp <- "HCC_FASNneg"
#!!!Client options: Top DE genes number
topnum <- 100

data <- read.csv("./branch/benchmark/input/gene_tidy.CSV")
firstline <- scan("./branch/benchmark/input/gene_tidy.CSV", 
                  what = "character", nlines = 1, sep = ",", quote = "\"", na.strings = c("N/A", "NA"))
sampleList <- read.csv("./branch/benchmark/input/sampleList.CSV")
rownames(data) <- data[, 1]
colnames(data) <- firstline
eset <- subset(data, select = sampleList$samples)

## Normalization
pdf(file = paste0(datadir, gse, "/raw_box.pdf"))
p <- boxplot(eset,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
title(main = list('Before normalization',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()
eset <- normalizeBetweenArrays(eset) 
pdf(file = paste0(datadir, gse, "/normalized_box.pdf"))
p1 <- boxplot(eset,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F)
title(main = list('Normalization',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

batremvOpt <- ifelse("batch" %in% colnames(sampleList), T, F)
if(batremvOpt){
  batch <- sampleList$batch
  eset <- removeBatchEffect(eset, batch)
}
f <- factor(sampleList$conditions)
designFormula <- eval(parse(text = ifelse(batremvOpt, "~ batch + f", "~ 0+f")))
design <- model.matrix(designFormula)
colnames(design) <- levels(f)

fit <- lmFit(eset, design)
constrast.matrix <- makeContrasts(eval(parse(text = paste0(experGrp, '-', controlGrp))),
                                  levels = design)
fit2 <- contrasts.fit(fit, constrast.matrix)
fit2 <- eBayes(fit2)
data_top <- topTable(fit2, adjust="BH", number = topnum)



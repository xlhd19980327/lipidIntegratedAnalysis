library(MetaboAnalystR)
library(tidyverse)
##After MARpreproc, get dataSet, mSet

option_list <- list( 
  make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-p", "--metreg_output"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))

controlGrp <- dataSet$controlGrp
groupsLevel <- dataSet$groupsLevel
data_norm <- mSet$dataSet$norm
data_norm2 <- mSet$dataSet$row.norm
## Volcano analysis & plot
paired <- FALSE
equal.var <- TRUE
#!!!Client options: "raw" or "fdr"
pval.type <- "raw"
cprind <- mSet[["dataSet"]][["cls"]]  == controlGrp
exper <- groupsLevel[groupsLevel != controlGrp]
expind <- mSet[["dataSet"]][["cls"]]  == exper
p.value <- apply(as.matrix(data_norm), 2, 
                 function(x) {
                   tmp <- tryCatch(t.test(x[cprind], x[expind], paired = paired, 
                                          var.equal = equal.var), 
                                   error = function(e){cat("Note: A NA generate!\n")}
                   )
                   if (is.null(tmp)) {
                     return(NA)
                   }
                   else {
                     return(tmp$p.value)
                   }
                 })
if(pval.type == "fdr"){
  p.value <- p.adjust(p.value, "fdr")
}
p.log <- -log10(p.value)
cprmean <- colMeans(data_norm2[cprind, ])
expmean <- colMeans(data_norm2[expind, ])
ratio <- expmean/cprmean
fc.log <- signif(log2(ratio), 5)
result <- cbind(p.log, fc.log)
#!!!Client options: Fold change threshold
fcthresh <- 1.2
#!!!Client options: Fold change p.value
pthresh <- 0.1
target.data.up <- result %>%
  as.data.frame() %>%
  filter(
    fc.log > log(fcthresh, 2) & p.log > -log10(pthresh)
  )
target.data.down <- result %>%
  as.data.frame() %>%
  filter(
    fc.log < -log(fcthresh, 2) & p.log > -log10(pthresh)
  )
#background.data <- t(data_norm)

#Get met names
cmpd.vec_up<-rownames(target.data.up)
write.csv(cmpd.vec_up, paste0(opt$metreg_output, "up.csv"), row.names = F)
cmpd.vec_down<-rownames(target.data.down)
write.csv(cmpd.vec_down, paste0(opt$metreg_output, "down.csv"), row.names = F)

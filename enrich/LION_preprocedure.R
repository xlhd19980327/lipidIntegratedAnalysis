##After MARpreproc, get dataSet, mSet

library(dplyr)
library(tibble)
source("./utilityFunc/tidyLip_forLION.R")
dataType <- dataSet[["dataType"]]

#'Target-list mode'
controlGrp <- dataSet$controlGrp
groupsLevel <- dataSet$groupsLevel
data_norm <- mSet$dataSet$norm
#data_norm2 <- mSet$dataSet$row.norm
data_norm2 <- qs::qread("row_norm.qs")
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
fcthresh <- 2.0
#!!!Client options: Fold change p.value
pthresh <- 0.1
target.data.up <- result %>%
  as.data.frame() %>%
  #rownames_to_column(var = "lipid") %>%
  filter(
    fc.log > log(fcthresh, 2) & p.log > -log10(pthresh)
  )
target.data.down <- result %>%
  as.data.frame() %>%
  #rownames_to_column(var = "lipid") %>%
  filter(
    fc.log < -log(fcthresh, 2) & p.log > -log10(pthresh)
  )
background.data <- t(data_norm)
# rownames(target.data.up) <- tidy_lipid(rownames(target.data.up))
# rownames(target.data.down) <- tidy_lipid(rownames(target.data.down))
# rownames(background.data) <- tidy_lipid(rownames(background.data))
# write.csv(target.data.up, "~/temp/tar_data_up.csv")
# write.csv(target.data.down, "~/temp/tar_data_down.csv")
# write.csv(background.data, "~/temp/bg_data.csv")
lipid_target_list_up <- tidy_lipid(rownames(target.data.up), dataType = dataType)
lipid_target_list_down <- tidy_lipid(rownames(target.data.down), dataType = dataType)
lipid_background_list <- tidy_lipid(rownames(background.data), dataType = dataType)
write.csv(lipid_target_list_up, paste0(opt$lipreg_output, "up.csv"), row.names = F)
write.csv(lipid_target_list_down, paste0(opt$lipreg_output, "down.csv"), row.names = F)
write.csv(lipid_background_list, paste0(opt$lipreg_output, "background.csv"), row.names = F)

#'Ranking mode'
#"p_value"/"log2FC"
rankingVar <- "log2FC"
source("./utilityFunc/lipVolcanoPlot.R")
data_input <- lipVolcanoPlot(dataSet = dataSet, mSet = mSet, stat = T)
data_input$lipid <- tidy_lipid(lipids = data_input$lipid, dataType = dataType)
rankarg <- switch (rankingVar,
  log2FC = "fc.log", 
  p_value = "p.value"
)

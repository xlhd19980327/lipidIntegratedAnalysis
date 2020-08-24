### For LipidSearch data ###

## Source will offer the following variables:
## allgroups, groupsLevel, mSet
source("D:/myLearning/lipGroup/riverGroup/integerateOmics/lipidPathways/newInspirationWork/procedure/github/lipidpreproc/LipidSearchData_preprocess.R")

## Data filtering
# The purpose of the data filtering is to identify and remove variables that are unlikely to be of use when modeling the data. 
# Default methods:
# None for targeted lipidomics and most other lipidomics data
mSet<-FilterVariable(mSet, 
                     filter = "none", qcFilter = "F", rsd = 25)
## Normalization / Scaling
# Default methods: 
# 1. Normalization method: Probabilistic Quotient Normalization(PQN) without using a reference sample
# 2. Transformation method: None
# 3. Scaling method: Auto scaling(mean-centered and divided by the standard deviation of each variable)
mSet<-PreparePrenormData(mSet)
#!!!Client options
mSet<-Normalization(mSet, 
                    rowNorm = "ProbNormT", transNorm = "NULL", scaleNorm = "AutoNorm", ref = NULL, 
                    ratio=FALSE, ratioNum=20)


data_norm <- mSet$dataSet$norm
data_norm2 <- mSet$dataSet$row.norm
data_prenorm <- mSet$dataSet$prenorm

### Statistics and visualization with MAR ###
## Destination code
# mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, 0.75,F, 0.1, TRUE, "raw")
# mSet<-PlotVolcano(mSet, "volcano_0_",1, "png", 72, width=NA)
# mSet<-PCA.Anal(mSet)
# mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,0,0)
# mSet<-PlotHeatMap(mSet, "heatmap_3_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "detail", F, T, NA, T, F)
# mSet<-PlotSubHeatMap(mSet, "heatmap_2_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "tanova", 75, "detail", F, T, T, F)
## Volcano analysis & plot
paired <- FALSE
equal.var <- TRUE
pval.type <- "raw"
cprind <- mSet[["dataSet"]][["cls"]]  == controlGrp
#!!!!!Control flow WARNING: group_by_group analysis will not use this control flow
#!!!!!Control flow WARNING: Only all_together analysis will use this control flow
for(i in groupsLevel[groupsLevel != controlGrp]){
  expind <- mSet[["dataSet"]][["cls"]]  == i
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
}

fcthresh <- 2.0
pthresh <- 0.1
volcano.data <- result %>%
  as.data.frame() %>%
  rownames_to_column(var = "lipid") %>%
  mutate(
    up = ifelse(fc.log > log(fcthresh, 2) & p.log > -log10(pthresh), 1, 0),
    down = ifelse(fc.log < -log(fcthresh, 2) & p.log > -log10(pthresh), 2, 0),
    regState = sapply(up + down, function(x) switch(x, "upreg", "downreg")), 
    regState = sapply(regState, function(x) ifelse(is.null(x), 0, x)), 
    Class = gsub("(.*?)\\(.*", "\\1", lipid)
  )
volcano.data_reg <- subset(volcano.data, subset = regState != 0)
volcano.data_unreg <- subset(volcano.data, subset = regState == 0)
volcano.plot <- ggplot() +
  geom_point(data = subset(volcano.data, subset = regState != 0), 
             aes(x = fc.log, y = p.log, color = factor(regState)))+ 
  geom_point(data = subset(volcano.data, subset = regState == 0), 
             aes(x = fc.log, y = p.log), color = "gray") +
  theme_bw() +
  geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
  labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "regulation", 
       title = paste0(i, " vs ", controlGrp)) +
  theme(plot.title = element_text(hjust = 0.5, size = 10))+ 
  scale_color_aaas()
volcano.plot
volcano.plot_class <- ggplot() +
  geom_point(data = subset(volcano.data, subset = regState != 0), 
             aes(x = fc.log, y = p.log, color = reorder(Class, -abs(fc.log), mean)))+ 
  geom_point(data = subset(volcano.data, subset = regState == 0), 
             aes(x = fc.log, y = p.log), color = "gray") +
  theme_bw() +
  geom_vline(xintercept = c(-log(fcthresh, 2), log(fcthresh, 2)), linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = -log10(pthresh), linetype = "dashed", size = 0.5) + 
  labs(x = "log2(Fold change)", y = "-log10(p.value)", color = "Class", 
       title = paste0(i, " vs ", controlGrp)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))#+ scale_color_aaas()
volcano.plot_class
## PCA & plot
pcaAnal <- prcomp(data_norm, center = TRUE, scale = F)
pclevel <- c(1, 2)
pcaScore <- as.data.frame(pcaAnal$x)[, pclevel] %>%
  rownames_to_column(var = "sampleLabel") %>%
  mutate(group = mSet$dataSet$cls)
pcaVar <- summary(pcaAnal)$importance[2, ]
ellipse_coor_all <- data.frame()
for(i in groupsLevel){
  groupVar <- var(subset(pcaScore, subset = group == i)[, 2:3], na.rm = T)
  groupMean <- cbind(mean(subset(pcaScore, subset = group == i)[, 2], na.rm = T), 
                     mean(subset(pcaScore, subset = group == i)[, 3], na.rm = T))
  ellipse_coor <- as.data.frame(ellipse::ellipse(groupVar, centre = groupMean, 
                                                 level = 0.95, npoints = 100))
  colnames(ellipse_coor) <- c("x", "y")
  ellipse_coor_all <- rbind(ellipse_coor_all, 
                            cbind(ellipse_coor, group = i))
}
pcaScorePlot <- ggplot() +
  geom_polygon(ellipse_coor_all, mapping = aes(x = x, y = y, color = group, fill = group)) +
  geom_point(pcaScore, mapping =aes(x = PC1, y = PC2, color = group)) +
  scale_color_aaas() +
  scale_fill_aaas()
colorpars_all <- ggplot_build(pcaScorePlot)$data[[2]]$colour
names(colorpars_all) <- pcaScore$group
colorpars <- unique(colorpars_all)
names(colorpars) <- unique(pcaScore$group)
pcaScorePlot <- pcaScorePlot + 
  scale_fill_manual(values = alpha(colorpars, 0.3)) +
  theme_bw() + 
  labs(x = paste0("PC", pclevel[1], "(", round(100 * pcaVar[pclevel[1]], 1), "%)"), 
       y = paste0("PC", pclevel[2], "(", round(100 * pcaVar[pclevel[2]], 1), "%)"), 
       color = "group", 
       title = "PCA Score Plot") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
ggsave("PCA_score_plot_all.pdf", plot = pcaScorePlot, 
       device = "pdf", width = 9, height = 9)
## Heatmap
#!!!Client options
plotHeatmap <- function(mSetObj = mSet, 
                        norm = T, topind = NA){
  if(norm == T){
    heatmapdata <- data_norm
  }else{
    heatmapdata <- data_prenorm
  }
  if(!all(is.na(topind))){
    heatmapdata <- heatmapdata[topind]
  }
  datagroup <- mSetObj$dataSet$cls
  x <-pheatmap::pheatmap(mat = t(heatmapdata), 
                         annotation = data.frame(group = datagroup, row.names = rownames(heatmapdata)), 
                         fontsize = 8, 
                         fontsize_row = 8, 
                         clustering_distance_rows = "euclidean", 
                         clustering_distance_cols = "euclidean", 
                         clustering_methods = "ward.D", 
                         clustering_rows = T, 
                         clustering_cols = F, 
                         scale = "row", 
                         annotation_colors = list(group = colorpars)
  )
  pdf("test.pdf", width=(nrow(heatmapdata)*25+300)/72, height=(ncol(heatmapdata)*18+150)/72)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#!!!Client options
plotSubHeatmap <- function(mSetObj = mSet, 
                           norm = T, topnum = 75){
  if(length(groupsLevel) == 2){
    mSetObj <- Ttests.Anal(mSetObj)
    toplip <- names(sort(mSetObj$analSet$tt$p.value))[1:topnum]
    toplip_ind <- match(toplip, colnames(data_norm))
  }else if(length(groupsLevel) > 2){
    mSetObj <- ANOVA.Anal(mSetObj)
    toplip <- names(sort(mSetObj$analSet$aov$p.value))[1:topnum]
    toplip_ind <- match(toplip, colnames(data_norm))
  }
  plotHeatmap(mSet, topind = toplip_ind)
}
plotHeatmap(mSet)
plotSubHeatmap(mSet)



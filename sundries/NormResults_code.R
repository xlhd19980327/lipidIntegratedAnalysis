### Normalization Assessment Methods:
### D:\myLearning\lipGroup\riverGroup\Literature\integratedOmicsModeling\normalization
### Paper: https://kns.cnki.net/KCMS/detail/detail.aspx?dbname=CMFD201901&filename=1018853644.nh
### page: 41

## Orignial code for "Samples" view about normalization
## Function: PlotSampleNormSummary
#Before normalization
boxplot(t(mSetObj$dataSet$proc[pre.inx, , drop = FALSE]), 
        names = namesVec, ylim = rangex.pre, las = 2, col = "lightgreen", 
        horizontal = T)
plot(density(apply(mSetObj$dataSet$proc, 1, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
#After normalization
boxplot(t(mSetObj$dataSet$norm[norm.inx, , drop = FALSE]), 
        names = namesVec, ylim = rangex.norm, las = 2, col = "lightgreen", 
        ylab = "", horizontal = T)
plot(density(apply(mSetObj$dataSet$norm, 1, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")

## Orignial code for "Features" view about normalization
## Function: PlotNormSummary
#Before normalization
boxplot(mSetObj$dataSet$proc[, pre.inx, drop = FALSE], names = namesVec, 
        ylim = rangex.pre, las = 2, col = "lightgreen", horizontal = T, 
        show.names = T)
plot(density(apply(mSetObj$dataSet$proc, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")
#After normalization
boxplot(mSetObj$dataSet$norm[, norm.inx, drop = FALSE], 
        names = namesVec, ylim = rangex.norm, las = 2, col = "lightgreen", 
        horizontal = T, show.names = T)
plot(density(apply(mSetObj$dataSet$norm, 2, mean, na.rm = TRUE)), 
     col = "darkblue", las = 2, lwd = 2, main = "", xlab = "", 
     ylab = "")

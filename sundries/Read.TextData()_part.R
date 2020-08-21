mSet$dataSet$cls.type <- "disc"
  mSet$dataSet$format <- "colu"
dat <- read.csv()
msg <- c(msg, "Samples are in columns and features in rows.")
      var.nms <- dat[-1, 1]
      dat[, 1] <- NULL
      smpl.nms <- colnames(dat)
      cls.lbl <- dat[1, ]
      conc <- t(dat[-1, ])
mSet$dataSet$type.cls.lbl <- class(cls.lbl)
empty.inx <- is.na(smpl.nms) | smpl.nms == ""
empty.inx <- is.na(cls.lbl) | cls.lbl == ""

#smpl.nms <- CleanNames(smpl.nms, "sample_name")
  orig.var.nms <- var.nms
#  var.nms <- CleanNames(var.nms, "var_name")
  names(orig.var.nms) <- var.nms
#  cls.lbl <- ClearStrings(as.vector(cls.lbl))
  rownames(conc) <- smpl.nms
  colnames(conc) <- var.nms
  
      mSet$dataSet$orig.cls <- mSet$dataSet$cls <- as.factor(as.character(cls.lbl))

mSet$dataSet$cmpd <- var.nms
mSet$dataSet$mumType <- "table"
  mSet$dataSet$orig.var.nms <- orig.var.nms
  mSet$dataSet$orig <- conc
  mSet$msgSet$read.msg <- c(msg, paste("The uploaded data file contains ", 
    nrow(conc), " (samples) by ", ncol(conc), " (", tolower(GetVariableLabel(mSet$dataSet$type)), 
    ") data matrix.", sep = ""))
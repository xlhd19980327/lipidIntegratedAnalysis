##!!!!!WARNING:: this p_value calculation may be not appropriate?
#Totally equal to t.test() one by one

getPValue <- function(x, calcType, controlGrpin = controlGrp){
  #Statistics in different situations
  if(calcType == "headgroup"){
    alls <- unique(x$Class)
    name <- "Class"
  }
  if(calcType == "FAchain"){
    alls <- unique(x$subclass)
    name <- "subclass"
  }
  #Totally equal to t.test() one by one
  p <- pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value
  p1 <- p[rownames(p) == controlGrpin, ]
  if(length(p1) != 0){
    names(p1) <- colnames(p)
  }
  p2 <- p[, colnames(p) == controlGrpin]
  if(length(p2) != 0){
    names(p2) <- rownames(p)
  }
  p_tidy <- c(p1, p2)
  p_tidy <- p_tidy[!is.na(p_tidy)]
  if(length(p_tidy) == 0){
    result <- NULL
  }else {
    result <- data.frame(
      group = c(names(p_tidy), controlGrpin),
      p = c(p_tidy, NA), 
      alls
    )
    colnames(result)[3] <- name
  }
  return(result)
}

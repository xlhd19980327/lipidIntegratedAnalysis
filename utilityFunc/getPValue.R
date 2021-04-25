##!!!!!WARNING:: this p_value calculation may be not appropriate?
#Totally equal to t.test() one by one

getPValue <- function(x, calcType, controlGrpin){
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
  p <- tryCatch(pairwise.t.test(x$lipidsum, x$group, p.adjust.method = "none", pool.sd = F)$p.value, 
                error = function(e){
                  return(NULL)
                })
  if(is.null(p)){
    allgroups <- unique(x$group)
    sing_t <- c()
    for(i in allgroups[allgroups != controlGrpin]){
      tempt <- tryCatch(t.test(x$lipidsum[allgroups == controlGrpin], 
                               x$lipidsum[allgroups == i])$p.value, 
                        error = function(e){
                          return(NA)
                        })
      sing_t <- c(sing_t, tempt)
    }
    result <- data.frame(group = c(allgroups[allgroups != controlGrpin], controlGrpin), 
                         p = c(sing_t, NA), 
                         V3 = alls)
    colnames(result)[3] <- name
  }else{
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
  }
  return(result)
}

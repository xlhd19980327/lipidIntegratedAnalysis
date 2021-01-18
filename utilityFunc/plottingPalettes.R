plottingPalettes <- function(n, type){
  if(type == "continuous"){
    #retired, cuz ugly
    color <- pal_gsea("default", n = n, reverse = F)(n)
  }
  if(type == "discrete"){
    if(n <= 10){
      color <- pal_npg("nrc")(n)
    }else if(n <= 20){
      color <- c(pal_npg("nrc")(10), pal_cosmic("hallmarks_light")(n-10))
    }else{
      #library("scales")
      #color <- c(pal_npg("nrc")(10), pal_cosmic("hallmarks_light")(n-10))
      
      color <- c(pal_npg("nrc")(10), pal_cosmic("hallmarks_light")(10), 
                 colorRampPalette(brewer.pal(9, "Set1"))(n-20))
      # ind <- which(nchar(color) == 9)
      # color[ind] <- gsub("FF$", "", color[ind])
    }
  }
  return(color)
}

plottingPalettes <- function(n, type){
  if(type == "continuous"){
    color <- pal_gsea("default", n = n, reverse = F)(n)
  }
  if(type == "discrete"){
    if(n <= 10){
      color <- pal_aaas("default")(n)
    }else{
      color <- colorRampPalette(brewer.pal(9, "Set1"))(n)
    }
  }
  return(color)
}

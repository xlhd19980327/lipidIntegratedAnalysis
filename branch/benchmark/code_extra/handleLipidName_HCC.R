handleLipidName_HCC <- function(x){
  if(grepl(".*?\\(.*\\)", x)){
    Class <- gsub("(.*?)\\(.*", "\\1", x)
    otherinfo <- gsub(".*?(\\(.*)", "\\1", x)
    otherinfo <- gsub("/", "_", otherinfo)
    if(Class == "CE"){
      Class <- "ChE"
    }
    if(Class == "CER"){
      Class <- "Cer"
    }
    if(Class == "DAG"){
      Class <- "DG"
    }
    if(Class == "FFA"){
      Class <- "FA"
    }
    if(Class == "CE"){
      Class <- "ChE"
    }
    combineName <- paste0(Class, otherinfo)
  }
  if(grepl("\\-FA", x)){
    if(grepl("TAG", x)){
      Class <- "TG"
      tgall <- gsub("TAG(.*)\\-FA.*", "\\1", x)
      fa <- gsub(".*FA(.*)", "\\1", x)
      combineName <- paste0("TG", "(", tgall, "_", fa, ")")
    }
  }
  return(combineName)
}

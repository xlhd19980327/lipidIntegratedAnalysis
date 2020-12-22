tidy_lipid <- function(lipids, dataType){
  if(dataType == "LipidSearch"){
    lipids <- lipids
  }
  if(dataType == "MS_DIAL"){
    lipids <- gsub("^.*\\|", "", lipids)
  }
  return(lipids)
}
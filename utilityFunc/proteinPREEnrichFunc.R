proteinPREEnrichFunc <- function(protein){
  ##See https://www.uniprot.org/help/api_idmapping
  library(httr)
  myprotein <- paste(protein, collapse = " ")
  print(myprotein)
  dataList <- list(
    from = 'ACC+ID',
    to = 'P_ENTREZGENEID',
    format = 'tab',
    query = myprotein
  )
  datapost <- POST('https://www.uniprot.org/uploadlists/', 
                   body = dataList)
  dataraw <- httr::content(datapost, as = "text", encoding = "UTF-8")
  data <- readr::read_delim(dataraw, delim = '\t')
  print(data$To)
  return(data$To)
}
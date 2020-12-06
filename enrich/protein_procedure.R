library(httr)
dataList <- list(
  from = 'ACC+ID',
  to = 'ENSEMBL_ID',
  format = 'tab',
  query = 'P40925 P40926 O43175 Q9UM73 P97793'
)
datapost <- POST('https://www.uniprot.org/uploadlists/', 
                 body = dataList)
dataraw <- content(datapost, as = "text", encoding = "UTF-8")
data <- read_delim(dataraw, delim = '\t')
data$To

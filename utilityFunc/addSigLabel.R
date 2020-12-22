addSigLabel <- function(x){
  x[is.na(x)] <- 2
  label <- c('ns', "*", "**", "***", "****")
  allLable <- ifelse(x == 2, '', NA)
  allLable <- ifelse(x > 0.05 & x != 2, label[1], allLable)
  allLable <- ifelse(x <= 0.05 & x > 0.01, label[2], allLable)
  allLable <- ifelse(x <= 0.01 & x > 0.001, label[3], allLable)
  allLable <- ifelse(x <= 0.001 & x > 0.0001, label[4], allLable)
  allLable <- ifelse(x <= 0.0001, label[5], allLable)
}

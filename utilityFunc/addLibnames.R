addLibnames <- function(x, allcolnames){
  libInfo <- read.csv("./patternHunting/lipidReatomeStat/headgroupsvslibnames.csv")
  result <- c()
  names(x) <- allcolnames
  ind <- (libInfo$Headgroup == x[["Class"]] & libInfo$otherInfo == "") | libInfo$otherInfo == x[["subclass"]]
  result <- c(result, libInfo$Libnames[ind])
  ind_amb <- grep("xx", libInfo$otherInfo)
  info_amb <- libInfo$otherInfo[ind_amb]
  ##!!!!!WARNING: SM may cause something wrong, simply handle
  if(x[["Class"]] == "SM" & !grepl("\\([a-z][0-9]+:", x[["subclass"]])){
    result <- c(result, "SM(d18:1/xx:xx)")
  }
  ##Now the ambiguous module only contain FA, simply use following procedure
  getambInfo <- function(y){
    chains_l <- gsub(".*?([0-9]*)_([0-9]*).*", "\\1", y)
    chains_r <- gsub(".*?([0-9]*)_([0-9]*).*", "\\2", y)
    chains <- c(as.integer(chains_l), as.integer(chains_r))
    return(chains)
  }
  ambInfo <- lapply(info_amb, getambInfo)
  testChainLength <- function(z, input){
    chain <- as.integer(gsub(".*?([0-9]+):.*", "\\1", input))
    z[1] <- ifelse(is.na(z[1]), 1, z[1])
    z[2] <- ifelse(is.na(z[2]), 9999999, z[2])
    res <- ifelse(z[1] <= chain & z[2] >= chain, T, F)
    return(res)
  }
  if(x[["Class"]] == "FA"){
    ind_fa <- sapply(ambInfo, testChainLength, x[["subclass"]])
    ind_all <- ind_amb[ind_fa]
    result <- c(result, libInfo$Libnames[ind_all])
  }
  result <- paste(result, collapse = ";")
  return(result)
}

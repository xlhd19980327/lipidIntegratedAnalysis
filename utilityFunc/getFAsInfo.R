##!!!!!WARNING: Class statistics will only contain the following lipid class:
## Cer, SM, SPH, FA, MG, DG, TG, PA, PC, PE, PG, PI, PS, LPA, LPC, LPE, LPG, LPI, LPS, ChE, CL

getFAsInfo <- function(i){
  Class <- gsub("(.*?)\\(.*", "\\1", i, perl = T)
  if(Class %in% c("Cer", "SM", "SPH")){ 
    # Ignore "+O" info
    if(grepl("_", i)){
      spbase <- gsub(".*\\(([mdt][0-9]+:[0-9]+)_.*", "\\1", i, perl = T)
      fa <- gsub(".*_([0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- c(spbase, fa)
    } else {
      spall <- gsub(".*\\(([mdt][0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- spall
    }
  } else if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                         "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                         "LPS", "ChE", "CL")){
    # All ignore e/p(O-/P-) connection
    if(grepl("_", i)){
      fainfo <- gsub(".*\\((.*)\\)", "\\1", i)
      fas <- strsplit(fainfo, "[_/]")
      m <- gregexpr("[0-9]+:[0-9]+", fas)
      fas <- unlist(regmatches(as.character(fas), m))
    } else {
      m <- gregexpr("[0-9]+:[0-9]+", i)
      fas <- unlist(regmatches(i, m))
      ## ChE may have some 0:0 info 
      if(length(fas) == 0){
        fas <- "0:0"
      }
    }
  }else{
    result <- list(otherClass = "UnknownPattern", ms1 = NA)
    names(result) <- c(Class, "ms1")
    return(result)
  }
  n <- length(fas)
  if(n > 1){
    ms1 <- F
  } else{
    if(Class %in% c("Cer", "DG", "PA", "PC", "PE", "PG", "PI", "PS", "SM", "TG", "CL")){
      ms1 <- T
    } else {
      ms1 <- F
    }
  }
  result <- list(fas, ms1)
  names(result) <- c(Class, "ms1")
  return(result)
}
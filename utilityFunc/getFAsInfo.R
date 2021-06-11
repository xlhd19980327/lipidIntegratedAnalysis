##!!!!!WARNING: Class statistics will only contain the following lipid class:
## Cer, SM, SPH, FA, MG, DG, TG, PA, PC, PE, PG, PI, PS, LPA, LPC, LPE, LPG, LPI, LPS, ChE, CL
## & its O/P Class
## New: "Hex1Cer", "Hex2Cer", "AcCa", "MLCL" 

getFAsInfo <- function(i, ignorei){
  ## Source will offer the following contents:
  ## Function(s): getClassInfo
  source("./utilityFunc/getClassInfo.R")
  
  Class <- getClassInfo(i, "LipidSearch", ignore = ignorei)
  spclasses <- c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                 "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                 "LPS", "ChE", "CL", "AcCa", "MLCL")
  epclasses <- c(paste0(spclasses, "(O)"), paste0(spclasses, "(P)"))
  spclasses2 <- c("DG", "PA", "PC", "PE", "PG", "PI", "PS", "TG", "CL", "MLCL")
  epclasses2 <- c(paste0(spclasses2, "(O)"), paste0(spclasses2, "(P)"))
  if(Class %in% c("Cer", "SM", "SPH", "Hex1Cer", "Hex2Cer")){ 
    # Ignore "+O" info
    if(grepl("_", i)){
      spbase <- gsub(".*\\(([mdt][0-9]+:[0-9]+)_.*", "\\1", i, perl = T)
      fa <- gsub(".*_([0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- c(spbase, fa)
    } else {
      spall <- gsub(".*\\(([mdt][0-9]+:[0-9]+)(\\+O)*.*", "\\1", i, perl = T)
      fas <- spall
    }
  } else if(Class %in% c(spclasses, epclasses)){
    # All ignore e/p connection
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
    if(Class %in% c("Cer", "SM", spclasses2, epclasses2)){
      ms1 <- T
    } else {
      ms1 <- F
    }
  }
  result <- list(fas, ms1)
  names(result) <- c(Class, "ms1")
  return(result)
}

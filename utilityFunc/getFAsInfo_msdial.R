##!!!!!WARNING: Class statistics will only contain the following lipid class:
## "Cer", "SM", "FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS"

getFAsInfo_msdial <- function(i){
  Class <- gsub("(.*?) .*", "\\1", i, perl = T)
  if(Class %in% c("Cer", "SM")){ #may use Class %in% c("Cer", "SM", "SPH")
    if(grepl("\\|", i)){
      spbase <- gsub(".*\\|.*?([0-9]+:[0-9]+;[0-9]O).*", "\\1", i, perl = T)
      spbase_OH <- gsub("[0-9]+:[0-9]+;([0-9])O", "\\1", spbase, perl = T)
      spbase_OH_label <- switch (spbase_OH,
                                 '1' = "m",
                                 '2' = "d", 
                                 '3' = "t"
      )
      spbase_carbon <- gsub("([0-9]+:[0-9]+);[0-9]O", "\\1", spbase, perl = T)
      spbase_tidy <- paste0(spbase_OH_label, spbase_carbon)
      fa <- gsub(".*/([0-9]+:[0-9]+).*", "\\1", i, perl = T)
      fas <- c(spbase_tidy, fa)
    } else {
      spall <- gsub(".*?([0-9]+:[0-9]+;[0-9]O).*", "\\1", i, perl = T)
      spall_OH <- gsub("[0-9]+:[0-9]+;([0-9])O", "\\1", spall, perl = T)
      spall_OH_label <- switch (spall_OH,
                                '1' = "m",
                                '2' = "d", 
                                '3' = "t"
      )
      spall_carbon <- gsub("([0-9]+:[0-9]+);[0-9]O", "\\1", spall, perl = T)
      spall_tidy <- paste0(spall_OH_label, spall_carbon)
      fas <- spall_tidy
    }
  } else if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                         "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                         "LPS")){#will not contain "CL" or "ChE"
    # All ignore e/p(O-/P-) connection
    if(grepl("\\|", i)){
      fainfo <- gsub(".*\\|(.*)", "\\1", i)
      fas <- strsplit(fainfo, "[_/]")
      m <- gregexpr("[0-9]+:[0-9]+", fas)
      fas <- unlist(regmatches(as.character(fas), m))
    } else {
      m <- gregexpr("[0-9]+:[0-9]+", i)
      fas <- unlist(regmatches(i, m))
    }
  }else{
    result <- list(otherClass = Class, ms1 = NA)
    return(result)
  }
  n <- length(fas)
  if(n > 1){
    ms1 <- F
  } else{
    if(Class %in% c("Cer", "DG", "PA", "PC", "PE", "PG", "PI", "PS", "SM", "TG")){
      ms1 <- T
    } else {
      ms1 <- F
    }
  }
  result <- list(fas, ms1)
  names(result) <- c(Class, "ms1")
  return(result)
}

getClassInfo <- function(i, datatype, ignore = T){
  Class = switch(datatype,
                 LipidSearch = gsub("(.*?)\\(.*", "\\1", i), 
                 MS_DIAL = gsub("(.*?) .*", "\\1", i))
  if(!ignore){
    if(Class %in% c("FA", "MG", "DG", "TG", "PA", "PC", "PE", "PG", 
                    "PI", "PS", "LPA", "LPC", "LPE", "LPG", "LPI", 
                    "LPS", "ChE", "CL")){
      if(datatype == "LipidSearch"){
        if(grepl("[0-9]+:[0-9]+([ep])", i)){
          type <- gsub(".*?[0-9]+:[0-9]+([ep]).*", "\\1", i)
          type <- switch (type,
                          e = "(O)",
                          p = "(P)"
          )
          Class <- paste0(Class, type)
        }
      }
      if(datatype == "MS_DIAL"){
        if(grepl("[OP]\\-[0-9]+:[0-9]+", i)){
          type <- gsub(".*?([OP])\\-[0-9]+:[0-9]+.*", "\\1", i)
          type <- switch (type,
                          O = "(O)",
                          P = "(P)"
          )
          Class <- paste0(Class, type)
        }
      }
    }
  }
  if(nchar(Class) > 10){
    Class <- substr(Class, 1, 10)
  }
  return(Class)
}
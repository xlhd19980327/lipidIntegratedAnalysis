library(optparse)

option_list <- list( 
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  #make_option(c("-c", "--control_group"), action="store", default = ""), 
  make_option(c("-t", "--data_type"), action="store"), 
  #make_option(c("-f", "--feature_field"), action="store", default = NA, type = "character"), 
  make_option(c("-l", "--lipid_odd_chain_deletion"), action="store", default = F), 
  make_option(c("-n", "--NA_string"), action="store", default = NULL, type = "character"), 
  
  make_option(c("-p", "--output_temp"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

datafile <- opt$input_file
sampleList <- opt$description_file
#controlGrp = "", 
dataType <- opt$data_type
#lipField = NA, 
delOddChainOpt <- opt$lipid_odd_chain_deletion
na.char <- opt$NA_string

sampleInfo <- read.csv(sampleList)
allgroups <- sampleInfo$conditions
groupsLevel <- unique(allgroups)
write.csv(groupsLevel, paste0(opt$output_temp, "groupsLevel.csv"), row.names = F)
firstline <- scan(datafile, what = "character", nlines = 1, sep = ",", quote = "\"", 
                  na.strings = c("N/A", "NA"))
if(dataType %in% c("LipidSearch")){
  firstline <- firstline[which(firstline != "LipidIon")]
  firstline <- c("LipidIon", firstline)
}
if(dataType %in% c("MS_DIAL")){
  firstline <- firstline[which(firstline != "Metabolite name")]
  firstline <- c("Metabolite name", firstline)
}
write.csv(firstline, paste0(opt$output_temp, "firstline.csv"), row.names = F)
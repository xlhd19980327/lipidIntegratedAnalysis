library(optparse)

option_list <- list( 
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"), 
  make_option(c("-s", "--subgroups"), action="store"), 
  make_option(c("-p", "--output_files"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

datafile <- opt$input_file
sampleList <- opt$description_file
options(stringsAsFactors = F)
selGrps <- read.table(paste0(opt$subgroups, "subgroup.txt"))[, 1, drop = T]
data <- read.csv(datafile, na.strings = "", comment.char = "", check.names = F)
sampleInfo <- read.csv(sampleList, colClasses = "character")
descOut <- subset(sampleInfo, subset = conditions %in% selGrps)
dataOut <- subset(data, select = c(colnames(data)[1], descOut$samples))
write.csv(dataOut, paste0(opt$output_files, "input.csv"), row.names = F)
write.csv(descOut, paste0(opt$output_files, "description.csv"), row.names = F)

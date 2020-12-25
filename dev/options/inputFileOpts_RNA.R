library(optparse)

option_list <- list( 
  #make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-d", "--description_file"), action="store"),
  #make_option(c("-c", "--control_group"), action="store", default = ""),
  make_option(c("-p", "--output_temp"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sampleList <- opt$description_file
sampleInfo <- read.csv(sampleList)
groupsLevel <- unique(sampleInfo$conditions)
write.csv(groupsLevel, paste0(opt$output_temp, "groupsLevel_RNA.csv"), row.names = F)


library(optparse)
#source("./utilityFunc/proteinPREEnrichFunc.R")
source("./utilityFunc/geneEnrichFunc.R")

option_list <- list( 
  make_option(c("-i", "--input_file"), action="store"),
  make_option(c("-o", "--output_loc"), action="store"), 
  
  make_option(c("-t", "--species"), action="store", default = "mmu"), 
  #make_option(c("-g", "--gene_type"), action="store", default = "ENSEMBL"), 
  make_option(c("-s", "--show_num"), action="store", default = 20), 
  make_option(c("-c", "--go_term"), action="store", default = "Biological_Process")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

up <- read.csv(paste0(opt$input_file, "up.csv"))$x
down <- read.csv(paste0(opt$input_file, "down.csv"))$x
if(length(up) == 0){
  cat("No UP regulation proteins found! Try reduce the threshold of Fold change/p.value!\n")
}else{
  #data <- proteinPREEnrichFunc(up)
  geneEnrichFunc(genes = up, spe = opt$species, gene_type = "UNIPROT", 
                 shownum = opt$show_num, gocat = opt$go_term, reg = "up", loc = opt$output_loc)
}
if(length(down) == 0){
  cat("No DOWN regulation proteins found! Try reduce the threshold of Fold change/p.value!\n")
}else{
  #data <- proteinPREEnrichFunc(down)
  geneEnrichFunc(genes = down, spe = opt$species, gene_type = "UNIPROT", 
                 shownum = opt$show_num, gocat = opt$go_term, reg = "down", loc = opt$output_loc)
}

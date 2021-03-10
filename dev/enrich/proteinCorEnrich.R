library(optparse)
options(stringsAsFactors = F)
source("./utilityFunc/geneEnrichFunc.R")

option_list <- list( 
  #make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-i", "--data_file"), action="store"),
  #make_option(c("-j", "--row"), action="store", type = integer), 
  make_option(c("-k", "--colum"), action="store", type = "integer"), 
  
  make_option(c("-t", "--species"), action="store", default = "mmu"), 
  #make_option(c("-g", "--gene_type"), action="store", default = "SYMBOL"), 
  make_option(c("-s", "--show_num"), action="store", default = 20), 
  make_option(c("-c", "--go_term"), action="store", default = "Biological_Process"), 
  
  make_option(c("-o", "--output_loc"), action="store"),
  make_option(c("-p", "--output_temp"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dataloc <- paste0(opt$data_file, 
                  dir(opt$data_file, pattern = paste0("genes_", opt$colum, "\\.csv")))
#load(paste0(opt$rdata_file, "data.RData"))
#Input
genes <- read.csv(dataloc)$x
geneEnrichFunc(genes = genes, spe = opt$species, gene_type = "UNIPROT", 
               shownum = opt$show_num, gocat = opt$go_term, reg = "none", loc = opt$output_loc)
save(go, file = paste0(opt$output_temp, "data_circos.RData"))

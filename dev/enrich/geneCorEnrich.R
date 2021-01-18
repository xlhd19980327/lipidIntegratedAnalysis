library(tidyverse)
library(clusterProfiler)
library(cowplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(optparse)

option_list <- list( 
  #make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-i", "--data_file"), action="store"),
  #make_option(c("-j", "--row"), action="store", default = 4), 
  make_option(c("-k", "--colum"), action="store", default = 4), 
  
  make_option(c("-t", "--species"), action="store", default = "mmu"), 
  make_option(c("-g", "--gene_type"), action="store", default = "ENSEMBL"), 
  make_option(c("-s", "--show_num"), action="store", default = 50), 
  make_option(c("-c", "--go_term"), action="store", default = "Biological_Process"), 
  
  make_option(c("-o", "--output_loc"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dataloc <- paste0(opt$data_file, 
                  dir(opt$data_file, pattern = paste0("genes_", opt$colum, "\\.csv")))
#load(paste0(opt$rdata_file, "data.RData"))
#Input
gene_up <- read.csv(dataloc)$x
#!!!Client options: gene type(sugg: "ENSEMBL"/"SYMBOL")(ref: keytypes([spe db]))
gene_type <- opt$gene_type
# if(!data_type %in% c("RNAseq", "MiAr")){
#   library(httr)
#   dataList <- list(
#     from = 'ACC+ID',
#     to = 'ENSEMBL_ID',
#     format = 'tab',
#     query = paste(gene_up, collapse = " ")
#   )
#   datapost <- POST('https://www.uniprot.org/uploadlists/', 
#                    body = dataList)
#   dataraw <- content(datapost, as = "text", encoding = "UTF-8")
#   data <- read_delim(dataraw, delim = '\t')
#   gene_up <- data$To
#   gene_type <- "ENSEMBL"
# }
#!!!Client options: species("hsa"/"mmu")
spe <- opt$species
#!!!Client options: show ontology item number
shownum <- opt$show_num
#!!!Client options: select go term catalog("BP"/"CC"/"MF"/"ALL")
gocat <- opt$go_term
orgdb <- switch (spe,
                 hsa = org.Hs.eg.db, 
                 mmu = org.Mm.eg.db
)
goopt <- switch (gocat,
                 Biological_Process = "BP", 
                 Cellular_Component = "CC",
                 Molecular_Function = "MF", 
                 ALL = "ALL"
)
if(goopt == "BP"){
  go_BP <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  p <- barplot(go_BP,showCategory=shownum,drop=T) +
    ggtitle("Biological Process") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(opt$output_loc, "GOenrich_Biological_Process.pdf"), plot = p, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(goopt == "CC"){
  go_CC <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  p <- barplot(go_CC,showCategory=shownum,drop=T) +
    ggtitle("Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(opt$output_loc, "GOenrich_Cellular_Component.pdf"), plot = p, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(goopt == "MF"){
  go_MF <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  p <- barplot(go_MF,showCategory=shownum,drop=T) +
    ggtitle("Molecular Function") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(opt$output_loc, "GOenrich_Molecular_Function.pdf"), plot = p, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
# if(goopt == "ALL"){
#   go_BP <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
#                     qvalueCutoff = 0.2,keyType = gene_type)
#   go_CC <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
#                     qvalueCutoff = 0.2,keyType = gene_type)
#   go_MF <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
#                     qvalueCutoff = 0.2,keyType = gene_type)
#   p1 <- barplot(go_BP,showCategory=shownum,drop=T) +
#     ggtitle("Biological Process") +
#     theme(plot.title = element_text(hjust = 0.5, size = 20))
#   p2 <- barplot(go_CC,showCategory=shownum,drop=T) +
#     ggtitle("Cellular Component") +
#     theme(plot.title = element_text(hjust = 0.5, size = 20))
#   p3 <- barplot(go_MF,showCategory=shownum,drop=T) +
#     ggtitle("Molecular Function") +
#     theme(plot.title = element_text(hjust = 0.5, size = 20))
#   p <- cowplot::plot_grid(p1, p2, p3, nrow=3)
#   ggsave("~/temp/GOenrich_All.pdf", plot = p, 
#          device = "pdf", width = 20, height = 3*15/50*shownum, limitsize = FALSE)
# }


library(tidyverse)
library(clusterProfiler)
library(cowplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(optparse)

option_list <- list( 
  make_option(c("-r", "--rdata_file"), action="store"),
  make_option(c("-f", "--fc_thresh"), action="store", default = 2.0), 
  make_option(c("-p", "--p_thresh"), action="store", default = 0.05), 
  
  make_option(c("-t", "--species"), action="store", default = "mmu"), 
  make_option(c("-g", "--gene_type"), action="store", default = "ENSEMBL"), 
  make_option(c("-s", "--show_num"), action="store", default = 50), 
  make_option(c("-c", "--go_term"), action="store", default = "Biological_Process"), 
  
  make_option(c("-o", "--output_loc"), action="store")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

load(paste0(opt$rdata_file, "data.RData"))
#After running DEAnalysis(using [expr group] option), get DEAresult
resLFC <- DEAresult$resLFC
experGrp <- DEAresult$experGrp
controlGrp <- DEAresult$dataProc$dataSet$controlGrp
if(data_type == "MiAr"){
  resLFC <- resLFC %>%
    dplyr::rename(log2FoldChange = logFC, 
                  padj = adj.P.Val) 
}
#!!!Client options: Fold change threshold
fcthresh <- opt$fc_thresh
#!!!Client options: Fold change p.value
pthresh <- opt$p_thresh
target.data.up <- resLFC %>%
  as.data.frame() %>%
  #rownames_to_column(var = "lipid") %>%
  filter(
    log2FoldChange > log(fcthresh, 2) & padj < pthresh
  )
target.data.down <- resLFC %>%
  as.data.frame() %>%
  #rownames_to_column(var = "lipid") %>%
  filter(
    log2FoldChange < -log(fcthresh, 2) & padj < pthresh
  )
gene_up <- target.data.up$gene
gene_down <- target.data.down$gene

#!!!Client options: species("hsa"/"mmu")
spe <- opt$species
#!!!Client options: gene type(sugg: "ENSEMBL"/"SYMBOL")(ref: keytypes([spe db]))
gene_type <- opt$gene_type
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
  go_BP_up <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  go_BP_down <- enrichGO(gene = gene_down, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.2,keyType = gene_type)
  p1 <- barplot(go_BP_up,showCategory=shownum,drop=T) +
    ggtitle("Up reg genes - Biological Process") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  p2 <- barplot(go_BP_down,showCategory=shownum,drop=T) +
    ggtitle("Down reg genes - Biological Process") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(opt$output_loc, "upreggenes_GOenrich_Biological_Process.pdf"), plot = p1, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
  ggsave(paste0(opt$output_loc, "downreggenes_GOenrich_Biological_Process.pdf"), plot = p2, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(goopt == "CC"){
  go_CC_up <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  go_CC_down <- enrichGO(gene = gene_down, OrgDb = orgdb, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2,keyType = gene_type)
  p1 <- barplot(go_CC_up,showCategory=shownum,drop=T) +
    ggtitle("Up reg genes - Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  p2 <- barplot(go_CC_down,showCategory=shownum,drop=T) +
    ggtitle("Down reg genes - Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(opt$output_loc, "upreggenes_GOenrich_Cellular_Component.pdf"), plot = p1, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
  ggsave(paste0(opt$output_loc, "downreggenes_GOenrich_Cellular_Component.pdf"), plot = p2, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(goopt == "MF"){
  go_MF_up <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  go_MF_down <- enrichGO(gene = gene_down, OrgDb = orgdb, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2,keyType = gene_type)
  p1 <- barplot(go_MF_up,showCategory=shownum,drop=T) +
    ggtitle("Molecular Function") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  p2 <- barplot(go_MF_down,showCategory=shownum,drop=T) +
    ggtitle("Molecular Function") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(opt$output_loc, "upreggenes_GOenrich_Molecular_Function.pdf"), plot = p1, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
  ggsave(paste0(opt$output_loc, "downreggenes_GOenrich_Molecular_Function.pdf"), plot = p2, 
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

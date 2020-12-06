library(clusterProfiler)
library(cowplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#Input
gene_up <- data$To
#!!!Client options: species("hsa"/"mmu")
spe <- "mmu"
#!!!Client options: gene type(sugg: "ENSEMBL"/"SYMBOL")(ref: keytypes([spe db]))
gene_type <- "ENSEMBL"
#!!!Client options: show ontology item number
shownum <- 50
#!!!Client options: select go term catalog("BP"/"CC"/"MF"/"ALL")
gocat <- "BP"
orgdb <- switch (spe,
                 hsa = org.Hs.eg.db, 
                 mmu = org.Mm.eg.db
)
goopt <- switch (gocat,
                 BP = "BP", 
                 CC = "CC",
                 MF = "MF", 
                 ALL = "ALL"
)
if(gocat == "BP"){
  go_BP <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  p <- barplot(go_BP,showCategory=shownum,drop=T) +
    ggtitle("Biological Process") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave("~/temp/GOenrich_Biological_Process.pdf", plot = p, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(gocat == "CC"){
  go_CC <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  p <- barplot(go_CC,showCategory=shownum,drop=T) +
    ggtitle("Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave("~/temp/GOenrich_Cellular_Component.pdf", plot = p, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(gocat == "MF"){
  go_MF <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = gene_type)
  p <- barplot(go_MF,showCategory=shownum,drop=T) +
    ggtitle("Molecular Function") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave("~/temp/GOenrich_Molecular_Function.pdf", plot = p, 
         device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)
}
if(gocat == "ALL"){
  go_BP <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,keyType = gene_type)
  go_CC <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,keyType = gene_type)
  go_MF <- enrichGO(gene = gene_up, OrgDb = orgdb, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,keyType = gene_type)
  p1 <- barplot(go_BP,showCategory=shownum,drop=T) +
    ggtitle("Biological Process") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  p2 <- barplot(go_CC,showCategory=shownum,drop=T) +
    ggtitle("Cellular Component") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  p3 <- barplot(go_MF,showCategory=shownum,drop=T) +
    ggtitle("Molecular Function") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  p <- cowplot::plot_grid(p1, p2, p3, nrow=3)
  ggsave("~/temp/GOenrich_All.pdf", plot = p, 
         device = "pdf", width = 20, height = 3*15/50*shownum, limitsize = FALSE)
}

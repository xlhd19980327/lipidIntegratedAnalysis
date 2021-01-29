geneEnrichFunc <- function(genes, spe, gene_type, shownum, gocat,
                           reg = "none", loc){
  library(clusterProfiler)
  library(cowplot)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggplot2)
  #!!!Client options: species("hsa"/"mmu")
  #spe <- opt$species
  #!!!Client options: gene type(sugg: "ENSEMBL"/"SYMBOL")(ref: keytypes([spe db]))
  #gene_type <- opt$gene_type
  #!!!Client options: show ontology item number
  #shownum <- opt$show_num
  #!!!Client options: select go term catalog("BP"/"CC"/"MF"/"ALL")
  #gocat <- opt$go_term
  regname <- ifelse(reg == "none", "", paste0(reg, "_"))
  orgdb <- switch (spe,
                   hsa = org.Hs.eg.db, 
                   mmu = org.Mm.eg.db
  )
  goopt <- switch (gocat,
                   Biological_Process = "BP", 
                   Cellular_Component = "CC",
                   Molecular_Function = "MF"
  )
  gotitle <- switch (gocat,
                     Biological_Process = "Biological Process", 
                     Cellular_Component = "Cellular Component",
                     Molecular_Function = "Molecular Function"
  )
  go <- enrichGO(gene = genes, OrgDb = orgdb, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.2,keyType = gene_type)
  if(is.null(go)){
    regstat <- switch (reg,
      up = "UP", 
      down = "DOWN", 
      none = ""
    )
    cat(paste0("No ", regstat, 
               " genes enriched! Try check your data!\n"))
    return(1)
  }else{
    p1 <- barplot(go,showCategory=shownum,drop=T, font.size=15) +
      ggtitle(gotitle) +
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    ggsave(paste0(loc, regname, "GOenrich_", gocat, ".pdf"), plot = p1, 
           device = "pdf", width = 12, height = 15/50*shownum, limitsize = FALSE)
    
  }
}
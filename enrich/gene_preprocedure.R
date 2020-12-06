#After running DEAnalysis(using [expr group] option), get DEAresult
resLFC <- DEAresult$resLFC
experGrp <- DEAresult$experGrp
controlGrp <- DEAresult$dataProc$dataSet$controlGrp
if(type == "MiAr"){
  resLFC <- resLFC %>%
    dplyr::rename(log2FoldChange = logFC, 
                  padj = adj.P.Val) 
}
#!!!Client options: Fold change threshold
fcthresh <- 2.0
#!!!Client options: Fold change p.value
pthresh <- 0.05
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


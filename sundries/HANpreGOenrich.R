library(clusterProfiler)
library(cowplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
a <- read.csv("~/temp/cor/genes_2.csv")
a2 <- read.csv("~/temp/cor2/genes_5.csv")
al <- rbind(a, a2)
write.csv(al, "~/temp/cor/genes_100.csv", row.names = F)
library(org.Mm.eg.db)
go <- clusterProfiler::enrichGO(a2$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = "SYMBOL")
shownum <- 20
p <- barplot(go,showCategory=shownum,drop=T) +
  ggtitle("Biological Process") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
ggsave("~/temp/GOenrich_Biological_Process.pdf", plot = p, 
       device = "pdf", width = 20, height = 15/50*shownum, limitsize = FALSE)


rownames(correlation)
rvd2 <- correlation[rownames(correlation) == "RvD2 ", ]
up <- names(rvd2[rvd2 > 0])
down <- names(rvd2[rvd2 < 0])
library(org.Mm.eg.db)
go <- clusterProfiler::enrichGO(down, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.2,keyType = "SYMBOL")
shownum <- 20
p <- barplot(go,showCategory=shownum,drop=T) +
  ggtitle("Biological Process") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
  
ggsave("~/temp/GOenrich_Biological_Process.pdf", plot = p, 
       device = "pdf", width = 7, height = 4, limitsize = FALSE)
write.csv(down, "~/temp/cor/genes_100.csv", row.names = F)

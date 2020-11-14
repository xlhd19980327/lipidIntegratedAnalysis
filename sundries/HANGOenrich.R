
install.packages("gprofiler2")
library(gprofiler2)
a <- read.csv("~/temp/up_regup_gene.csv")
b <- gost(query = a$x, organism = 'mmusculus')
d <- gostplot(b, interactive = F)
e <- publish_gostplot(d)
f <- publish_gosttable(b)

library(clusterProfiler)
a <- read.csv("~/temp/up_regup_gene.csv")
go <- enrichGO(a$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = "ENSEMBL")
barplot(go,showCategory=50,drop=T)

a2 <- read.csv("~/temp/up_regdown_gene.csv")
go2 <- enrichGO(a2$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = "ENSEMBL")
barplot(go2,showCategory=20,drop=T)

a3 <- read.csv("~/temp/down_regup_gene.csv")
go3 <- enrichGO(a3$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                qvalueCutoff = 0.2,keyType = "ENSEMBL")
barplot(go3,showCategory=50,drop=T)

a4 <- read.csv("~/temp/down_regdown_gene.csv")
go4 <- enrichGO(a4$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                qvalueCutoff = 0.2,keyType = "ENSEMBL")
barplot(go4,showCategory=50,drop=T)

a5 <- read.csv("~/temp/correlation_rvd2_regup_gene.csv")
go5 <- enrichGO(a5$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                qvalueCutoff = 0.2,keyType = "ENSEMBL")
barplot(go5,showCategory=50,drop=T)

a6 <- read.csv("~/temp/correlation_rvd2_regdown_gene.csv")
go6 <- enrichGO(a6$x, OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                qvalueCutoff = 0.2,keyType = "ENSEMBL")
barplot(go6,showCategory=50,drop=T)

## For visualize the reactomeDB pathway
# options(stringsAsFactors = F)
# results <- read.csv("./testData/SVF191222/output/regStat_D7D.csv")
# db <- read.csv("./patternHunting/lipidReatomeStat/hsaDB/hsa_all_integData.csv")
# db2 <- read.csv("./patternHunting/lipidReatomeStat/hsaDB/hsa_lipidreact.csv")
# db <- db %>% 
#   rename(gene = gene_symbol) %>%
#   select(reaction_id, gene) %>%
#   unique()
# db2 <- db2 %>%
#   select(reaction_id, firstLev, secondLev) %>%
#   unique()
# integ <- inner_join(db2, db)
# integ2 <- inner_join(results, integ) %>%
#   select(gene, P_R, p_value, firstLev, secondLev) %>%
#   unique()
# write.csv(integ2, "~/temp/integ2.csv", row.names = F)


#library(ggplot2)
#setwd("D:/downloadedTemp")
#library(ggsci)
library(scales)
#source("./utilityFunc/plottingPalettes.R")
options(stringsAsFactors = F)
data <- read.csv("~/temp/vis.csv")
data <- unique(data)
data2 <- cbind(data, log10p = -log10(data$p_value))
#color <- plottingPalettes(100, type = "continuous")
ggplot(data = data2, aes(x = P_R, y = reorder(gene, P_R), fill = log10p)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_gradient2(low = muted("blue"),
                       mid = "white",
                       high = muted("red"),
                       midpoint = 0)+
  labs(x = "log2FC", y = "gene", fill = "-log10PValue")
  #scale_fill_gradientn(colours = color)

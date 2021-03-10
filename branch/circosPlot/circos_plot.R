options(stringsAsFactors = F)
go <<- go

cordata <- read.csv("~/temp/cor2/correlation_1_4.csv", 
                    row.names = 1)
lipids <- rownames(cordata)
genes <- colnames(cordata)
godata <- go_BP@result
thresh1 <- 0.6
subdatas <- lapply(cordata, function(x) lipids[abs(x) > thresh1])
subcors <- lapply(cordata, function(x) x[abs(x) > thresh1])
thresh2 <- 0.05
topnum <- 20
subgodata <- subset(godata, p.adjust < thresh2)
subgodata <- subgodata[order(subgodata$p.adjust), ][1:topnum, ]
subgogenes <- strsplit(subgodata$geneID, "/")

library(paletteer)
library(circlize)
gene_corfilt <- unlist(mapply(rep, names(subcors), ifelse(sapply(subcors, length) == 0, 0, 1)))
gene_gofilt <- unlist(subgogenes)
gene_filt <- gene_corfilt[gene_corfilt %in% gene_gofilt]
lipid_filt <- subdatas[gene_filt]
gl_link <- mapply(cbind, names(lipid_filt), lipid_filt)
gl_link <- do.call(rbind, gl_link)
gl_link_ind <- data.frame(
  geneind = match(gl_link[, 1], gene_filt), 
  lipind = match(gl_link[, 2], unique(unlist(lipid_filt))),
  value = mapply(function(x, y){cordata[x, y]}, gl_link[, 2], gl_link[, 1])
)
colind1 <- match(gl_link_ind$value, sort(unique(gl_link_ind$value)))
col1 <- as.character(paletteer_c("ggthemes::Classic Area Green", n = length(unique(gl_link_ind$value))))
#col1 <- gsub("FF", "", col1)
gl_link_ind <- cbind(gl_link_ind, 
                     col = col1[colind1])
go_link <- lapply(subgogenes, function(x){
  ind <- match(x, gene_filt)
  ind <- ind[!is.na(ind)]
  return(ind)
})
subgodata_item <- subgodata$Description
subgo_filt <- subgodata_item[unlist(mapply(rep, 1:length(subgodata_item), sapply(go_link, function(x) ifelse(length(x) == 0, 0, 1))))]
indtemp <- sapply(go_link, length)
indtemp <- indtemp[indtemp != 0]
go_link_ind <- data.frame(
  goind = unlist(mapply(rep, 1:length(subgo_filt), indtemp)), 
  geneind = unlist(go_link)
)
go_link_ind <- cbind(go_link_ind, 
                     value = subgodata$p.adjust[go_link_ind$goind])
colind2 <- match(go_link_ind$value, sort(unique(go_link_ind$value), decreasing = T))
col2 <- as.character(paletteer_c("ggthemes::Classic Area Red", n = length(unique(go_link_ind$value))))
#col2 <- gsub("FF", "", col2)
go_link_ind <- cbind(go_link_ind, 
                     col = col2[colind2])

sectors <- c("Genes", "Lipids", "Gene_Ontology")
n1 <- length(gene_filt)
n2 <- length(unique(unlist(lipid_filt)))
n3 <- length(subgo_filt)
df <- data.frame(
  sectors = sectors, 
  xlim = c(n1, n2, n3)
)
circos.clear()
pdf("~/temp/aaa.pdf")
circos.initialize(factor(df$sectors, levels = df$sectors), 
                  xlim = matrix(c(rep(0,3),df$xlim),ncol=2))
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  item=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim),mean(ylim),item,cex=0.9,col="grey40",
              facing="bending.inside",niceFacing=TRUE, font = 2)
},bg.col="grey90",bg.border=F,track.height=0.08)
circos.par("track.margin" = c(0, 0.1))
circos.par("cell.padding" = c(0, 1, 0, 1))
circos.track(ylim = c(0,10), panel.fun=function(x,y) {
  features <- CELL_META$xrange
  data <- switch(CELL_META$sector.index, 
                 Genes = gene_filt, 
                 Lipids = unique(unlist(lipid_filt)), 
                 Gene_Ontology = subgodata_item)
  for(i in 1:features){
    circos.text(i, 0, data[i], facing = "clockwise", adj = c(0, 0.5),
                cex = 0.5, niceFacing = T, font = 2)
  }
}, bg.border=F)
circos.par("track.margin" = c(0, 0))
for(i in 1:nrow(gl_link_ind)){
  circos.link("Genes", gl_link_ind$geneind[i], "Lipids", gl_link_ind$lipind[i], 
              col = gl_link_ind$col[i])
}
for(i in 1:nrow(go_link_ind)){
  circos.link("Gene_Ontology", go_link_ind$goind[i], "Genes", go_link_ind$geneind[i], 
              col = go_link_ind$col[i])
}
dev.off()







# Libraries
library(ggraph) #### Useful for network plot but not rounded for circos plot
library(igraph)
library(tidyverse)

# create a data frame giving the hierarchical structure of your individuals
d1 <- data.frame(from = "origin", to = c("Lipids", "Genes", "Gene_Ontology"))
d2 <- data.frame(from = unlist(mapply(rep, d1$to, c(length(lipids), length(genes), nrow(godata)))), 
                 to = c(lipids, genes, godata$Description))
hierarchy <- rbind(d1, d2)
# create a dataframe with connection between leaves (individuals)
connect1 <- data.frame(
  from = unlist(mapply(rep, names(subcors), sapply(subcors, length))), 
  to = unlist(subdatas), 
  value = unlist(subcors)
)
connect2 <- data.frame(
  from = unlist(subgogenes), 
  to = unlist(mapply(rep, subgodata$Description, sapply(subgogenes, length))), 
  value = unlist(mapply(rep, subgodata$p.adjust, sapply(subgogenes, length)))
)
rownames(connect2) <- NULL
connect <- rbind(connect1, connect2)
# create a vertices data.frame. One line per object of our hierarchy
vertices  <-  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))) , 
  value = 1
) 
# Let's add a column with the group of each name. It will be useful later to color points
vertices$group  <-  hierarchy$from[ match( vertices$name, hierarchy$to ) ]
vertices$value <- c(c('O', 'A', 'B', 'C'), unlist(lapply(vertices$group, function(x){switch(x, 
                                               Lipids = 'A',
                                               Genes = 'B', 
                                               Gene_Ontology = 'C')})))
# Create a graph object
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )
# The connection object must refer to the ids of the leaves:
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)
mycolor <- 
# Basic graph
p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", tension = .5) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  theme_void()

p +  
  geom_conn_bundle(data = get_con(from = from, to = to), aes(colour=value)) +
  scale_edge_colour_discrete()







sectors <- c(unique(unlist(lipid_filt)), gene_filt, subgodata_item)
n <- length(sectors)
df <- data.frame(
  sectors = rep(sectors, 2),
  x = c(rep(-1, n), rep(1, n))
)
circos.clear()
circos.par("cell.padding" = c(0, 0,
                            0, 0))
#circos.par("track.height" = 0.1)
#circos.par("gap.degree" = 0.1)
circos.par("circle.margin" = 3)
circos.initialize(df$sectors, x = df$x)
circos.track(ylim = c(-1,1), 
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(5), 
                           CELL_META$sector.index, 
                           facing = "clockwise", 
                           adj = c(0, 0))
               circos.axis(labels.cex = 0.6)
             })
circos.link("TG(18:1_24:1_24:1)", 0, "Bet1", 0, h = 0.4)

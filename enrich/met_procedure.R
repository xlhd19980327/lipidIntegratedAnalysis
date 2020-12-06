cmpd.vec <- read.csv("~/temp/cmpd.csv")$x
## cmpd.vec will be the input

library(MetaboAnalystR)
mSet2<-InitDataObjects("conc", "mSet2ora", FALSE)
mSet2<-Setup.MapData(mSet2, cmpd.vec);
mSet2<-CrossReferencing(mSet2, "name");
mSet2<-CreateMappingResultTable(mSet2)
queryn <- nrow(mSet2[["dataSet"]][["map.table"]])
queryNA <- which(mSet2[["dataSet"]][["map.table"]][, 2] == "NA")
if(length(queryNA) / queryn > 0.2 | queryn < 40){
  queryNAnames <- mSet2[["dataSet"]][["map.table"]][, 1][queryNA]
  for(i in queryNAnames){
    ## Change incorrect mapping names
    mSet2<-PerformDetailMatch(mSet2, i);
    mSet2<-GetCandidateList(mSet2);
    mSet2<-SetCandidate(mSet2, i, 
                        ifelse(mSet2[["name.map"]][["hits.candidate.list"]][1, 1] != '', 
                               mSet2[["name.map"]][["hits.candidate.list"]][1, 1], "NA"));
    
  }
}
## Enrichment
mSet2<-SetMetabolomeFilter(mSet2, F);
#kegg(may add other enrichment options)
mSet2<-SetCurrentMsetLib(mSet2, "kegg_pathway", 2);
#output a chart
mSet2<-CalculateHyperScore(mSet2)
#bar plot
mSet2<-PlotORA(mSet2, "ora_0_", "net", "png", 72, width=NA)
#dot plot
mSet2<-PlotEnrichDotPlot(mSet2, "ora", "ora_dot_0_", "png", 72, width=NA)

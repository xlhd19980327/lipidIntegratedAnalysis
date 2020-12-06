##### Output the result #####
#!!!Client options: output file location(using absolute file loc)
fileLoc <- "/home/lifs/temp/"

### LION input
write.csv(myoutput$value, paste0(fileLoc, "LION_input.csv"), row.names = F)

### LION enrichment table(enrichment summary)
write.csv(myoutput$ontology.table, paste0(fileLoc, "LION_enrichment.csv"), row.names = F)

### LION enrichment table(enrichment report)
write.csv(myoutput$downloadDetailTable$to_display_detailed, 
          file = paste0(fileLoc, "LION-enrichment-detailed-job.csv"), row.names = FALSE,quote = TRUE)
write.csv(myoutput$downloadDetailTable$lipidReport, 
          file = paste0(fileLoc, "LION-lipid-associations.csv"), row.names = FALSE,quote = TRUE)
write.csv(myoutput$downloadDetailTable$lipidsInTerms, 
          file = paste0(fileLoc, "LION-term-associations.csv"), row.names = FALSE,quote = TRUE)
# zip(zipfile="~/temp/LION-enrichment-report.zip", 
#     files=c(paste0(fileLoc, "LION-enrichment-detailed-job.csv"), 
#             paste0(fileLoc, "LION-lipid-associations.csv"),
#             paste0(fileLoc, "LION-term-associations.csv")))

### LION enrichment graph
ggsave(paste0(fileLoc, "LION-enrichment-plot.png"),
       width = 2.5*5.52, height = 2.5*3.66,
       myoutput$ontology.graph+
         theme(plot.subtitle = element_text(hjust = 0.5),
               plot.title = element_text(hjust = 0.5)),   
       device = 'png', dpi=300)
ggsave(paste0(fileLoc, "LION-enrichment-plot.svg"),
       width = 2.5*5.52, height = 2.5*3.66,
       myoutput$ontology.graph+
         theme(plot.subtitle = element_text(hjust = 0.5),
               plot.title = element_text(hjust = 0.5)),   
       device = svg,   scale=1.5)

### LION network view
# net = graph_from_data_frame(d = edges ,vertices = nodes,directed = T)
# layout <- layout_as_tree(net, flip.y = F, root = 1 , rootlevel = nodes$level)[,c(2,1)]
# layout[,1] <- unlist(coords_base[,1]) / 200
# layout[,2] <- unlist(coords_base[,2]) / 200
# plot(net, layout = layout,
#      vertex.label.cex=.4,
#      vertex.label.degree = 1.5*pi, edge.arrow.size = .40,
#      vertex.label.color="black")
visSave(myoutput$ontology.plot, paste0(fileLoc, "LION-network.html"))

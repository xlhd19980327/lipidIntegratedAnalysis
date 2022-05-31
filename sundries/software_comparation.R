options(stringsAsFactors = F)
orig <- read.csv("./testData/SVFmultiomics_210118/input/lipids.csv")
#if(dataType == "LipidSearch"){
source("./utilityFunc/getFAsInfo.R")
#}
lipids <- orig$LipidIon
fainfo <- lapply(lipids, getFAsInfo, F)

## For lipidr ##
#data conversion
lipids_tidy <- sapply(fainfo, function(x){
  fa <- x[[1]]
  Class <- names(x[1])
  lip <- paste0(Class, " ", paste(fa, collapse = "/"))
  return(lip)
})
data <- orig
data$LipidIon <- lipids_tidy
library(lipidr)
dataObj <- as_lipidomics_experiment(data)
#calculate lipids not in db
def <- lipidr:::.myDataEnv$lipidDefaults$clean_mols
molecules <- unique(as.character(lipids_tidy))
not_in_db <- molecules[!molecules %in% def$Molecule]
clean_ <- lipidr:::.clean_molecule_name(not_in_db)
not_in_db_num <- length(clean_$Molecule[clean_$not_matched])
not_in_db_num
dataanno <- read.csv("./testData/SVFmultiomics_210118/input/sampleList.csv")
dataanno_tidy <- data.frame(Sample = dataanno$samples, 
                            Condition = dataanno$conditions)
dataObj <- add_sample_annotation(dataObj, dataanno_tidy)
write.csv(data, "~/temp/out.csv", row.names = F)
#enrichment
d_normalized <- normalize_pqn(dataObj, measure = "Area", exclude = "blank", log = TRUE)
de_results = de_analysis(
  data=d_normalized, 
  Day10 - Day0,
  measure="Area"
)
enrich_results = lsea(de_results, rank.by = "logFC")
significant_lipidsets(enrich_results)
write.csv(enrich_results[, -9], "~/temp/out5.csv", row.names = F)
#plot_enrichment(de_results, significant_lipidsets(enrich_results), annotation="class")

## For LIPEA ##
lipids_tidy <- sapply(fainfo, function(x){
  fa <- x[[1]]
  Class <- names(x[1])
  peinfo <- ifelse(grepl("\\(", Class), gsub(".*?\\((O|P)\\)", "\\1- ", Class), "")
  Class <- ifelse(grepl("\\(", Class), gsub("\\([O|P]\\)", "", Class), Class)
  Class <- ifelse(grepl("DG", Class), gsub("DG", "DAG", Class), Class)
  Class <- ifelse(grepl("TG", Class), gsub("TG", "TAG", Class), Class)
  Class <- ifelse(grepl("MG", Class), gsub("MG", "MAG", Class), Class)
  Class <- ifelse(grepl("ChE", Class), gsub("ChE", "Cholesterol", Class), Class)
  lip <- paste0(Class, " ", peinfo, paste(fa, collapse = "/"))
  return(lip)
})
write.csv(lipids_tidy, "~/temp/out2.csv")
#enrichment
volcano_data <- read.csv("./testData/SVFmultiomics_210118/output/MARresults/volcano_data.csv")
volcano_data <- volcano_data %>%
  mutate(p = 10 ^ -p.log, 
         fc = 2 ^ fc.log) %>%
  filter(p < 0.1 & fc >2)
lipids <- volcano_data$X
fainfo <- lapply(lipids, getFAsInfo, F)
lipids_volcanofilt <- sapply(fainfo, function(x){
  fa <- x[[1]]
  Class <- names(x[1])
  peinfo <- ifelse(grepl("\\(", Class), gsub(".*?\\((O|P)\\)", "\\1- ", Class), "")
  Class <- ifelse(grepl("\\(", Class), gsub("\\([O|P]\\)", "", Class), Class)
  Class <- ifelse(grepl("DG", Class), gsub("DG", "DAG", Class), Class)
  Class <- ifelse(grepl("TG", Class), gsub("TG", "TAG", Class), Class)
  Class <- ifelse(grepl("MG", Class), gsub("MG", "MAG", Class), Class)
  Class <- ifelse(grepl("ChE", Class), gsub("ChE", "Cholesterol", Class), Class)
  lip <- paste0(Class, " ", peinfo, paste(fa, collapse = "/"))
  return(lip)
})
write.csv(lipids_volcanofilt, "~/temp/out4.csv")

## For MAR ##
lipids_tidy <- gsub("_", "/", lipids)
write.csv(lipids_tidy, "~/temp/out3.csv")
mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, lipids_tidy);
mSet<-CrossReferencing(mSet, "name", lipid = T);
mSet<-CreateMappingResultTable(mSet)
in_db_num <- sum(mSet[["dataSet"]][["map.table"]][, 7] == "1")
nrow(orig) - in_db_num


## For BioPAN ##
#See "~/myWork/github/lipidIntegratedAnalysis/sundries/bioPAN.R"
#can use Day8 vs Day0
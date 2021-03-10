## Lipidr comparation
#options(stringsAsFactors = F)
#library(lipidr)
#data <- read.csv("./sundries/comparationData/HANlipid_tidy_lipidr.csv")
#dataObj <- as_lipidomics_experiment(data)
#dataanno <- read.csv("./sundries/comparationData/HANsampleList_lipid_lipidr.CSV")
#dataObj <- add_sample_annotation(dataObj, dataanno)

## Prepare lipids(from intra-omics analysis with strong correlation)
data1 <- read.csv("./branch/benchmark/input/HANlipid_onlyforallclass.csv")




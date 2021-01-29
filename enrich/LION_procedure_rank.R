##!!!Client options: rank parameter("p.value"/"fc.log")
rankarg <- "p.value"

### Preparation
#### libraries
require(shiny)
require(visNetwork)
require(data.table)
require(igraph)
require(ggplot2)
require(ggthemes)
library(shinyTree)
library(shinyWidgets)
library(shinyBS)
library(httr)
library(formattable)
library(jsonlite)
library(ggrepel)
library(shinycssloaders)
## loading lipid ontology data
require(RSQLite)
require(topOnto)
require('topOnto.LION.db')
topOnto::initONT('LION')
source("./../LION-web-master/OntologyApp/data/20191008 LIONweb functions.R")
associationFile  <-  "./../LION-web-master/OntologyApp/data/20190704 LION_association.txt"
## read associations
lipidID2TERM <- readMappings(file = associationFile)   ## topOnto function, but removes spaces
LIONterms_rules <- read.csv(file = './../LION-web-master/OntologyApp/data/20191008 LIONterms_rules.csv', header = T)
LIONterms_rules$RULE1[LIONterms_rules$RULE1 == ""] <- "-"
LIONterms_FAs <- read.csv(file = './../LION-web-master/OntologyApp/data/20191008 LIONterms_FAs.csv', header = T)

input <- list()
output <- list()
myoutput <- list()
input$RemoveRedundantTerms <- T
input$SmartMatching <- T
input$FAprediction <- F
input$EmailMissingLipids <- F
input$LIONselection <- F
input$ranking <- "ascending"
input$ks_sided <- "ks"

###813
if(isolate(input$ranking)=="ascending"){direction <- 1} else {direction <- -1}

###858-892
# pValueList <- list()
# pValueList$input <- lapply(data_input$lipid, function(x) x[1])
# pValueList$input <- paste(sprintf("[#%04d]", 1:length(pValueList$input)),pValueList$input)
lipidadd <- paste(sprintf("[#%04d]", 1:length(data_input$lipid)),data_input$lipid)
# pValueList$pvalues <- lapply(data_input[[rankarg]], function(x) x[1])
pValueList <- data.frame(input = lipidadd, pvalues = data_input[[rankarg]])
lipidExistance <- data.frame(input = pValueList$input)
lipidExistance$'simplified input' <- convertLipidNames(gsub("^\\[#\\d+\\] ","",lipidExistance$input))

###893-1079
if (length(lipidExistance$input) < 1) {
  
} else {
  
  ##### new 20191106
  lipidExistance_list <-
    lapply(1:dim(lipidExistance)[1], function(row_i) {
      
      lipid_i <- lipidExistance[row_i, 2]   ## lipid_i is lipid of this iteration
      
      LION_ID <-
        unlist(lipidID2TERM[lipid_i == names(lipidID2TERM)])
      
      if (is.null(LION_ID)) {
        ### if input is LION:xxxx
        LION_ID <-
          unlist(lipidID2TERM[lipid_i == unlist(lipidID2TERM)])
        LION_ID <-
          LION_ID[!grepl("^SLM:|^LM..\\d+", names(LION_ID))]   ## remove SwissLipids/LIPIDMAPS IDs
      }
      
      if (!is.null(LION_ID)) {
        output_df <-
          data.frame(
            input = lipidExistance[row_i, 1],
            'simplified input' = lipid_i,
            name = names(LION_ID),
            LION = LION_ID,
            match = "direct"
          )
      } else {
        ### no direct matching is found
        if (isolate(input$SmartMatching)) {
          ### 'smartmatching' is on
          lipid_index <-
            list(
              generalized = simplifyLipidnames(lipid_i),
              headgroup = getLipidHeadgroup(lipid_i),
              linkage = getLipidLinkage(lipid_i),
              FAs = getFAs(lipid_i)
            )
          
          generalized_LION_ID <-
            unlist(lipidID2TERM[lipid_index$generalized == names(lipidID2TERM)])
          if (!is.null(generalized_LION_ID)) {
            terms <-
              data.frame(name = names(lipidID2TERM)[lipid_index$generalized == names(lipidID2TERM)],
                         LION = generalized_LION_ID)
          } else {
            ## no generalized lipid found
            LIONterms_rules_i <-
              LIONterms_rules[lipid_index$headgroup == LIONterms_rules$RULE1 ,]
            
            if (all(LIONterms_rules_i$RULE2 == "")) {
              terms <- LIONterms_rules_i[, 1:2]
            } else {
              terms <-
                LIONterms_rules_i[lipid_index$linkage == LIONterms_rules_i$RULE2,][, 1:2]
            }
          }
          
          
          terms <-
            rbind(terms, LIONterms_FAs[LIONterms_FAs$name %in% lipid_index$FAs,])
          if (dim(terms)[1] > 0) {
            output_df <-
              data.frame(input = lipidExistance[row_i, 1],
                         'simplified input' = lipid_i,
                         terms,
                         match = "smart matching")
            
          } else {
            ## no match by smart matching
            output_df <-
              data.frame(
                input = lipidExistance[row_i, 1],
                'simplified input' = lipid_i,
                name = "not found",
                LION = "not found",
                match = ""
              )
          }
          
          
        } else {
          ### 'smartmatching' is off, and no match found
          output_df <-
            data.frame(
              input = lipidExistance[row_i, 1],
              'simplified input' = lipid_i,
              name = "not found",
              LION = "not found",
              match = ""
            )
        }
        
        
      }
      if(isolate(input$FAprediction)){    ## predict FAs if nescerrary
        
        lipid_index <-
          list(
            generalized = simplifyLipidnames(lipid_i),
            headgroup = getLipidHeadgroup(lipid_i),
            FAs = getFAs(lipid_i)
          )
        
        if(lipid_index$generalized == lipid_i &   length(lipid_index$FAs) == 1){   ## is FA prediction applicable?
          
          predicted_FAs <- predict_FAs(headgroup = lipid_index$headgroup, 
                                       summedFA = gsub("C","",lipid_index$FAs), 
                                       composition_table = FA_composition_table )
          
          if(dim(LIONterms_FAs[LIONterms_FAs$name %in% predicted_FAs,])[1]>0){    ## any result?
            output_df <-
              rbind(output_df,
                    data.frame(input = lipidExistance[row_i, 1],
                               'simplified input' = lipid_i,
                               LIONterms_FAs[LIONterms_FAs$name %in% predicted_FAs,],
                               match = "FA-prediction"))
          } 
          
          
        }
        
        
      }
      
      return(output_df)
    })
  
  matching_statistics <- data.frame(total = length(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})),
                                    matched = sum(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})),
                                    percent = mean(sapply(lipidExistance_list, function(lipid_i){!any(lipid_i$LION == "not found")})    ) * 100)
  
  
  lipidExistance <- do.call("rbind",lipidExistance_list)
  
  colnames(lipidExistance) <- c('input','simplified input','LION name','LION ID', "match")
  
  
  color_df <- data.frame(match = c("direct","smart matching", "FA-prediction",""),
                         color = c("#4f4d4d","#9c3c2d","#c4942b","#bdbbbf"))
  
  lipidExistance_toview <-
    do.call("rbind", lapply(lipidExistance_list, function(lipidExistance_i) {
      color_pattern <- color_df$color[match(lipidExistance_i$match, color_df$match)]
      
      data.frame(
        input = unique(lipidExistance_i$input),
        name = paste("<font color='",color_pattern,"'>",lipidExistance_i$name,"</font>" , sep ="", collapse = "<br>"),
        #LION = paste("<font color='",color_pattern,"'>",lipidExistance_i$LION,"</font>" , sep ="",  collapse = "<br>")
        LION = paste("<a href='",
                     'http://bioportal.bioontology.org/ontologies/LION/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F',
                     gsub(":","_",lipidExistance_i$LION),
                     "' style='color: ", color_pattern,
                     "' target='_blank'>",
                     lipidExistance_i$LION,
                     "</a>",
                     sep = "", collapse = "<br>")
        
        
      )
    }))
  
  
  colnames(lipidExistance_toview) <- c('input','LION name','LION ID')
  
  lipidExistance_toview[[3]] <- gsub(" href=.+Fnot found'","",lipidExistance_toview[[3]])
  
  ##### end new 20191006
  
  lipidExistance$'simplified input' <- NULL        ## remove this column, not interesting for user
  
  
  ## use matched lipids as assocation table 
  
  lipidExistance_feasable <- lipidExistance[lipidExistance$`LION ID` != "not found",]
  
  matched_lipidID2TERM <- 
    sapply(unique(lipidExistance_feasable$input), function(input){
      lipidExistance_feasable$`LION ID`[lipidExistance_feasable$input == input]
    }, simplify = FALSE)
  
  
}

###1083:1259-1448
pValueList$pvalues <- as.numeric(pValueList$pvalues)

lipidIDrank <- rank(pValueList$pvalues * direction)    ## direction by radiobutton input$ranking
names(lipidIDrank) <- pValueList$input   

mySel <- function(allScore) {
  return(rep(TRUE, length(lipidIDrank)))
}


if(any(names(lipidIDrank) %in% names(matched_lipidID2TERM))){   ### can at least one lipid be matched?
  ONTdata <- new(
    "topONTdata",
    ontology = "LION",
    allGenes = lipidIDrank,
    annot = annFUN.gene2GO,
    gene2GO = matched_lipidID2TERM,
    geneSelectionFun = mySel
  )
} else {       ### make a fake object
  names(lipidIDrank)[1] <- 'PC(32:0)'
  ONTdata <- new(
    "topONTdata",
    ontology = "LION",
    allGenes = lipidIDrank,
    annot = annFUN.gene2GO,
    gene2GO = lipidID2TERM,
    geneSelectionFun = mySel
  )
}

resultFis <-
  runTest(ONTdata,
          algorithm = "classic",
          statistic = isolate(input$ks_sided))

resultFis@score <- abs(resultFis@score)

ES <- runTest(ONTdata,  algorithm = "classic",
              statistic = "ks.score")@score

sign_i_df <- data.frame(LION = names(ES), 
                        sign = sign(ES),
                        ES = ES)




#incProgress(.7, detail = paste("enrichment statistics"))

to_display <-
  GenTable(ONTdata,
           'p-value' = resultFis,
           topNodes = 2000)

to_display <-
  to_display[grep("LION", to_display$TERM.ID), ]


LUT <- data.frame(ID = names(ONTdata@termName),
                  Discription = ONTdata@termName)

to_display$Term <- LUT$Discription[match(to_display$TERM.ID, LUT$ID)]        ### otherwise, some names are abbrev.
to_display <- to_display[to_display$Annotated > 2, ]        ### remove terms with 1 or 2 lipids

####  limiting by LION-term selection
if(isolate(input$LIONselection)){    ### switch for LION-selection
  TermsOfInterest <- get_selected(isolate(input$tree), format = "names") 
  TermsOfInterest <- unique(c(unlist(TermsOfInterest),
                              unlist(sapply(TermsOfInterest, function(element){attr(element,"ancestry")}))
  )) ## vector with LION-names
  to_display <- to_display[to_display$Term %in% TermsOfInterest,]  
}


## limiting redundant/similar parent terms

## limiting redundant/similar parent terms

if(isolate(input$RemoveRedundantTerms)){    ### switch for LION-selection
  ONT_DAG <- ONTdata@graph
  ONT_DAG_lev <- buildLevels(ONT_DAG)
  DAG.env <- ONT_DAG_lev$nodes2level
  
  DAGlevel <-   sapply(to_display$TERM.ID, function(id) {
    DAG.env[[id]]
  })
  
  lipidsInTerms <- genesInTerm(ONTdata)
  
  TermContent <- 
    sapply(to_display$TERM.ID, function(id) {
      lipidsInTerms[[id]]
    })
  
  test_similarity <-
    sapply(names(TermContent), function(term_i) {
      sapply(names(TermContent), function(term_j) {
        if (term_i != term_j) {
          if (length(TermContent[[term_i]]) == length(TermContent[[term_j]])) {
            same <- all(TermContent[[term_i]] %in% TermContent[[term_j]])
          } else {
            same <- FALSE
          }
          
          
          outputList <- list(
            term_a = term_i,
            term_b = term_j,
            isSame = same,
            term_a_level = DAGlevel[names(DAGlevel) == term_i],
            term_b_level = DAGlevel[names(DAGlevel) == term_j]
          )
          
          
          outputList$remove <- ""
          
          if ((outputList$term_a_level - outputList$term_b_level)==1) {
            outputList$remove <- outputList$term_b
          } 
          if ((outputList$term_a_level - outputList$term_b_level)==-1) {
            outputList$remove <- outputList$term_a
          }
          
          outputList
          
        } else {
          list(
            term_a = term_i,
            term_b = term_j,
            remove = "",
            isSame = FALSE
          )
        }
      }, simplify = FALSE)
    })
  
  test_similarity <-
    test_similarity[unlist(lapply(test_similarity, function(n) {
      n[['isSame']]
    }))]
  
  test_similarity <- 
    unique(unlist(lapply(test_similarity, function(n) {
      n[['remove']]
    })))
  
  to_display <- to_display[!(to_display$TERM.ID %in% test_similarity),]  
  
}


###


to_display$`p-value` <- gsub("< 1e","< 1.0e",to_display$`p-value`)    ### '< 1e-30' cannot be understood
to_display$`p-value` <- gsub("< ","",to_display$`p-value`)    ### '< 1e-30' cannot be understood

to_display$'FDR q-value' <-
  format(p.adjust(as.numeric(to_display$`p-value`), "fdr"), digits = 3)

colnames(to_display) <-
  c(
    "Term ID",
    "Discription",
    "Annotated",
    "Significant",
    "Expected",
    "p-value",
    "FDR q-value"
  )
to_display <- to_display[, c(1, 2, 3, 6, 7)]

if(!is.null(sign_i_df) & isolate(input$ks_sided) == "ks2"){   ### show scores when 2-sided
  to_display$ES <- sign_i_df$ES[match(to_display$`Term ID`, sign_i_df$LION)]
  to_display$Regulated <- as.character(factor(sign_i_df$sign[match(to_display$`Term ID`, sign_i_df$LION)], levels = c(1,-1), labels = c("UP","DOWN")))
}

###1083:1449-1543
lengthOfTable <- length(to_display$"Term ID")

### make detailed table with lipids per term

to_display_detailed <- to_display

lipidsInTerms <- genesInTerm(ONTdata)

identifiers <- lapply(to_display_detailed$'Term ID', function(ID){
  lipidsInTerms[[ID]]
  #genesInTerm(ONTdata)[[ID]]
})

####  

if(length(to_display_detailed[,1]) > 0){     ### if this is zero, script goes wrong
  
  output <- lapply(1:length(to_display_detailed[,1]), function(row){   ### add rows with lipid info
    output <- as.data.frame(matrix(ncol = length(to_display_detailed[1,]),
                                   nrow = length(identifiers[[row]])+1,
                                   data=""))
    output[1,] <-   to_display_detailed[row,]
    output[-1,1] <- paste("   ",identifiers[[row]])
    
    output
  })
  
  to_display_detailed <- as.data.frame(rbindlist(output))
  colnames(to_display_detailed) <- colnames(to_display)
} else {to_display_detailed <- to_display}


### lipid associations report

lipidReport <- associatedTerms(lipid = ONTdata@allGenes, ontologyObject = ONTdata, reformat = TRUE)
colnames(lipidReport) <- lipidReport[1,]
lipidReport <- lipidReport[-1,]

### term associations report

#lipidsInTerms <- genesInTerm(ONTdata)

lipidsInTerms <- data.frame(ID = names(lipidsInTerms),
                            Discription = LUT$Discription[match(names(lipidsInTerms), LUT$ID)],
                            nr_lipids = sapply(lipidsInTerms, function(term){length(term)}),
                            lipids = sapply(lipidsInTerms, function(term){ paste(gsub("^\\[#\\d+\\] ","",term), collapse = '; ')   })
)

lipidsInTerms <- subset(lipidsInTerms, Discription != lipids)
colnames(lipidsInTerms) <- c("LION-term","LION-name", "nr_lipids", "lipid identifiers")

### ggplot of data
to_ggplot_display <- to_display


to_ggplot_display$color <-
  -log(as.numeric(to_ggplot_display$`FDR q-value`), base = 10)
to_ggplot_display$color[to_ggplot_display$color > 6] <- 6


if(!is.null(isolate(input$split_direction_in_output))){
  graph_option <- isolate(input$split_direction_in_output)
} else {
  graph_option <-"combine"
}


main_title <- substitute(bold("LION enrichment analysis")~
                           italic(x), list(x = tolower("Ranking mode")))

sub_title <- 
  c(ID = ifelse(
    isolate(input$file_input == "load external dataset") & isolate(input$method) == "Ranking mode",
    paste(isolate(input_data()$studyID)," ", sep =""),
    ""
  ),
  A = ifelse(isolate(input$method) == "Target-list mode","",paste("",isolate(input$conditionA), sep = "")),
  B = ifelse(isolate(input$method) == "Target-list mode","",paste("",isolate(input$conditionB), sep = "")),
  vs  = ifelse(is.null(isolate(input$conditionA)) | is.null(isolate(input$conditionB)) | isolate(input$method) == "Target-list mode",
               "", "vs."),
  mode = isolate(input$method))

sub_title <- 
  substitute(bold(studyID) ~ A ~ italic(vs)~ B, 
             list(studyID = sub_title["ID"],
                  A = sub_title["A"],
                  B = sub_title["B"],
                  vs  = sub_title["vs"])
  )

if(!is.null(isolate(input$local_statistics))){
  if(isolate(input$local_statistics)=="3"){   ### in case of a F-test for local statistics
    sub_title <- ""
  }
}

###1083:1614-1641
to_ggplot_display <- 
  ggplot(data = to_ggplot_display[1:min(lengthOfTable,40), ],
         aes(
           y = -log(as.numeric(`FDR q-value`), base = 10),
           x = reorder(Discription, -log(as.numeric(
             `FDR q-value`
           ) , base =  10))
         )) +
  geom_hline(yintercept = -log(0.05, base = 10),
             alpha = .3) +
  geom_bar(stat = 'identity', aes(fill = color), width = .70) +
  scale_fill_gradient2(
    limits = c(0, 6),
    midpoint = -log(0.05, base = 10),
    low = "grey",
    mid = "grey",
    high = "red"
  ) +
  labs(title = main_title)+
  xlab("") +
  ylab("-LOG(FDR q-value)") +
  guides(fill = "none") +
  coord_flip() +
  theme_pander()

###1083:1644-1796

  resultTable <-
    GenTable(ONTdata,  'p-value' = resultFis, topNodes = 2000)
  

if (exists("resultFis")) {
  
  resultTable$`p-value` <- gsub("< ","",resultTable$`p-value`)    ### '< 1e-30' cannot be understood
  
  resultTable$'FDR q-value' <-
    format(p.adjust(as.numeric(resultTable$`p-value`), "fdr"), digits = 3)
  nrOfSignNodes <-
    max(c(sum(
      as.numeric(resultTable$'FDR q-value') < 0.05
    ), 5))   ## with a minimum of 5
  
  ### score(resultFis) >>> 0 scores are not tolerated
  resultFis_score <- score(resultFis)
  resultFis_score[resultFis_score==0] <- 0.001
  
  #nr of nodes by FDR qvalue
  enrichmentGraph <-
    showSigOfNodes(
      ONTdata,
      resultFis_score, #k, #abs(resultFis_score), #score(resultFis),
      firstSigNodes = nrOfSignNodes,
      useInfo = 'all',
      swPlot = FALSE
    )
  
  enrichment_iGraph <-
    graph_from_graphnel(enrichmentGraph$dag)
  
  nodes <-
    data.frame(id = unique(c(
      get.data.frame(enrichment_iGraph)$from,
      get.data.frame(enrichment_iGraph)$to
    )))
  nodes$label <-
    resultTable$Term[match(nodes$id, resultTable$TERM.ID)]
  nodes$shape <- "dot"
  
  #val_col <- colorRampPalette(c("#BEBEBE","#BEBEBE","#E8AD43","#FFC100","#FFFF00","red"))(100)
  val_col <-
    colorRampPalette(c("#BEBEBE", "#BEBEBE", "yellow", "orange", "red"))(100)
  logPs <-
    -log(as.numeric(resultTable$`p-value`[match(nodes$id, resultTable$TERM.ID)]), 10)
  logPs[logPs < 1] <- 1
  logPs[logPs > 10] <- 10
  nodes$color <- val_col[logPs * 10]
  
  nodes$color[grep("CAT", nodes$id)] <- "#EEEEEE"
  nodes$shape[grep("CAT", nodes$id)] <- "box"
  nodes$color[grep("all", nodes$id)] <- "black"
  nodes$shape[grep("all", nodes$id)] <- "triangle"
  
  ## size as function of nr of annotions
  nodes$size <-
    0.2 * as.numeric(resultTable$Annotated[match(nodes$id, resultTable$TERM.ID)])
  ##
  
  nodes$size[nodes$size < 4] <- 4
  nodes$size[grep("all", nodes$id)] <- 20
  nodes$label[grep("all", nodes$id)] <- "ontology root"
  
  ## get levels
  ONT_DAG <- ONTdata@graph
  ONT_DAG <- buildLevels(ONT_DAG)
  DAG.env <- ONT_DAG$nodes2level
  DAG.env$`LION:0080986`
  nodes$level <-   sapply(nodes$id, function(id) {
    DAG.env[[id]]
  })
  
  edges <- get.data.frame(enrichment_iGraph)
  edges$color <- "grey"
  edges$arrows <- "from"
  
  ####  limiting by LION-term selection
  if(isolate(input$LIONselection)){    ### switch for LION-selection
    #TermsOfInterest ## vector with LION-names from previous section
    nodes <- nodes[nodes$label %in% c(TermsOfInterest,"ontology root"),]
    ## delete nodes not selected for analysis
  }
  ###
  
  ### delete dead ends of CATs w/o LION-terms
  CATs <-
    edges$from[paste(gsub(":\\d+", "", edges$from),
                     gsub(":\\d+", "", edges$to),
                     sep = "") %in% c("CATLION", "CATCAT")]
  CATs <- nodes$id[!(nodes$id %in% CATs)]
  CATs <- CATs[grepl("CAT", CATs)]
  
  nodes <- nodes[!nodes$id %in% CATs, ]
  
  to_display_network <- visNetwork(nodes, edges) %>%
    visPhysics(stabilization = TRUE) %>%
    visHierarchicalLayout(direction = "LR", levelSeparation = 250) %>%
    visInteraction(navigationButtons = TRUE)   %>% 
    visOptions(highlightNearest = TRUE) %>%
    visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(gravitationalConstant = -20)
    )
  
  #incProgress(.9, detail = paste("network"))
  
  ### return as list object
  
}

if(isolate(input$EmailMissingLipids)){   ## email missing lipids when option is set
  #sendEmail(subject = "missing annotations",mail_message = paste(subset(lipidExistance,`LION ID`=="not found")$input, collapse = "\n"))
  sendEmail(subject = "missing annotations", mail_message = paste(apply(lipidExistance, 1, function(row){paste(row, collapse = "\t")}), collapse = "\n"))
}



### now LUT is available, match LION names correctly
lipidExistance$'LION name' <- LUT$Discription[match(lipidExistance$'LION ID', LUT$ID)]
lipidExistance$'LION name'[is.na(lipidExistance$'LION name')] <- "not found"



list(
  ONTdata = ONTdata,
  lipidExistance = lipidExistance,
  lipidExistance_toview = lipidExistance_toview,
  matching_statistics = matching_statistics,
  to_display = to_display,
  to_display_detailed = to_display_detailed,
  lipidReport = lipidReport,
  lipidsInTerms = lipidsInTerms,
  to_ggplot_display = to_ggplot_display,
  edges_nodes = list(edges, nodes),
  to_display_network = to_display_network
  
)

###1819-2064
myoutput$value <- lipidExistance
myoutput$ontology.table <- to_display
myoutput$downloadDetailTable <- list(
  to_display_detailed = to_display_detailed,
  lipidReport = lipidReport, 
  lipidsInTerms = lipidsInTerms
)
myoutput$ontology.graph <- to_ggplot_display
myoutput$ontology.plot <- to_display_network

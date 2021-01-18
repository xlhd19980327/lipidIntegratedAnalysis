metEnrichFunc <- function(csvfile, filename, fileLoc){
  SetCandidate_cus <- function (mSetObj = NA, query_nm, can_nm) 
  {
    query_inx <- which(mSetObj$name.map$query.vec == query_nm)
    can_mat <- mSetObj$dataSet$candidates
    if (!is.null(can_mat)) {
      cmpd.db <- qs::qread("compound_db.qs")
      can_inx <- which(can_mat[, 2] == can_nm)
      if (can_inx <= nrow(can_mat)) {
        can_inx <- which(cmpd.db$name == can_nm)[1]
        hit <- cmpd.db[can_inx, , drop = F]
        mSetObj$name.map$hit.inx[query_inx] <- can_inx
        mSetObj$name.map$hit.values[query_inx] <- hit[, 
                                                      2]
        mSetObj$name.map$match.state[query_inx] <- 1
        csv.res <- mSetObj$dataSet$map.table
        if (ncol(csv.res) > 7) {
          csv.res[query_inx, ] <- c(csv.res[query_inx, 
                                            1], mSetObj$name.map$hit.values[query_inx], 
                                    hit$hmdb_id, hit$pubchem_id, hit$chebi_id, 
                                    hit$kegg_id, hit$metlin_id, hit$smiles, 1)
        }
        else {
          csv.res[query_inx, ] <- c(csv.res[query_inx, 
                                            1], mSetObj$name.map$hit.values[query_inx], 
                                    hit$hmdb_id, hit$pubchem_id, hit$kegg_id, 
                                    hit$smiles, 1)
        }
        write.csv(csv.res, file = "name_map.csv", row.names = F)
        mSetObj$dataSet$map.table <- csv.res
      }
      else {
        mSetObj$name.map$hit.inx[query_inx] <- 0
        mSetObj$name.map$hit.values[query_inx] <- ""
        mSetObj$name.map$match.state[query_inx] <- 0
        print("No name matches found.")
      }
    }
    return(mSetObj)
  }
  
  CalculateHyperScore_cus <- function (mSetObj = NA) 
  {
    mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
    nm.map <- MetaboAnalystR:::GetFinalNameMap(mSetObj)
    valid.inx <- !(is.na(nm.map$hmdb) | duplicated(nm.map$hmdb))
    ora.vec <- nm.map$hmdb[valid.inx]
    q.size <- length(ora.vec)
    if (is.na(ora.vec) || q.size == 0) {
      MetaboAnalystR:::AddErrMsg("No valid HMDB compound names found!")
      return(0)
    }
    #if (!.on.public.web & grepl("kegg", mSetObj$analSet$msetlibname)) {
    mSetObj$api$oraVec <- ora.vec
    if (mSetObj$api$filter) {
      mSetObj$api$filterData <- mSetObj$dataSet$metabo.filter.kegg
      toSend <- list(libNm = mSetObj$api$libname, filter = mSetObj$api$filter, 
                     oraVec = mSetObj$api$oraVec, filterData = mSetObj$api$filterData, 
                     excludeNum = mSetObj$api$excludeNum)
    }
    else {
      toSend <- list(libNm = mSetObj$api$libname, filter = mSetObj$api$filter, 
                     oraVec = mSetObj$api$oraVec, excludeNum = mSetObj$api$excludeNum)
    }
    MetaboAnalystR:::load_httr()
    base <- api.base
    endpoint <- "/enrichmentora"
    call <- paste(base, endpoint, sep = "")
    query_results <- httr::POST(call, body = toSend, encode = "json")
    query_results_text <- content(query_results, "text")
    query_results_json <- RJSONIO::fromJSON(query_results_text, 
                                            flatten = TRUE)
    if (is.null(query_results_json$enrichRes)) {
      MetaboAnalystR:::AddErrMsg("Error! Enrichment ORA via api.metaboanalyst.ca unsuccessful!")
      return(0)
    }
    oraDataRes <- do.call(rbind.data.frame, query_results_json$enrichRes)
    colnames(oraDataRes) <- query_results_json$enrichResColNms
    rownames(oraDataRes) <- query_results_json$enrichResRowNms
    MetaboAnalystR:::fast.write.csv(oraDataRes, 
                                    file = paste0(fileLoc, filename, "_msea_ora_result.csv"))
    mSetObj$analSet$ora.mat <- oraDataRes
    mSetObj$api$guestName <- query_results_json$guestName
    return(MetaboAnalystR:::.set.mSet(mSetObj))
    #}
    current.mset <- current.msetlib$member
    if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)) {
      current.mset <- lapply(current.mset, function(x) {
        x[x %in% mSetObj$dataSet$metabo.filter.hmdb]
      })
      mSetObj$dataSet$filtered.mset <- current.mset
      ora.vec <- ora.vec[ora.vec %in% unique(unlist(current.mset, 
                                                    use.names = FALSE))]
      q.size <- length(ora.vec)
    }
    uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)))
    set.size <- length(current.mset)
    if (set.size == 1) {
      AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!")
      return(0)
    }
    hits <- lapply(current.mset, function(x) {
      x[x %in% ora.vec]
    })
    hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE)
    if (sum(hit.num > 0) == 0) {
      AddErrMsg("No match was found to the selected metabolite set library!")
      return(0)
    }
    set.num <- unlist(lapply(current.mset, length), use.names = FALSE)
    res.mat <- matrix(NA, nrow = set.size, ncol = 6)
    rownames(res.mat) <- names(current.mset)
    colnames(res.mat) <- c("total", "expected", "hits", "Raw p", 
                           "Holm p", "FDR")
    for (i in 1:set.size) {
      res.mat[i, 1] <- set.num[i]
      res.mat[i, 2] <- q.size * (set.num[i]/uniq.count)
      res.mat[i, 3] <- hit.num[i]
      res.mat[i, 4] <- phyper(hit.num[i] - 1, set.num[i], 
                              uniq.count - set.num[i], q.size, lower.tail = F)
    }
    res.mat[, 5] <- p.adjust(res.mat[, 4], "holm")
    res.mat[, 6] <- p.adjust(res.mat[, 4], "fdr")
    res.mat <- res.mat[hit.num > 0, ]
    ord.inx <- order(res.mat[, 4])
    mSetObj$analSet$ora.mat <- signif(res.mat[ord.inx, ], 3)
    mSetObj$analSet$ora.hits <- hits
    fast.write.csv(mSetObj$analSet$ora.mat, file = "msea_ora_result.csv")
    return(.set.mSet(mSetObj))
  }
  
  cmpd.vec <- read.csv(csvfile)$x
  ## cmpd.vec will be the input
  
  if(length(cmpd.vec) == 0){
    print(paste0("No significant ", filename, " regulation features found!"))
    return(1)
  }else{
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
        mSet2<-SetCandidate_cus(mSet2, i, 
                                ifelse(mSet2[["name.map"]][["hits.candidate.list"]][1, 1] != '', 
                                       mSet2[["name.map"]][["hits.candidate.list"]][1, 1], "NA"));
      }
    }
    ## Enrichment
    mSet2<-SetMetabolomeFilter(mSet2, F);
    #kegg(may add other enrichment options)
    mSet2<-SetCurrentMsetLib(mSet2, "kegg_pathway", 2);
    #output a chart
    mSet2<-CalculateHyperScore_cus(mSet2)
    #bar plot
    mSet2<-PlotORA(mSet2, paste0(fileLoc, filename, "_ora_"), "net", "png", 72, width=NA)
    #dot plot
    mSet2<-PlotEnrichDotPlot(mSet2, "ora", paste0(fileLoc, filename, "_ora_dot_"), "png", 72, width=NA)
    return(0)
  }
    
}

statFAChains_pathAna <- function(lipid_subclass_tidyStat, fileLoc, lipsample, spe){
  source("./utilityFunc/addLibnames.R")
  #source("./utilityFunc/FAchainStat.R")
  
  ## Data preparation
  # Use the following code to get lipid_subclass_tidyStat
  #lipid_subclass_integStat <- FAchainStat(..., stat = F)
  
  # Lipid reaction library file input
  if(spe == "hsa"){
    #Just need contain these colums: "reaction_id", "gene_symbol"
    allReact <- read.csv("./patternHunting/lipidReatomeStat/hsaDB/hsa_all_integData.csv")
    
    lipReact <- read.csv("./patternHunting/lipidReatomeStat/hsaDB/hsa_lipidreact.csv")
  }
  if(spe == "mmu"){
    #Just need contain these colums: "reaction_id", "gene_symbol"
    allReact <- read.csv("./patternHunting/lipidReatomeStat/mmuDB/mmu_integInfo_modified.csv")
    
    lipReact <- read.csv("./patternHunting/lipidReatomeStat/mmuDB/mmu_lipidreact.csv")
  }
  #use unique/distinct, other info will not consider here
  allReact <- unique(subset(allReact, select = c("reaction_id", "gene_symbol")))
  integReact <- inner_join(lipReact, allReact, by = "reaction_id")
  integReact <- subset(integReact, 
                       select = c("reaction_id", "firstLev", "gene_symbol", "R_lipid", "P_lipid"))
  
  # !!!!!WARNING: Some libnames will be wrong, cause carbon base info add
  lipidsubclass <- cbind(lipid_subclass_tidyStat, 
                         libnames = apply(lipid_subclass_tidyStat, 1, addLibnames, colnames(lipid_subclass_tidyStat)))
  lipsubname <- levels(factor(lipidsubclass$libnames))
  lipsubname <- unlist(strsplit(lipsubname, ";"))
  lipsubname <- levels(factor(lipsubname))
  # Only use the reactions that lipid Class both corresponding to Reactant and Product
  ind <- integReact$R_lipid %in% lipsubname & integReact$P_lipid %in% lipsubname
  integReact_filter <- integReact[ind, ]
  # Calculate library lipid names ind from data 
  cal_libInd <- function(x, y){
    libname <- unlist(strsplit(x, ";"))
    getInd <- function(i, lib = y){
      return(which(lib == i))
    }
    ind <- lapply(libname, getInd)
    ind <- unlist(ind)
    return(ind)
  }
  R_ind <- lapply(lipidsubclass$libnames, cal_libInd, integReact_filter$R_lipid)
  P_ind <- lapply(lipidsubclass$libnames, cal_libInd, integReact_filter$P_lipid)
  # Transfer and integrate library lipid names ind to data ind
  integ_ind <- function(libInd, type = c("R_dataInd", "P_dataInd")){
    n <- length(libInd)
    result <- data.frame()
    for(i in 1:n){
      dataIndi <- i
      libIndi <- libInd[[i]]
      result <- rbind(result, 
                      data.frame(libInd = libIndi, rep(dataIndi, length(libIndi))))
    }
    colnames(result)[2] <- type
    return(result)
  }
  R_ind <- integ_ind(R_ind, "R_dataInd")
  P_ind <- integ_ind(P_ind, "P_dataInd")
  reaction_ind <- inner_join(R_ind, P_ind)
  reaction_ind_list <- split(reaction_ind, reaction_ind$libInd)
  getReactionInfo <- function(ind, data = lipidsubclass){
    getInd <- function(x, y = ind){
      res <- y[[x]][!is.na(y[[x]])]
      res <- unique(res)
      return(res)
    }
    R_lipid_ind <- getInd("R_dataInd")
    R_lipid <- data[R_lipid_ind, ]
    P_lipid_ind <- getInd("P_dataInd")
    P_lipid <- data[P_lipid_ind, ]
    return(list(
      R_lipid = R_lipid, 
      P_lipid = P_lipid
    ))
  }
  reaction_info <- lapply(reaction_ind_list, getReactionInfo)
  gene_info <- data.frame(
    pathway = integReact_filter[as.numeric(names(reaction_info)), ]$firstLev,
    libInd = names(reaction_info), 
    gene = integReact_filter[as.numeric(names(reaction_info)), ]$gene_symbol
  )
  libIndList <- split(gene_info$libInd, gene_info$pathway)
  bindGeneList <- function(ind){
    ind <-  names(reaction_info) %in% ind
    data <- reaction_info[ind]
    R_lipid <- data.frame()
    P_lipid <- data.frame()
    for(i in data){
      R_lipid <- rbind(R_lipid, data.frame(i[["R_lipid"]]))
      P_lipid <- rbind(P_lipid, data.frame(i[["P_lipid"]]))
    }
    return(list(
      R_lipid = unique(R_lipid), 
      P_lipid = unique(P_lipid)
    ))
  }
  reaction_info <- lapply(libIndList, bindGeneList)
  
  ## Statistics use our methods, use t-test to test the statistics(log2FC difference)
  getRegState <- function(data){
    R_info <- subset(data[["R_lipid"]], subset = group == lipsample)
    P_info <- subset(data[["P_lipid"]], subset = group == lipsample)
    getInfo <- function(info, type = c("R", "P")){
      n <- length(info$log2FC)
      up <- sum(info$log2FC > log2(1.5))
      down <- sum(info$log2FC < -log2(1.5))
      zero <- n - up - down 
      log2FC_avg <- sum(info$log2FC) / n
      up_pct <- up / n
      down_pct <- down / n
      if(up > down){
        result <- "up"
      } else if(up == down){
        result <- "zero"
      } else{
        result <- "down"
      }
      data_result <- data.frame(log2FC_avg, n, up_pct, down_pct)
      colnames(data_result) <- paste0(type, "_", colnames(data_result))
      return(data_result)
    }
    R_result <- getInfo(R_info, type = "R")
    P_result <- getInfo(P_info, type = "P")
    P_R <- P_result$P_log2FC_avg - R_result$R_log2FC_avg
    p_value <- ifelse(P_result$P_n >= 3 & R_result$R_n >= 3, 
                      t.test(R_info$log2FC, P_info$log2FC)$p.value, 
                      NA)
    regState <- cbind(R_result, P_result, P_R, p_value = p_value)
    return(regState)
  }
  regState <- lapply(reaction_info, getRegState)
  regState <- cbind(pathway = gsub("(.*)\\(.*?\\)$", "\\1", names(regState)), 
                    do.call(rbind, regState))
  write.csv(regState, 
            paste0(fileLoc, "regStat_", lipsample, ".csv"), row.names = F)
  gene_info$pathway <- gsub("(.*)\\(.*?\\)$", "\\1", gene_info$pathway)
  return(list(regState = regState, path_info = gene_info))
}

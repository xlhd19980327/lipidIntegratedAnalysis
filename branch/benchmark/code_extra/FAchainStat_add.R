## plotInfo == "FA_info": Statistics fatty acid chain&unsaturated info in that lipid class
## plotInfo == "all_info": Statistics all chain&unsaturated info in that lipid class
FAchainStat <- function(dataSet, mSet,  
                        fileLoc, plotInfo, 
                        #use this for our statistics method, client cannot modify
                        stat = F){
  allgroups <- dataSet$allgroups
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  dataType <- dataSet$dataType
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  
  ## Source will offer the following contents:
  ## Function(s): getFAsInfo
  if(dataType == "LipidSearch"){
    source("./utilityFunc/getFAsInfo.R")
  }
  if(dataType == "MS_DIAL"){
    source("./utilityFunc/getFAsInfo_msdial.R")
  }
  if(dataType == "HCC"){
    source("./branch/benchmark/code_extra/getFAInfo_HCC.R")
  }
  ## Source will offer the following contents:
  ## Function(s): getPValue
  source("./utilityFunc/getPValue.R")
  ## Source will offer the following contents:
  ## Function(s): addSigLabel
  source("./utilityFunc/addSigLabel.R")
  
  data_tidy <- as.data.frame(t(mSet[["dataSet"]][["preproc"]])) %>%
    rownames_to_column(var = "lipidName") %>%
    mutate(Class = switch(dataSet$dataType,
                          LipidSearch = gsub("(.*?)\\(.*", "\\1", lipidName), 
                          MS_DIAL = gsub("(.*?) .*", "\\1", lipidName), 
                          HCC = gsub("(.*?)\\(.*", "\\1", lipidName)))

  ## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
  lipids <- data_tidy$lipidName
  # Get the FA info and MS1 info
  fasInfo <- lapply(lipids, getFAsInfo)
  ms1Info <- sapply(fasInfo, function(x) x$ms1)
  data_tidy <- cbind(data_tidy, ms1 = ms1Info)
  
  ### Tidy for and do Visualization ###
  # Tidy and integrate itensity of lipid class containing FA chain info
  lipid_subclass_handle <- data.frame()
  if(plotInfo == "all_info"){
    # Use MS1 & MS2 to do the later statistics
    for(i in 1:nrow(data_tidy)){
      # e,p connection info will be ignored
      fa <- fasInfo[[i]][[1]]
      Class <- names(fasInfo[[i]][1])
      if(any(fa %in% "UnknownPattern")){
        subclass <- paste0(Class, "(", fa, ")")
      }else{
        chains <- sum(as.numeric(gsub(".*?([0-9]+):.*", "\\1", fa)))
        unsaturate <- sum(as.numeric(gsub(".*?:([0-9]+).*", "\\1", fa)))
        subclass <- paste0(Class, "(", chains, ":", unsaturate, ")")
      }
      
      lipid_subclass_handle <- rbind(lipid_subclass_handle, 
                                     cbind(subclass = subclass, 
                                           bind_rows(replicate(length(subclass), data_tidy[i, ], simplify = FALSE))))
    }
  }else if(plotInfo == "FA_info"){
    # Use MS2 to do the later statistics only
    for(i in 1:nrow(data_tidy)){
      # e,p connection info will be ignored
      fa <- fasInfo[[i]][[1]]
      Class <- names(fasInfo[[i]][1])
      subclass <- paste0(Class, "(", fa, ")")
      lipid_subclass_handle <- rbind(lipid_subclass_handle, 
                                     cbind(subclass = subclass, 
                                           bind_rows(replicate(length(subclass), data_tidy[i, ], simplify = FALSE))))
    }
    lipid_subclass_handle <- subset(lipid_subclass_handle, subset = ms1 == F)
  }
  
  ## Use lipid_subclass_handle to calculate itensity of lipid class containing FA chain info(subclass)
  ## Visualize with the facet plot
  divnum <- function(x){
    ##!!!!!WARNING: Class statistics will only contain the lipid class in getFAsInfo function recognize
    Class <- unique(x)
    if(Class %in% c("Cer", "SM", "SPH", "FA", "MG", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS", "ChE")){
      div <- 1
    }else if(Class %in% c("DG", "PA", "PC", "PE", "PG", "PI", "PS")){
      div <- 2
    }else if(Class %in% c("TG")){
      div <- 3
    }else if(Class %in% c("CL")){
      div <- 4
    }
    return(div)
  }
  lipid_subclass_stat <- lipid_subclass_handle %>%
    ungroup() %>%
    select(-ms1, -lipidName) %>%
    filter(!grepl("UnknownPattern", subclass)) %>%
    gather(key = "case", value = "lipidsum", -subclass, -Class) %>%
    group_by(subclass, case, Class) %>% 
    filter(!is.na(lipidsum)) %>%    #delete low abundance of the lipid signal
    #!!!!!WARNING: this algorithm may not appropriate
    summarise(lipidsum = switch(plotInfo, 
                                FA_info = sum(lipidsum) / divnum(Class), 
                                all_info = sum(lipidsum))) %>% #calc abundance of the lipid in a sample
    mutate(group = allgroups[match(case, names(allgroups))]) %>%
    ungroup() %>%
    group_by(group, subclass) %>%
    filter(n() >= 3) #if the lipid signal of one group are detected in no more than 2 samples, the lipid will be dropped
  lipid_subclass_stat_output <- lipid_subclass_handle %>%
    ungroup() %>%
    select(-ms1, -lipidName) %>%
    filter(!grepl("UnknownPattern", subclass)) %>%
    gather(key = "case", value = "lipidsum", -subclass, -Class) %>%
    group_by(subclass, case) %>% 
    filter(!is.na(lipidsum)) %>%    #delete low abundance of the lipid signal
    #!!!!!WARNING: this algorithm may not appropriate
    summarise(lipidsum = switch(plotInfo, 
                                FA_info = sum(lipidsum) / divnum(Class), 
                                all_info = sum(lipidsum)))  %>%
    spread(key = case, value = lipidsum)
  lipid_subclass_stat2 <- lipid_subclass_stat %>%
    group_by(Class, subclass, group) %>%
    summarise(realmean = mean(lipidsum) ,
              sd = sd(lipidsum))
  lipid_subclass_stat_p <- split(lipid_subclass_stat, lipid_subclass_stat$subclass)
  lipid_subclass_stat3 <- lapply(lipid_subclass_stat_p, getPValue, "FAchain", controlGrp)
  lipid_subclass_stat3 <- do.call(rbind, lipid_subclass_stat3)
  lipid_subclass_integStat <- left_join(lipid_subclass_stat2, lipid_subclass_stat3)
  sigLabel2 <- addSigLabel(lipid_subclass_integStat$p)
  lipid_subclass_integStat <- cbind(lipid_subclass_integStat, sigLabel = sigLabel2)
  if(stat == T){
    #Use for "statFAChains"
    lipid_subclass_tidyStat <- lipid_subclass_integStat %>%
      ungroup() %>%
      select(subclass, Class, group, realmean) %>%
      distinct() %>%
      group_by(subclass) %>%
      mutate(ind = ifelse(group == controlGrp, "A", "B")) %>% 
      arrange(ind, .by_group = T) %>% ## put controlGrp into ind[1] so that we can recognize the control group
      mutate(log2FC = log2(realmean / realmean[1])) %>% 
      filter(group != controlGrp) %>%
      select(-ind)
    return(lipid_subclass_tidyStat)
  }
  write.csv(lipid_subclass_stat_output, 
            paste0(fileLoc, "lipid_subclass_stat_", plotInfo, "_", pname, ".csv"), 
            row.names = F)  
  for(i in unique(lipid_subclass_integStat$Class)){
    oneLipClassData <- subset(lipid_subclass_integStat, 
                              subset = Class == i)
    oneLipClassData2 <- subset(lipid_subclass_stat, 
                               subset = Class == i)
    oneLipClassData$group <- 
      factor(oneLipClassData$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
    oneLipClassData2$group <- 
      factor(oneLipClassData2$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
    if(nrow(oneLipClassData) != 0 & nrow(oneLipClassData2) != 0){
      nsubclass <- length(unique(oneLipClassData$subclass))
      h_plot <- ifelse(nsubclass > 3, nsubclass, 3) / 3 * (25/12)
      ggplot() +
        geom_bar(data = oneLipClassData, aes(x = group, y = realmean, fill = group), stat = "identity") +
        geom_errorbar(data = oneLipClassData, aes(x = group, ymin = realmean, ymax = realmean + sd, color = group), width = 0.2) +
        geom_dotplot(data = oneLipClassData2, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
        geom_text(data = oneLipClassData, aes(x = group, y = realmean+2*sd, label = sigLabel),
                  size = 3, fontface = "bold", color = "red") +
        theme_bw() +
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank()) +
        facet_wrap(~subclass, scales="free", ncol = 3) +
        scale_fill_aaas() +
        scale_color_aaas() +
        labs(x = "group",
             y = "concentration", 
             color = "group", 
             title = "Lipid subclass statistics") +
        theme(plot.title = element_text(hjust = 0.5, size = 20)) 
      ggsave(paste0(fileLoc, "integPlot_", i, "_", pname, ".pdf"), 
             dpi = 300, width = 9, height = h_plot, limitsize = FALSE)
    }
  }
  
  ## Visualize with the tile plot
  realmean <- lipid_subclass_stat2$realmean
  names(realmean) <- lipid_subclass_stat2$group
  getFC <- function(x){
    cpr <- x[names(x) == controlGrp]
    if(length(cpr) != 0){
      res <- x/cpr
    }else{
      res <- rep(NA, length(x))
    }
  }
  FC <- unlist(tapply(realmean, lipid_subclass_stat2$subclass, getFC))
  lipid_subclass_stat_tile <- lipid_subclass_stat2 %>%
    left_join(lipid_subclass_stat3) %>%
    ungroup() %>%
    add_column(FC = round(FC, 2)) %>%
    #!!!!!WARNING: drop subclass that have other info(i.e. d/t symbol)
    filter(!grepl("\\([a-z][0-9]+:", subclass)) %>%
    mutate(chain = as.numeric(gsub(".*?([0-9]+):.*", "\\1", subclass)), 
           unsaturate = as.numeric(gsub(".*?:([0-9]+).*", "\\1", subclass)), 
           Class = gsub("(.*?)\\(.*", "\\1", subclass))
  sigLabel3 <- addSigLabel(lipid_subclass_stat_tile$p)
  lipid_subclass_stat_tile <- cbind(lipid_subclass_stat_tile, sigLabel = sigLabel3)
  for(i in groupsLevel[groupsLevel != controlGrp]){
    oneGrpdata <- subset(lipid_subclass_stat_tile, 
                         subset = group == i & !is.na(FC))
    oneGrpdata <- oneGrpdata %>%
      mutate(regState = apply(oneGrpdata, 1, function(x){
        names(x) <- colnames(oneGrpdata)
        if(as.numeric(x[names(x) == "FC"]) > 1){
          if(x[names(x) == "sigLabel"] != ""){
            return("up-sig")
          } else{
            return("up-unsig")
          }
        } else if(as.numeric(x[names(x) == "FC"]) < 1){
          if(x[names(x) == "sigLabel"] != ""){
            return("down-sig")
          } else{
            return("down-unsig")
          }
        } else{
          return("notreg")
        }
      }))
    grid_df <- expand.grid(x = 0:max(oneGrpdata$chain), y = 0:max(oneGrpdata$unsaturate))
    if(plotInfo == "all_info"){
      w_plot <- 28
    }
    if(plotInfo == "FA_info"){
      w_plot <- 18
    }
    h_plot <- w_plot * max(oneGrpdata$unsaturate) / max(max(oneGrpdata$chain), 1) * length(unique(oneGrpdata$Class)) +6
    grid_plot <- ggplot(mapping = aes(x = chain, y = unsaturate)) + 
      geom_point(data = subset(oneGrpdata, regState != "notreg"), 
                 mapping = aes(fill = factor(regState)), color = "black", 
                 size = 9, shape = 21) + 
      geom_point(data = subset(oneGrpdata, regState == "notreg"), 
                 fill = "transparent", color = "black", 
                 size = 9, shape = 21) + 
      geom_text(data = oneGrpdata, aes(label = FC),
                size = 3, fontface = "bold", color = "black") +
      geom_tile(data = grid_df, aes(x = x, y = y), 
                fill = "transparent", color = alpha("grey50", 0.3)) +
      scale_fill_manual(values = c("up-unsig" = alpha("red", 0.3), 
                                   "down-unsig" = alpha("blue", 0.3), 
                                   "up-sig" = alpha("red", 0.6), 
                                   "down-sig" = alpha("blue", 0.6))) +
      theme(
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank(),
        axis.text.x = element_text(vjust = 6),
        axis.text.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5, size = 20)
      ) +
      guides(fill = guide_legend(
        ncol = 2, nrow = 2,
        label.position = "bottom", 
        title.hjust = 0.5
      )) +
      facet_wrap(~Class, ncol = 1, scales="free") + 
      scale_x_continuous(breaks = seq(0, max(oneGrpdata$chain), 2), 
                         expand = c(0,0)) +
      scale_y_continuous(breaks = seq(0, max(oneGrpdata$unsaturate), 1)) +
      labs(x = "Number of FA chain carbon",
           y = "Number of FA double-bonds", 
           title = "Lipid FA chains statistics", 
           fill = "Regulation State") 
    ggsave(paste0(fileLoc, "tilePlot_", i, "_", plotInfo, ".pdf"), 
           grid_plot, 
           dpi = 300, width = w_plot, height = h_plot, limitsize = FALSE)
  }
}

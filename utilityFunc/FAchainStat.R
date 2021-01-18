## plotInfo == "FA_info": Statistics fatty acid chain&unsaturated info in that lipid class
## plotInfo == "all_info": Statistics all chain&unsaturated info in that lipid class
FAchainStat <- function(dataSet, mSet,  
                        fileLoc, plotInfo = "FA_info", ignore = T, topnum = 70,
                        #use this for our statistics method, client cannot modify
                        stat = F, stat2 = F){
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
  ## Source will offer the following contents:
  ## Function(s): getPValue
  source("./utilityFunc/getPValue.R")
  ## Source will offer the following contents:
  ## Function(s): addSigLabel
  source("./utilityFunc/addSigLabel.R")
  ## Source will offer the following contents:
  ## Function(s): getClassInfo
  source("./utilityFunc/getClassInfo.R")
  
  data_tidy <- dataSet[["data"]] %>%
    mutate(lipidName = dataSet[["lipidName"]], 
           Class = switch(dataSet$dataType,
                          LipidSearch = sapply(lipidName, getClassInfo, "LipidSearch", ignore = ignore), 
                          MS_DIAL = sapply(lipidName, getClassInfo, "MS_DIAL", ignore = ignore)))
  
  ## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
  lipids <- data_tidy$lipidName
  # Get the FA info and MS1 info
  fasInfo <- lapply(lipids, getFAsInfo, ignore)
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
    class1 <- c("FA", "MG", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS", "ChE")
    epclass1 <- c(paste0(class1, "(O)"), paste0(class1, "(P)"))
    class2 <- c("DG", "PA", "PC", "PE", "PG", "PI", "PS")
    epclass2 <- c(paste0(class2, "(O)"), paste0(class2, "(P)"))
    Class <- unique(x)
    if(Class %in% c("Cer", "SM", "SPH", class1, epclass1)){
      div <- 1
    }else if(Class %in% c(class2, epclass2)){
      div <- 2
    }else if(Class %in% c("TG", "TG(O)", "TG(P)")){
      div <- 3
    }else if(Class %in% c("CL", "CL(O)", "CL(P)")){
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
  
  ## Visualize with the heatmap plot
  datagroup <- factor(allgroups, 
                      levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
  colorpars <- plottingPalettes(n = length(groupsLevel), type = "discrete")
  names(colorpars) <- c(controlGrp, groupsLevel[groupsLevel != controlGrp])
  
  lipid_subclass_heatmap <- lipid_subclass_stat_output %>%
    mutate(Class = gsub("\\(.*\\)$", "", subclass)) %>%
    column_to_rownames(var = "subclass")
  lipid_subclass_heatmap <- lipid_subclass_heatmap[apply(lipid_subclass_heatmap, 1, function(x) ifelse(is.na(sd(x[-length(x)], na.rm = T)), F, ifelse(sd(x[-length(x)], na.rm = T) == 0, F, T))), ]
  if(!is.na(topnum)){
    if(topnum > nrow(lipid_subclass_heatmap)){
      cat("Not enough features found! Using all the lipid items.\n")
      topnum <- nrow(lipid_subclass_heatmap)
      pname2 <- "alllipids"
    }else{
      pname2 <- paste0("top", topnum)
    }
    data_mad <- apply(lipid_subclass_heatmap, 1, function(x) mad(as.numeric(x[-length(x)]), na.rm = T))
    lipid_subclass_heatmap <- lipid_subclass_heatmap[order(data_mad, decreasing=T), ][1:topnum, ]
  }else{
    pname2 <- "alllipids"
  }
  # if(nrow(data_Class_det) < nrow(data_tidy)){
  #   cat("Some lipids in low variation will not show in the plot.\n")
  # }
  lipid_subclass_heatmap <- arrange(lipid_subclass_heatmap, Class) 
  ClassInfo <- lipid_subclass_heatmap$Class
  colorpars2 <- plottingPalettes(n = length(unique(ClassInfo))+length(groupsLevel), type = "discrete")[-1:-length(groupsLevel)]
  names(colorpars2) <- unique(ClassInfo)
  lipid_subclass_heatmap <- select(lipid_subclass_heatmap, -Class) 
  lipid_subclass_heatmap <- lipid_subclass_heatmap[, match(colnames(lipid_subclass_heatmap), names(allgroups))]
  x <- pheatmap::pheatmap(mat = lipid_subclass_heatmap,
                          annotation_col = data.frame(group = datagroup, row.names = colnames(lipid_subclass_heatmap)), 
                          annotation_row = data.frame(lipidClass = ClassInfo, row.names = rownames(lipid_subclass_heatmap)),
                          # fontsize_col = 20, 
                          # fontsize_row = 20, 
                          fontsize = 20,
                          clustering_distance_rows = "euclidean", 
                          clustering_distance_cols = "euclidean", 
                          clustering_methods = "ward.D", 
                          cluster_rows = F, 
                          cluster_cols = F, 
                          scale = "row", 
                          show_rownames = T, 
                          annotation_colors = list(group = colorpars, lipidClass = colorpars2))
  w_heat <- 12/7.65*(4.95/9*length(allgroups)+2.7)
  h_heat <- 28/17.8*(16.7/70*nrow(lipid_subclass_heatmap)+1.1)
  pdf(paste0(fileLoc, "heatmap_lipsubClass_", pname, "_", pname2, ".pdf"), width=w_heat, height=h_heat)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
  
  if(stat2 == T){
    return(lipid_subclass_stat_output)
  }
  lipid_subclass_stat2 <- lipid_subclass_stat %>%
    group_by(Class, subclass, group) %>%
    summarise(realmean = mean(lipidsum) ,
              sd = sd(lipidsum))
  lipid_subclass_stat_p <- split(lipid_subclass_stat, lipid_subclass_stat$subclass)
  lipid_subclass_stat3 <- lapply(lipid_subclass_stat_p, getPValue, "FAchain", controlGrp)
  lipid_subclass_stat3 <- do.call(rbind, lipid_subclass_stat3)
  #lipid_subclass_integStat <- left_join(lipid_subclass_stat2, lipid_subclass_stat3)
  #sigLabel2 <- addSigLabel(lipid_subclass_integStat$p)
  #lipid_subclass_integStat <- cbind(lipid_subclass_integStat, sigLabel = sigLabel2)
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
  
  ## Visualize with the subclass plot
  data_sub_classSum_stat_split <- split(lipid_subclass_stat, f = lipid_subclass_stat$Class)
  seq_gen <<-function(x, y) {
    myi <<- myi + 1L
    return(list(p.value = all_pvalues[myi]))
  }
  get_multi_pvalue <- function(x){
    ## Dunnett’s multiple comparison test in one-way ANOVA
    if(length(unique(x$group)) == 1){
      return(rep(NA, length(groupsLevel)-1))
    }
    data_ano <- aov(lipidsum ~ group, data = x)
    multcp <- glht(data_ano, linfct=mcp(group="Dunnett"),alternative="two.side") 
    multcp_sum <- summary(multcp)
    return(multcp_sum$test$pvalues)
  }
  for(i in 1:length(data_sub_classSum_stat_split)){
    all_pvalues <- c()
    cplist <- mapply(c, controlGrp, groupsLevel[groupsLevel != controlGrp], 
                     SIMPLIFY = F, USE.NAMES = F)
    data_sub_classSum_stat_split[[i]]$group <- 
      factor(data_sub_classSum_stat_split[[i]]$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
    Classi <- data_sub_classSum_stat_split[[i]]$Class[1]
    nclass <- length(unique(data_sub_classSum_stat_split[[i]]$subclass))
    
    w_plot <- 22/14.0*(1.5+((12.5/3-0.8)/3*length(groupsLevel)+0.8)*ifelse(nclass >= 3, 3, nclass))
    h_plot <- 100/29.4*(29.15/13*ifelse(nclass%%3, nclass%/%3+1, nclass%/%3)+0.25)
    
    p <- ggplot(data = data_sub_classSum_stat_split[[i]], aes(x = group, y = lipidsum, color = group)) + 
      geom_boxplot(width = 0.5, size = 1.2)+ 
      #geom_jitter(color = "black", position=position_jitter(0.2)) +
      scale_color_npg() +
      facet_wrap(~subclass, scales="free", ncol = 3)+
      labs(title = Classi, 
           x = "group",
           y = "total concentration") +
      theme_classic() +
      theme(
        strip.background = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"), 
        axis.title = element_text(size = 25), 
        axis.text = element_text(size = 25),
        legend.text = element_text(size = 25), 
        legend.title = element_text(size = 25), 
        strip.text = element_text(hjust = 0.5, size = 25, face = "bold")
      ) 
    if(length(groupsLevel) == 2){
      p <- p + ggsignif::geom_signif(comparisons = cplist, color = "black", map_signif_level = T, 
                                     step_increase = 0.1, tip_length = 0, textsize = 8, test = "t.test")
    }else if(length(groupsLevel) > 2){
      data_sub_classSum_stat_split_subspt <- split(data_sub_classSum_stat_split[[i]], f = data_sub_classSum_stat_split[[i]]$subclass)
      print(Classi)
      all_pvalues <- unlist(lapply(data_sub_classSum_stat_split_subspt, get_multi_pvalue))
      myi <<- 0
      p <- p + 
        ggsignif::geom_signif(comparisons = cplist, color = "black", map_signif_level = T, test = "seq_gen", 
                              step_increase = 0.1, tip_length = 0, textsize = 8)
    }
    ggsave(paste0(fileLoc, Classi, ".pdf"), p, width = w_plot, height = h_plot, limitsize = FALSE)
    # p <- ggboxplot(data_sub_classSum_stat_split[[i]], x = "group", y = "lipidsum",
    #                color = "group", add = "jitter")+ # Add global p-value
    #   stat_compare_means(aes(label = ..p.signif..),
    #                      method = "t.test", ref.group = controlGrp)+
    #   scale_color_npg() +
    #   facet_wrap(~subclass, scales="free", ncol = 3) +
    #   labs(title = paste0(Classi, " subclass statistics"), 
    #        x = "group",
    #        y = "concentration") +
    #   theme(legend.position = "right", 
    #         strip.background = element_blank(), 
    #         plot.title = element_text(hjust = 0.5, size = 20), 
    #         axis.title = element_text(size = 15), 
    #         legend.text = element_text(size = 12),
    #         legend.title = element_text(size = 12)
    #         # legend.text = element_blank(), 
    #         # legend.title = element_blank()
    #         ) +
    #   scale_y_continuous(expand = c(0.2,0))
    # ggsave(paste0(fileLoc, Classi, ".pdf"), p, width = w_plot, height = h_plot, limitsize = FALSE)
  }
  # w_plot <- length(groupsLevel) * (6/3)
  # for(i in unique(lipid_subclass_integStat$Class)){
  #   oneLipClassData <- subset(lipid_subclass_integStat, 
  #                             subset = Class == i)
  #   oneLipClassData2 <- subset(lipid_subclass_stat, 
  #                              subset = Class == i)
  #   oneLipClassData$group <- 
  #     factor(oneLipClassData$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  #   oneLipClassData2$group <- 
  #     factor(oneLipClassData2$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  #   if(nrow(oneLipClassData) != 0 & nrow(oneLipClassData2) != 0){
  #     nsubclass <- length(unique(oneLipClassData$subclass))
  #     h_plot <- ifelse(nsubclass > 3, nsubclass, 3) / 3 * (25/12)
  #     ggplot() +
  #       geom_bar(data = oneLipClassData, aes(x = group, y = realmean, fill = group), stat = "identity") +
  #       geom_errorbar(data = oneLipClassData, aes(x = group, ymin = realmean, ymax = realmean + sd, color = group), width = 0.2) +
  #       geom_dotplot(data = oneLipClassData2, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
  #       geom_text(data = oneLipClassData, aes(x = group, y = realmean+2*sd, label = sigLabel),
  #                 size = 3, fontface = "bold", color = "red") +
  #       theme_bw() +
  #       theme(axis.text.x = element_blank(),
  #             axis.ticks.x = element_blank()) +
  #       facet_wrap(~subclass, scales="free", ncol = 3) +
  #       scale_fill_npg() +
  #       scale_color_npg() +
  #       labs(x = "group",
  #            y = "concentration",
  #            color = "group",
  #            title = "Lipid subclass statistics") +
  #       theme(plot.title = element_text(hjust = 0.5, size = 20))
  #     ggsave(paste0(fileLoc, "integPlot_", i, "_", pname, ".pdf"), 
  #            dpi = 300, width = w_plot, height = h_plot, limitsize = FALSE)
  #   }
  # }
  
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
    left_join(lipid_subclass_stat3, copy = T) %>%
    ungroup() %>%
    add_column(FC = round(FC, 2)) %>%
    #!!!!!WARNING: drop subclass that have other info(i.e. d/t symbol)
    filter(!grepl("\\([a-z][0-9]+:", subclass)) %>%
    mutate(chain = as.numeric(gsub(".*?([0-9]+):.*", "\\1", subclass)), 
           unsaturate = as.numeric(gsub(".*?:([0-9]+).*", "\\1", subclass)), 
           Class = Class)
  sigLabel3 <- addSigLabel(lipid_subclass_stat_tile$p)
  lipid_subclass_stat_tile <- cbind(lipid_subclass_stat_tile, sigLabel = sigLabel3)
  for(i in groupsLevel[groupsLevel != controlGrp]){
    oneGrpdata <- subset(lipid_subclass_stat_tile, 
                         subset = group == i & !is.na(FC))
    oneGrpdata <- oneGrpdata %>%
      mutate(regState = apply(oneGrpdata, 1, function(x){
        names(x) <- colnames(oneGrpdata)
        if(as.numeric(x[names(x) == "FC"]) > 1){
          if(x[names(x) == "sigLabel"] != "ns"){
            return("up-sig")
          } else{
            return("up-unsig")
          }
        } else if(as.numeric(x[names(x) == "FC"]) < 1){
          if(x[names(x) == "sigLabel"] != "ns"){
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
      w2_plot <- 28
    }
    if(plotInfo == "FA_info"){
      w2_plot <- 18
    }
    h2_plot <- w2_plot * max(oneGrpdata$unsaturate) / max(max(oneGrpdata$chain), 1) * length(unique(oneGrpdata$Class)) +6
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
           dpi = 300, width = w2_plot, height = h2_plot, limitsize = FALSE)
  }
}

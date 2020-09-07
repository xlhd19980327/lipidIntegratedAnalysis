headgroupStat <- function(dataSet, mSet,  
                          fileLoc){
  allgroups <- dataSet$allgroups
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  
  ## Source will offer the following contents:
  ## Function(s): getFAsInfo
  source("./utilityFunc/getFAsInfo.R")
  ## Source will offer the following contents:
  ## Function(s): getPValue
  source("./utilityFunc/getPValue.R")
  ## Source will offer the following contents:
  ## Function(s): addSigLabel
  source("./utilityFunc/addSigLabel.R")
  
  data_tidy <- as.data.frame(t(mSet[["dataSet"]][["preproc"]])) %>%
    rownames_to_column(var = "lipidName") %>%
    mutate(Class = gsub("(.*?)\\(.*", "\\1", lipidName))
  
  ## Seperate MS1 and MS2 lipids & Calculate itensity of lipid class containing FA chain info
  lipids <- data_tidy$lipidName
  # Get the FA info and MS1 info
  fasInfo <- lapply(lipids, getFAsInfo)
  ms1Info <- sapply(fasInfo, function(x) x$ms1)
  data_tidy <- cbind(data_tidy, ms1 = ms1Info)
  data_sub_ms1 <- subset(data_tidy, subset = ms1 == T)
  data_sub_ms2 <- subset(data_tidy, subset = ms1 == F)
  # Tidy and integrate itensity of lipid class containing FA chain info
  # Use MS2 to do the later statistics only
  lipid_subclass_handle <- data.frame()
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
  
  ### Tidy for and do Visualization ###
  ## Use data_tidy(Delete duplication) to calculate itensity of each lipid class
  data_sub_classSum_stat <- data_tidy %>%
    ungroup() %>%
    select(-ms1, -lipidName) %>%
    gather(key = "case", value = "lipidsum", -Class) %>%
    group_by(Class, case) %>%
    summarise(lipidsum = sum(lipidsum, na.rm = T)) %>%
    mutate(group = allgroups[match(case, names(allgroups))]) 
  data_sub_classSum_stat2 <- data_sub_classSum_stat %>%
    group_by(Class, group) %>%
    summarise(realmean = mean(lipidsum),
              sd = sd(lipidsum))
  data_sub_classSum_p <- split(data_sub_classSum_stat, data_sub_classSum_stat$Class)
  data_sub_classSum_stat3 <- lapply(data_sub_classSum_p, getPValue, "headgroup", controlGrp)
  data_sub_classSum_stat3 <- do.call(rbind, data_sub_classSum_stat3)
  sigLabel <- addSigLabel(data_sub_classSum_stat3$p)
  data_sub_classSum_stat3 <- cbind(data_sub_classSum_stat3, sigLabel = sigLabel)
  data_sub_classSum_integStat <- left_join(data_sub_classSum_stat2, data_sub_classSum_stat3)
  nclass <- length(unique(data_sub_classSum_integStat$Class))
  h_plot <- nclass / 3 * (12/7)
  h2_plot <- nclass / 3 * (16/7)
  
  ## Visualization 
  data_sub_classSum_integStat$group <- 
    factor(data_sub_classSum_integStat$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  data_sub_classSum_stat$group <- 
    factor(data_sub_classSum_stat$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  plot_color <- ggplot() +
    geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean, fill = group), stat = "identity") +
    geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd, color = group), width = 0.2) +
    geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
    geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+1.5*sd, label = sigLabel),
              size = 3, fontface = "bold", color = "red") +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    facet_wrap(~Class, scales="free", ncol = 3) +
    scale_fill_aaas() +
    scale_color_aaas() +
    labs(x = "group",
         y = "total concentration", 
         color = "group", 
         title = "Lipid class statistics") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(fileLoc, "headgroup_color_", pname, ".pdf"), 
         plot = plot_color, device = "pdf", width = 9, height = h_plot)
  
  plot_black <- ggplot() +
    geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean), fill = "white", color = "black", stat = "identity", 
             width = 0.5, position=position_dodge(10)) +
    geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd), color = "black", width = 0.2) +
    geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center', 
                 binwidth = 5) +
    #!!!!!WARNING: y = realmean+1.5*sd may make some sigLabel be covered in this "black" plot
    geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+1.5*sd, label = sigLabel),
              size = 3, fontface = "bold", color = "red") +
    theme_bw() + 
    theme(axis.text.x=element_text(angle=45, hjust=1), 
          line = element_line(colour = "black", size = 1, 
                              linetype = 1, lineend = "butt")) +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~Class, scales="free", ncol = 3) +
    labs(x = "group",
         y = "total concentration", 
         title = "Lipid class statistics") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(paste0(fileLoc, "headgroup_black_", pname, ".pdf"), plot = plot_black, 
         device = "pdf", width = 6, height = h2_plot)
}

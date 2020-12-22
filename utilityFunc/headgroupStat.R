headgroupStat <- function(dataSet, mSet,  
                          fileLoc, ignore = T, stat = F){
  dataType <- dataSet$dataType
  allgroups <- dataSet$allgroups
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  
  ## Source will offer the following contents:
  ## Function(s): getClassInfo
  source("./utilityFunc/getClassInfo.R")
  ## Source will offer the following contents:
  ## Function(s): getPValue
  source("./utilityFunc/getPValue.R")
  ## Source will offer the following contents:
  ## Function(s): addSigLabel
  source("./utilityFunc/addSigLabel.R")
  
  data_tidy <- dataSet[["data"]] %>%
    mutate(lipidName = dataSet[["lipidName"]], 
           Class = switch(dataSet$dataType,
                          LipidSearch = sapply(lipidName, getClassInfo, "LipidSearch", ignore = ignore), 
                          MS_DIAL = sapply(lipidName, getClassInfo, "MS_DIAL", ignore = ignore)))
  
  ### Tidy for and do Visualization ###
  ## Use data_tidy(Delete duplication) to calculate itensity of each lipid class
  data_sub_classSum_stat <- data_tidy %>%
    ungroup() %>%
    select(-lipidName) %>%
    gather(key = "case", value = "lipidsum", -Class) %>%
    group_by(Class, case) %>%
    summarise(lipidsum = sum(lipidsum, na.rm = T)) %>%
    mutate(group = allgroups[match(case, names(allgroups))]) 
  if(stat == T){
    #Use for "statClassSum"
    data_classSum <- data_sub_classSum_stat %>%
      spread(Class, lipidsum) %>%
      select(c(-case, -group)) %>%
      t()
    return(data_classSum)
  }
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
  h_plot <- ifelse(nclass > 3, nclass, 3) / 3 * (30/9)
  h2_plot <- ifelse(nclass > 3, nclass, 3) / 3 * (16/7)
  w_plot <- length(groupsLevel) * (ifelse(nclass >= 3, 9, nclass*(9/3)+1)/3)
  w2_plot <- length(groupsLevel) * (6/3)
  w3_plot <- length(groupsLevel) * (ifelse(nclass >= 3, 6, nclass*(6/3)+1)/3)
  
  ## Visualization 
  data_sub_classSum_integStat$group <- 
    factor(data_sub_classSum_integStat$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  data_sub_classSum_stat$group <- 
    factor(data_sub_classSum_stat$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  # plot_color <- ggplot() +
  #   geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean, fill = group), stat = "identity") +
  #   geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd, color = group), width = 0.2) +
  #   geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
  #   geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+1.5*sd, label = sigLabel),
  #             size = 3, fontface = "bold", color = "red") +
  #   theme_bw() +
  #   theme(axis.text.x = element_blank(), 
  #         axis.ticks.x = element_blank()) +
  #   facet_wrap(~Class, scales="free", ncol = 3) +
  #   scale_fill_npg() +
  #   scale_color_npg() +
  #   labs(x = "group",
  #        y = "total concentration", 
  #        color = "group", 
  #        title = "Lipid class statistics") +
  #   theme(plot.title = element_text(hjust = 0.5, size = 20))
  headgroup_out <- data_sub_classSum_stat %>% 
    select(-group) %>% 
    spread(key = "case", value = "lipidsum")
  write.csv(headgroup_out, paste0(fileLoc, "lipidClass_conc.csv"), row.names = F)
  data_sub_classSum_stat_split <- split(data_sub_classSum_stat, f = data_sub_classSum_stat$Class)
  for(i in 1:length(data_sub_classSum_stat_split)){
    data_sub_classSum_stat_split[[i]]$group <- 
      factor(data_sub_classSum_stat_split[[i]]$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
    Classi <- data_sub_classSum_stat_split[[i]]$Class[1]
    p <- ggboxplot(data_sub_classSum_stat_split[[i]], x = "group", y = "lipidsum",
                   color = "group", add = "jitter")+ # Add global p-value
      stat_compare_means(aes(label = ..p.signif..),
                         method = "t.test", ref.group = controlGrp,
                         size = 5, fontface = "bold")+
      scale_color_npg() +
      labs(title = Classi, 
           x = "group",
           y = "total concentration") +
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5, size = 20), 
            axis.title = element_text(size = 15), 
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12)) 
    ggsave(paste0(fileLoc, Classi, ".pdf"), p, width = w2_plot, height = 12)
  }
  
  
  plot_cum <- ggplot(data = data_sub_classSum_integStat, aes(x = group, y = realmean, fill = Class)) + 
    geom_bar(stat="identity", width = 0.7)+
    scale_fill_npg() +
    theme_classic()+
    labs(title = "Cumulative lipid composition", 
         x = "group",
         y = "total concentration") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))
  ggsave(paste0(fileLoc, "headgroup_cum_", pname, ".pdf"), 
         plot = plot_cum, device = "pdf", 
         width = 4/9.2*(6.3/2*length(groupsLevel)+2.9), height = 10)
  
  plot_color <- ggboxplot(data_sub_classSum_stat, x = "group", y = "lipidsum",
            color = "group", add = "jitter")+ # Add global p-value
    stat_compare_means(aes(label = ..p.signif..),
                       method = "t.test", ref.group = controlGrp)+ 
    facet_wrap(~Class, scales="free", ncol = 3) +
    scale_color_npg()+
    labs(x = "group",
         y = "total concentration", 
         color = "group", 
         title = "Lipid class statistics") +
    theme(
      strip.background = element_blank(), 
      legend.position = "right", 
      plot.title = element_text(hjust = 0.5, size = 20), 
      axis.title = element_text(size = 15), 
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 12)
    )+
    scale_y_continuous(expand = c(0.2,0))
  ggsave(paste0(fileLoc, "headgroup_color_", pname, ".pdf"), 
         plot = plot_color, device = "pdf", width = w_plot, height = h_plot)
  
  plot_black <- ggplot() +
    geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean), fill = "white", color = "black", stat = "identity", 
             width = 0.5, position=position_dodge(10)) +
    geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd), color = "black", width = 0.2) +
    geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
    geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+1.5*sd, label = sigLabel),
              size = 3, fontface = "bold", color = "red") +
    theme_bw() + 
    theme(axis.text.x=element_text(angle=45, hjust=1), 
          line = element_line(colour = "black", size = 1, 
                              linetype = 1, lineend = "butt")) +
    facet_wrap(~Class, scales="free", ncol = 3) +
    labs(x = "group",
         y = "total concentration", 
         title = "Lipid class statistics") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.title = element_text(size = 15))
  ggsave(paste0(fileLoc, "headgroup_black_", pname, ".pdf"), plot = plot_black, 
         device = "pdf", width = w3_plot, height = h2_plot)
}

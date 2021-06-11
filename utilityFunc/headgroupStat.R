headgroupStat <- function(dataSet, mSet,  
                          fileLoc, ignore = T, stat = F){
  dataType <- dataSet$dataType
  allgroups <- dataSet$allgroups
  controlGrp <- dataSet$controlGrp
  groupsLevel <- dataSet$groupsLevel
  pname <- ifelse(length(groupsLevel) > 2, "all", paste0(groupsLevel[groupsLevel != controlGrp], 
                                                         "_vs_", controlGrp))
  
  source("./utilityFunc/plottingPalettes.R")
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
    dplyr::select(-lipidName) %>%
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
  # data_sub_classSum_p <- split(data_sub_classSum_stat, data_sub_classSum_stat$Class)
  # data_sub_classSum_stat3 <- lapply(data_sub_classSum_p, getPValue, "headgroup", controlGrp)
  # data_sub_classSum_stat3 <- do.call(rbind, data_sub_classSum_stat3)
  # sigLabel <- addSigLabel(data_sub_classSum_stat3$p)
  # data_sub_classSum_stat3 <- cbind(data_sub_classSum_stat3, sigLabel = sigLabel)
  # data_sub_classSum_integStat <- left_join(data_sub_classSum_stat2, data_sub_classSum_stat3)
  
  nclass <- length(unique(data_sub_classSum_stat2$Class))
  ##used for color plot
  #h_plot <- 65/15.3*(15.1/8*ifelse(nclass%%3, nclass%/%3+1, nclass%/%3)+0.2)
  h_plot <- 33/1264*(1224/4*ifelse(nclass%%3, nclass%/%3+1, nclass%/%3)+40)
  ##used for black plot
  #h2_plot <- ifelse(nclass > 3, nclass, 3) / 3 * (16/7)
  ##used for color plot
  #w_plot <- 21/17.2*(1.8+((15.4/3-0.9)/3*length(groupsLevel)+0.9)*ifelse(nclass >= 3, 3, nclass))
  w_plot <- 25/956*(82+(259/5*length(groupsLevel)+(288-259))*ifelse(nclass >= 3, 3, nclass))
  ##used for one class plot
  #w2_plot <- 8/15.9*(9.7/3*length(groupsLevel)+6.2)
  w2_plot <- 8/15.9*(9.7/3*length(groupsLevel)+6.2)
  ##used for black plot
  #w3_plot <- length(groupsLevel) * (ifelse(nclass >= 3, 6, nclass*(6/3)+1)/3)
  
  ## Visualization 
  data_sub_classSum_stat2$group <- 
    factor(data_sub_classSum_stat2$group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
  # data_sub_classSum_integStat$group <- 
  #   factor(data_sub_classSum_integStat$group, levels = c(controlGrp, unique(allgroups[allgroups != controlGrp])))
  data_sub_classSum_stat$group <- 
    factor(data_sub_classSum_stat$group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
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
  all_pvalues <- c()
  cplist <- mapply(c, controlGrp, groupsLevel[groupsLevel != controlGrp], 
                   SIMPLIFY = F, USE.NAMES = F)
  for(i in 1:length(data_sub_classSum_stat_split)){
    data_sub_classSum_stat_split[[i]]$group <- 
      factor(data_sub_classSum_stat_split[[i]]$group, levels = c(controlGrp, groupsLevel[groupsLevel != controlGrp]))
    Classi <- data_sub_classSum_stat_split[[i]]$Class[1]
    p <- ggplot(data = data_sub_classSum_stat_split[[i]], aes(x = group, y = lipidsum, color = group)) + 
      #use width = 0.5 before
      geom_boxplot(width = 0.8, size = 1.2)+ 
      #geom_jitter(color = "black", position=position_jitter(0.2)) +
      scale_color_npg() +
      labs(title = Classi, 
           x = "group",
           y = "total concentration") +
      theme_classic() +
      #use size = 25 before
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
            axis.title = element_text(size = 30), 
            axis.text = element_text(size = 30),
            legend.text = element_text(size = 30), 
            legend.title = element_text(size = 30)) 
    if(length(groupsLevel) == 2){
      # p <- ggboxplot(data_sub_classSum_stat_split[[i]], x = "group", y = "lipidsum",
      #                color = "group", add = "jitter")+ # Add global p-value
      #   stat_compare_means(aes(label = ..p.signif..),
      #                      method = "t.test", ref.group = controlGrp,
      #                      size = 5, fontface = "bold")+
      #   scale_color_npg() +
      #   labs(title = Classi, 
      #        x = "group",
      #        y = "total concentration") +
      #   theme(legend.position = "right", 
      #         plot.title = element_text(hjust = 0.5, size = 20), 
      #         axis.title = element_text(size = 15), 
      #         axis.text = element_text(size = 12),
      #         legend.text = element_text(size = 12), 
      #         legend.title = element_text(size = 12)) 
      p <- p + ggsignif::geom_signif(comparisons = cplist, color = "black", map_signif_level = T, 
                            step_increase = 0.05, tip_length = 0, textsize = 8, test = "t.test")
    }else if(length(groupsLevel) > 2){
      if(!controlGrp %in% unique(data_sub_classSum_stat_split[[i]]$group)){
        myp <- NULL
      }else if(length(unique(data_sub_classSum_stat_split[[i]]$group)) == 1){
        myp <- rep(NA, length(groupsLevel)-1)
      }else{
        ## Dunnett’s multiple comparison test in one-way ANOVA
        data_ano <- aov(lipidsum ~ group, data = data_sub_classSum_stat_split[[i]])
        #summary(data_ano)
        multcp <- glht(data_ano, linfct=mcp(group="Dunnett"),alternative="two.side") 
        multcp_sum <- summary(multcp)
        #gps <- levels(multcp_sum$model$model$group)
        myp <- multcp_sum$test$pvalues
        if(length(myp) != length(groupsLevel)-1){
          notins <- which(is.na(match(groupsLevel, unique(data_sub_classSum_stat_split[[i]]$group))))-1
          ins <- which(groupsLevel %in% unique(data_sub_classSum_stat_split[[i]]$group))[-1]-1
          myp2 <- numeric(length = length(groupsLevel)-1)
          myp2[ins] <- myp
          myp2[notins] <- NA
          myp <- myp2
        }
      }
      all_pvalues <- c(all_pvalues, myp)
      # stat.test <- tibble(
      #   group1 = controlGrp, 
      #   group2 = gps[gps != controlGrp], 
      #   p = multcp_sum$test$pvalues, 
      #   p.signif = addSigLabel(p),
      #   y.position = max(data_sub_classSum_stat_split[[i]]$lipidsum) + 
      #     sd(data_sub_classSum_stat_split[[i]]$lipidsum)
      # )
      myi <<- 0
      seq_gen <<-function(x, y) {
        myi <<- myi + 1L
        return(list(p.value = multcp_sum$test$pvalues[myi]))
      }
      p <- p + 
        ggsignif::geom_signif(comparisons = cplist, color = "black", map_signif_level = T, test = "seq_gen", 
                              step_increase = 0.05, tip_length = 0, textsize = 8)
    }
    ggsave(paste0(fileLoc, "others_", Classi, ".pdf"), p, width = w2_plot, height = 12)

  }
  
  nClass <- length(unique(data_sub_classSum_stat$Class))
  plot_cum <- ggplot(data = data_sub_classSum_stat2, aes(x = group, y = realmean, fill = Class)) + 
    geom_bar(stat="identity", width = 0.7)+
    theme_classic()+
    labs(title = "Cumulative lipid composition", 
         x = "group",
         y = "total concentration") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25), 
          axis.title = element_text(size = 25), 
          axis.text = element_text(size = 25),
          legend.text = element_text(size = 25), 
          legend.title = element_text(size = 25), 
          axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +
    scale_fill_manual(values = plottingPalettes(nClass, type = "discrete"))
  ggsave(paste0(fileLoc, "headgroup_cum_", pname, ".pdf"), 
         plot = plot_cum, device = "pdf", 
         #width = 6/7.7*(3.3/3*length(groupsLevel)+4.4), 
         width = 6/7.7*(3.3/3*length(groupsLevel)+8.8), 
         height = 10)
  
  # plot_color <- ggboxplot(data_sub_classSum_stat, x = "group", y = "lipidsum",
  #           color = "group", add = "jitter")+ # Add global p-value
  #   stat_compare_means(aes(label = ..p.signif..),
  #                      method = "t.test", ref.group = controlGrp)+ 
  #   facet_wrap(~Class, scales="free", ncol = 3) +
  #   scale_color_npg()+
  #   labs(x = "group",
  #        y = "total concentration", 
  #        color = "group", 
  #        title = "Lipid class statistics") +
  #   theme(
  #     strip.background = element_blank(), 
  #     legend.position = "right", 
  #     plot.title = element_text(hjust = 0.5, size = 20), 
  #     axis.title = element_text(size = 15), 
  #     axis.text = element_text(size = 12),
  #     legend.text = element_text(size = 12), 
  #     legend.title = element_text(size = 12)
  #   )+
  #   scale_y_continuous(expand = c(0.2,0))
  data_sub_classSum_stat$Class <- factor(data_sub_classSum_stat$Class, 
                                         levels = names(data_sub_classSum_stat_split)) 
  plot_color <- ggplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum, color = group)) + 
    geom_boxplot(width = 0.8, size = 1.2)+ 
    facet_wrap(~Class, scales="free", ncol = 3) +
    scale_color_npg()+
    labs(x = "group",
         y = "total concentration",
         color = "group",
         title = "Lipid class statistics") +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"), 
      axis.title = element_text(size = 30), 
      axis.text = element_text(size = 30),
      legend.text = element_text(size = 30), 
      legend.title = element_text(size = 30), 
      strip.text = element_text(hjust = 0.5, size = 30, face = "bold")
    ) 
  if(length(groupsLevel) == 2){
    plot_color <- plot_color + 
      ggsignif::geom_signif(comparisons = cplist, color = "black", map_signif_level = T,  
                            step_increase = 0.1, tip_length = 0, textsize = 8, test = "t.test")
  }else if(length(groupsLevel) > 2){
    ## Dunnett’s multiple comparison test in one-way ANOVA(complete above)
    myi <<- 0
    seq.gen <<-function(x, y) {
      myi <<- myi + 1L
      return(list(p.value = all_pvalues[myi]))
    }
    plot_color <- plot_color + 
      ggsignif::geom_signif(comparisons = cplist, color = "black", map_signif_level = T, test = "seq.gen", 
                            step_increase = 0.1, tip_length = 0, textsize = 8)
  }
  ggsave(paste0(fileLoc, "headgroup_color_", pname, ".pdf"), 
         plot = plot_color, device = "pdf", width = w_plot, height = h_plot, limitsize = FALSE)
  
  # plot_black <- ggplot() +
  #   geom_bar(data = data_sub_classSum_integStat, aes(x = group, y = realmean), fill = "white", color = "black", stat = "identity", 
  #            width = 0.5, position=position_dodge(10)) +
  #   geom_errorbar(data = data_sub_classSum_integStat, aes(x = group, ymin = realmean, ymax = realmean + sd), color = "black", width = 0.2) +
  #   geom_dotplot(data = data_sub_classSum_stat, aes(x = group, y = lipidsum), binaxis='y', stackdir='center') +
  #   geom_text(data = data_sub_classSum_integStat, aes(x = group, y = realmean+1.5*sd, label = sigLabel),
  #             size = 3, fontface = "bold", color = "red") +
  #   theme_bw() + 
  #   theme(axis.text.x=element_text(angle=45, hjust=1), 
  #         line = element_line(colour = "black", size = 1, 
  #                             linetype = 1, lineend = "butt")) +
  #   facet_wrap(~Class, scales="free", ncol = 3) +
  #   labs(x = "group",
  #        y = "total concentration", 
  #        title = "Lipid class statistics") +
  #   theme(plot.title = element_text(hjust = 0.5, size = 20),
  #         axis.title = element_text(size = 15))
  # ggsave(paste0(fileLoc, "headgroup_black_", pname, ".pdf"), plot = plot_black, 
  #        device = "pdf", width = w3_plot, height = h2_plot)
}

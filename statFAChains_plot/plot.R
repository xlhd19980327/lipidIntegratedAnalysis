library(plotly)
library(ggrepel)
#1 gene-P_R-fill=-log10PValue
regStat_gene_tidy <- regStat_gene %>%
  filter(!is.na(p_value)) %>%
  mutate(log10p = -log10(p_value)) #%>%
ggplot(data = regStat_gene_tidy, aes(x = P_R, y = reorder(gene, log10p), fill = log10p)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_gradient2(low = "yellow",
                       #mid = "white",
                       high = "red",
                       midpoint = 2-log10(5))+
  labs(x = "Change Ratio", y = "gene", fill = "-log10PValue")

#2 pathway-P_R-fill=-log10PValue
regStat_path_tidy <- regStat_path %>%
  filter(!is.na(p_value)) %>%
  mutate(log10p = -log10(p_value)) #%>%
ggplot(data = regStat_path_tidy, aes(x = P_R, y = reorder(pathway, log10p), fill = log10p)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_gradient2(low = "yellow",
                       #mid = "white",
                       high = "red",
                       midpoint = 2-log10(5))+
  labs(x = "Enrichment Ratio", y = "pathway", fill = "-log10PValue")

#3 p-FC-text=gene-size=pathway_count-(anno=pathways)
regStat_gene_info <- regStat_gene$gene_info %>%
  select(gene, pathway) %>%
  unique() %>% 
  group_by(gene) %>%
  summarise(
    count = n(),
    pathways = paste0(pathway, collapse = ";")
  )
regStat_gene_integ <- left_join(regStat_gene$regState, regStat_gene_info)
regStat_gene_integ_tidy <- regStat_gene_integ %>%
  filter(!is.na(p_value)) %>%
  mutate(log10p = -log10(p_value))
ggplot(data = regStat_gene_integ_tidy, aes(x = P_R, y = log10p)) +
  geom_point(mapping = aes(size = count), position = position_jitter()) +
  theme_classic() +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5) + 
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line.y = element_blank(),
    panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20")
  ) +
  geom_text_repel(mapping = aes(label = gene), 
                  hjust = 0, nudge_x = 0.05) +
  labs(x = "Change Ratio", y = "-log10PValue", size = "count") 
ggplotly(p3)

#4 p-pathway_count-text=gene-(anno=pathways)
ggplot(data = regStat_gene_integ_tidy, aes(x = ifelse(P_R > 0, count, -count), y = log10p)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5) + 
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line.y = element_blank(),
    panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20")
  ) +
  geom_text_repel(mapping = aes(label = gene), 
                  hjust = 0, nudge_x = 0.05) +
  labs(x = "Count", y = "-log10PValue") 


#5 p-FC-text=pathway-size=gene_count-(anno=genes)
regStat_pathway_info <- regStat_path$path_info %>%
  select(gene, pathway) %>%
  unique() %>% 
  group_by(pathway) %>%
  summarise(
    count = n(),
    genes = paste0(gene, collapse = ";")
  )
regStat_path_integ <- left_join(regStat_path$regState, regStat_pathway_info)
regStat_path_integ_tidy <- regStat_path_integ %>%
  filter(!is.na(p_value)) %>%
  mutate(log10p = -log10(p_value))
ggplot(data = regStat_path_integ_tidy, aes(x = P_R, y = log10p)) +
  geom_point(mapping = aes(size = count), position = position_jitter()) +
  theme_classic() +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5) + 
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line.y = element_blank(),
    panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20")
  ) +
  geom_text_repel(mapping = aes(label = pathway), 
                  hjust = 0, nudge_x = 0.05) +
  labs(x = "Enrichment Ratio", y = "-log10PValue", size = "count") 

#6 p-gene_count-text=pathway-(anno=genes)
ggplot(data = regStat_path_integ_tidy, aes(x = ifelse(P_R > 0, count, -count), y = log10p)) +
  geom_point() +
  theme_classic() +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.5) + 
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.line.y = element_blank(),
    panel.grid = element_line(colour = "grey92"), 
    panel.grid.minor = element_line(size = rel(0.5)), 
    strip.background = element_rect(fill = "grey85", 
                                    colour = "grey20")
  ) +
  geom_text_repel(mapping = aes(label = pathway), 
                  hjust = 0, nudge_x = 0.05) +
  labs(x = "Count", y = "-log10PValue") 

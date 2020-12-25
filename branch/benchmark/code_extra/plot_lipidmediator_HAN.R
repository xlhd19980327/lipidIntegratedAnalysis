library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
data <- read.csv("./branch/benchmark/input/HANLipidMediator_forPlotting.csv", 
                 na.strings = "")
data_tidy <- data %>%
  group_by(Compound, Class) %>%
  gather(key = "sample", value = "conc", -c("Compound", "Class")) %>%
  mutate(group = gsub("_[0-9]$", "", sample))
plot_cum <- ggplot(data = data_tidy, aes(x = group, y = conc, 
                                         fill = factor(Class, 
                                                       levels = c("Leukotrienes", "Prostaglandins", "Lipoxins", "D-series resolvins", "E-series resolvins")))) + 
  geom_bar(stat="identity", width = 0.7)+
  scale_fill_npg() +
  theme_classic()+
  labs(title = "Cumulative composition", 
       x = "group",
       y = "total concentration") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_blank())
ggsave("~/temp/lipidmediator_cum.pdf", 
       plot = plot_cum, device = "pdf", 
       width = 6, height = 10)

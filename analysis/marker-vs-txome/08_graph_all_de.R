library(tidyverse)
library(cowplot)

theme_set(theme_cowplot(font_size = 11))


sce_list <- readRDS("data/scelist/scesets.rds")
datasets <- names(sce_list)
rm(sce_list)

tpr_df <- read_csv("data/marker-vs-txome/tidy_de.csv")

n_additional_markers <- c(1,5,10,20,50,100,500,1000)
tpr_df <- mutate(tpr_df, n_additional_markers = factor(n_additional_markers[tpr_df$nam_index]))%>% 
  mutate(dataset = datasets[dataset])

tpr_df$algorithm <- stringr::str_to_title(as.character(tpr_df$algorithm))
tpr_df$algorithm <- plyr::mapvalues(tpr_df$algorithm, from = c("Pc1", "Tscan", "Dpt"),
                                    to = c("PC1", "TSCAN", "DPT"))
tpr_df$algorithm <- factor(as.character(tpr_df$algorithm), levels = c("Ouija", "PC1",
                                                                      "Monocle", "TSCAN", "DPT"))

tpr_df$dataset <- plyr::mapvalues(tpr_df$dataset,
                                  to = c("Trapnell", "Shin", "Zhou"),
                                  from = c("trapnell", "shin", "hsc"))



txwide_plot <- filter(tpr_df, !(grepl("marker", pst_str)), algorithm != "ouija") %>% 
  ggplot(aes(x = n_additional_markers, #group = n_additional_markers, 
             y = tpr, fill = algorithm, color = algorithm)) +
  geom_boxplot() + facet_wrap(~ dataset, ncol = 1) +
  scale_fill_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  ylab("Proportion overlap w/ full txome DE") +
  xlab("Number of additional genes")

marker_plot <- filter(tpr_df, grepl("marker", pst_str) | algorithm == "ouija") %>% 
  ggplot(aes(x = n_additional_markers, #group = n_additional_markers, 
             y = tpr, fill = algorithm, color = algorithm)) +
  geom_boxplot() + facet_wrap(~ dataset, ncol = 1) +
  scale_fill_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  ylab("Proportion overlap w/ marker DE") +
  xlab("Number of additional genes")

cowplot::plot_grid(txwide_plot, marker_plot, ncol = 1, labels = "AUTO")

ggsave("figs/marker-vs-txome-de-results.png", width=10,height=10)

## Line plots

# txwide_plot_line <- filter(tpr_df, !(grepl("marker", pst_str)), algorithm != "ouija") %>% 
#   group_by(dataset, algorithm, n_additional_markers) %>% 
#   summarise(mean_tpr = mean(tpr, na.rm = TRUE)) %>% 
#   ggplot(aes(x = n_additional_markers, 
#              y = mean_tpr, group = algorithm, color = algorithm)) +
#   geom_line() + facet_wrap(~ dataset, nrow = 1) +
#   scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
#   ylab("Proportion overlap w/ full txome DE") +
#   xlab("Number of additional genes")
# 
# marker_plot_line_de <- filter(tpr_df, grepl("marker", pst_str) | algorithm == "ouija") %>% 
#   group_by(dataset, algorithm, n_additional_markers) %>% 
#   summarise(mean_tpr = mean(tpr, na.rm = TRUE)) %>%  
#   ggplot(aes(x = n_additional_markers, #group = n_additional_markers, 
#              y = mean_tpr, group = algorithm, color = algorithm)) +
#   geom_line() + facet_wrap(~ dataset, nrow = 1) +
#   scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
#   ylab("Proportion overlap w/ marker DE") +
#   xlab("Number of additional genes")
# 
# save(txwide_plot_line, marker_plot_line_de, file = "de_plots.Rdata")


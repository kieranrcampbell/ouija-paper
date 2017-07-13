library(ggplot2)

heatmap <- readRDS("figs/trapnell_correlation.rds")
examples <- readRDS("figs/trapnell_example_genes.rds")
correlations <- readRDS("figs/marker_correlations.rds")


heatmap <- heatmap + labs(x = "Gene", y = "Gene") +
  theme(axis.title = element_text(size = 10),
        axis.title.y = element_text(margin = margin(0, -8, 0, 5)),
        axis.title.x = element_text(margin = margin(-8, 0, 15, 0)),
        legend.box.margin = margin(l = -20))
# last_plot()

examples <- examples + theme(axis.title = element_text(size = 10))
correlations <- correlations + theme(axis.title = element_text(size = 10))

left_grid <- cowplot::plot_grid(heatmap, correlations, ncol = 1, rel_heights = c(2,1),
                                labels = c("A", "B"), label_size = 12)

fig1 <- cowplot::plot_grid(left_grid, examples, nrow = 1)

ggsave("figs/fig1.png", width = 12, height = 8)

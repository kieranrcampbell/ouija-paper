
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
library(readr)

source("analysis/transient/common_simulation_functions.R")

# theme_set(theme_classic())


trans_dir <- "data/transient/results/"

trans_files <- dir(trans_dir, full.names = TRUE)

trans_df_list <- lapply(trans_files, read_csv)

trans_df <- bind_rows(trans_df_list)

trans_df$G <- factor(trans_df$G)
trans_df$prop_switch <- factor(trans_df$prop_switch)

ggplot(trans_df, aes(x = G, y = abscor, fill = prop_switch)) +
  geom_boxplot(width = 0.4) +
  ylab(expression("Pearson" ~ rho ~ "to true pseudotime")) +
  xlab("Number of genes") +
  scale_fill_brewer(palette = "Set2", name = "Proportion of switch-like genes") +
  theme(legend.position = "bottom",
        legend.background = element_rect(linetype = 'solid', color = 'grey30'),
        axis.title = element_text(size = 10))

trans_boxplot <- last_plot()


# Now make an example figure

mu0 <- 2
t0 <- 0.6
k <- 10
p <- 0.4
b <- 5


set.seed(1965)
pst <- runif(100, 0, 1)
switch_exprs <- rsinglecell_sigmoid(k, mu0, t0, pst)
trans_exprs <- rsinglecell_transient(mu0, b, p, pst)

switch_mean <- sigmoid(pst, k, mu0, t0)
trans_mean <- transient(pst, mu0, b, p)

switch_df <- data_frame(pst, expression = switch_exprs, mean = switch_mean, type = "Switch")
transient_df <- data_frame(pst, expression = trans_exprs, mean = trans_mean, type = "Transient")
df <- bind_rows(switch_df, transient_df)
df <- arrange(df, pst)

ggplot(df, aes(x = pst, color = type)) +
  geom_point(aes(y = expression), alpha = 0.7) +
  geom_line(aes(y = mean), size = 1.2) +
  scale_colour_brewer(palette = "Set1", name = "Behaviour") +
  labs(x = "Pseudotime", y = "Expression") +
  theme(legend.position = "bottom",
        legend.background = element_rect(linetype = 'solid', color = 'grey30'),
        axis.title = element_text(size = 10))

trans_exampleplot <- last_plot()


cor_plot <- readRDS("figs/li_cor.rds")
correct_genes <- readRDS("figs/li_correct_genes.rds")
incorrect_genes <- readRDS("figs/li_incorrect_genes.rds")

cor_plot <- cor_plot + theme(axis.title = element_text(size = 10))

# left_grid <- cowplot::plot_grid(correct_genes, incorrect_genes,
#                                 nrow = 2, labels = c("A", "B"),
#                                 label_size = 11, rel_heights = c(2,1))
# middle_grid <- cowplot::plot_grid(cor_plot, trans_exampleplot,
#                                   nrow = 2, labels = c("C", "D"),
#                                   label_size = 11, rel_heights = c(2,3))
# 
# right_grid <- cowplot::plot_grid(NULL, trans_boxplot, NULL,
#                                  rel_heights = c(1,4,1), ncol = 1)
# 
# cowplot::plot_grid(left_grid, middle_grid, right_grid, 
#                    nrow = 1, labels = c("", "", "E"),
#                    label_size = 11, rel_widths = c(2,3,3))
# 
# ggsave("~/Desktop/transient.png", width = 12, height = 6)
# 
# middle_grid2 <- cowplot::plot_grid(incorrect_genes, cor_plot,
#                                 nrow = 2, labels = c("B", "C"),
#                                 label_size = 11, rel_heights = c(1,1))
# 
# cowplot::plot_grid(correct_genes, middle_grid2, trans_exampleplot,
#                    trans_boxplot, nrow = 1, labels = c("A", "", "D", "E"),
#                    label_size = 11)
# 
# ggsave("~/Desktop/transient.png", width = 14, height = 4)


incorrect_genes <- incorrect_genes + scale_y_continuous(breaks = c(0, 2, 4, 6))
left_grid <- cowplot::plot_grid(correct_genes, incorrect_genes,
                                cor_plot, labels = "AUTO", label_size = 11,
                                ncol = 1, rel_heights = c(2, 1.3, 2))
right_grid <- cowplot::plot_grid(trans_exampleplot, trans_boxplot,
                                 labels = c("D", "E"), label_size = 10,
                                 ncol = 1)

cowplot::plot_grid(left_grid, right_grid, ncol = 2,
                    rel_widths = c(2.3,3))

ggsave("figs/fig_transient.png", width = 9, height = 7)


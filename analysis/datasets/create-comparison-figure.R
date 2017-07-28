library(tidyverse)
library(stringr)

data_dir <- "data/mvt_csv/"

data_files <- dir(data_dir)
data_sets <- sapply(strsplit(data_files, ".", fixed = TRUE), `[`, 1)

data_path <- file.path(data_dir, data_files)

read_file <- function(dpath, dset) {
  dframe <- read_csv(dpath, col_types = "dddd")
  dframe <- mutate(dframe, dataset = dset)
  return(dframe)
}

dframes <- Map(read_file, data_path, data_sets)
df <- bind_rows(dframes)

df$dataset <- paste(str_to_title(df$dataset), "et al.")

df_tidy <- gather(df, algorithm, correlation, -dataset, -max_cor)
df_tidy$algorithm <- plyr::mapvalues(df_tidy$algorithm,
                                     from = c("monocle", "dpt", "tscan"),
                                     to = c("Monocle 2", "DPT", "TSCAN"))

filter(df_tidy, dataset != "Li et al.") %>% 
  ggplot(aes(x = dataset, y = abs(correlation), fill = algorithm)) +
    geom_bar(stat = 'identity', position = 'dodge', color = 'grey30', width = 0.7) +
    ylim(c(0,1)) +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.title = element_blank()) +
    ylab("Correlation to Ouija\nmarker pseudotime") +
    xlab("Dataset") +
    labs(subtitle = "Comparison to transcriptome-wide pseudotimes")

marker_correlations <- last_plot()

filter(df_tidy, dataset != "Li et al.") %>% 
  ggplot(aes(x = dataset, y = abs(correlation) / max_cor, fill = algorithm)) +
    geom_bar(stat = 'identity', position = 'dodge', color = 'grey30', width = 0.7) +
    ylim(c(0,2.05)) +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.title = element_blank()) +
    ylab("Correlation to Ouija\nmarker pseudotime") +
    xlab("Dataset") +
    labs(subtitle = "Comparison to transcriptome-wide pseudotimes")

relative_marker_correlations <- last_plot()

saveRDS(marker_correlations, "figs/marker_correlations.rds")
saveRDS(relative_marker_correlations, "figs/relative_marker_correlations.rds")

filter(df, dataset != "Li et al.") %>% 
  write_csv("data/max_correlations.csv")


# Comparison figure -------------------------------------------------------

data_dir <- "data/cor_comp/"

data_files <- dir(data_dir)
data_sets <- sapply(strsplit(data_files, ".", fixed = TRUE), `[`, 1)

data_path <- file.path(data_dir, data_files)

read_file <- function(dpath, dset) {
  dframe <- read_csv(dpath)
  dframe <- mutate(dframe, dataset = dset)
  return(dframe)
}

dframes <- Map(read_file, data_path, data_sets)
df <- bind_rows(dframes)

df$dataset <- paste(str_to_title(df$dataset), "et al.")

df <- mutate(df, has_ouija = grepl("ouija", row_alg) | grepl("ouija", col_alg)) %>% 
  dplyr::filter(abs(correlation) < 1)

ggplot(df, aes(x = dataset, y = abs(correlation), color = has_ouija)) +
  #stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #            geom = "crossbar", width = 0.5, alpha = 0.5) +
  geom_point() +
  scale_color_manual(values = c("black", "red"),
                     name = "Comparison with\nOuija markers") +
  labs(x = "Dataset", y = "Absolute correlation") 

ggsave("figs/supp_correlation_comparison.png", width = 6, height = 4)

fit <- lm(abs(correlation) ~ has_ouija + dataset, data = df)
s <- summary(fit)
as_data_frame(s$coefficients) %>% 
  mutate(quantity = rownames(s$coefficients)) %>% 
  write_csv("data/correlation_comparison_linear_model.csv")

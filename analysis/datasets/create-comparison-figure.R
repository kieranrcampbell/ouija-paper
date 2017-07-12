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
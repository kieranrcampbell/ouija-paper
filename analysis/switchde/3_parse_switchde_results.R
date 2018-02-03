
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(aargh))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(cowplot))

plot_switchde_results <- function(input_directory = "../../data/switchde/sde",
                                  output_png = "../../figs/switchde/fig_switchde.png") {
                                  
  
  trace_csvs <- dir(input_directory, full.names = TRUE)
  
  df <- map_df(trace_csvs, read_csv)
  
  signif_genes <- group_by(df, gene, gsn) %>% 
    summarise(all_signif = all(qval < 0.05)) %>% 
    filter(all_signif) %>% 
    .$gene
  
  dfs <- filter(df, gene %in% signif_genes)
  
  # Now extract summary stats
  
  dfg <- group_by(dfs, gene, gsn) %>% 
    summarise(k_median = median(k),
              k_lower = quantile(k, probs = 0.05),
              k_upper = quantile(k, probs = 0.95),
              t0_median = median(t0),
              t0_lower = quantile(t0, probs = 0.05),
              t0_upper = quantile(t0, probs = 0.95)) %>% 
    ungroup()
  
  dfgr <- filter(dfg,
                 abs(k_median) < 10,
                 abs(t0_median) < 10,
                 abs(k_lower) < 20,
                 abs(k_upper) < 20,
                 abs(t0_lower) < 20,
                 abs(t0_upper) < 20)
  
  set.seed(123L)
  df_sample_k <- df_sample_t0 <- sample_n(dfgr, 100)
  
  df_sample_k$gsn <- factor(df_sample_k$gsn,
                           levels = df_sample_k$gsn[order(df_sample_k$k_median)])
  
  pltk <- ggplot(df_sample_k, aes(x = gsn, y = k_median)) +
    geom_errorbar(aes(ymin = k_lower, ymax = k_upper)) +
    geom_point() +
    coord_flip() +
    labs(x = "Gene", y = "Switch strength") +
    theme(axis.text.y = element_text(size = 5))
  
  df_sample_t0$gsn <- factor(df_sample_t0$gsn,
                            levels = df_sample_t0$gsn[order(df_sample_t0$t0_median)])
  
  pltt0 <- ggplot(df_sample_t0, aes(x = gsn, y = t0_median)) +
    geom_errorbar(aes(ymin = t0_lower, ymax = t0_upper)) +
    geom_point() +
    coord_flip() +
    labs(x = "Gene", y = "Switch time") +
    theme(axis.text.y = element_text(size = 5))
  
  
  pltg <- plot_grid(pltk, pltt0, ncol = 2, labels = "AUTO")
  
  ggsave(output_png, width = 6, height = 9)
}

aargh(plot_switchde_results)





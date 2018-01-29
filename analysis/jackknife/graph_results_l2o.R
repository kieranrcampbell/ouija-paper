library(scater)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(cowplot)
library(tidyr)
theme_set(theme_cowplot(font_size = 10))


results_dir <- "../../data/jackknife/l2o"

results_csv <- dir(results_dir, full.names = TRUE)

xx <- t(sapply(results_csv, function(r) as.numeric(strsplit(r, "_")[[1]][c(4,6)])))
to_keep <- xx[,1] < xx[,2] # Only keep pairs appearing once

results_csv <- results_csv[to_keep]

read_csv_and_annotate_gene <- function(r) {
  gene_lo <- as.numeric(strsplit(r, "_")[[1]][c(4,6)])
  read_csv(r) %>% 
    mutate(gene1 = gene_lo[1],
           gene2 = gene_lo[2])
}

df <- map_df(results_csv, read_csv_and_annotate_gene)

sce_list <- readRDS("../../data/scesets_with_marker_pseudotime.rds")
sce <- sce_list[['trapnell']]
sce_marker <- sce[fData(sce)$is_marker, ]

oui_pst <- pData(sce_marker)[["ouija_pseudotime"]]

df_o <- mutate(df, ouija_pseudotime = rep(oui_pst, length(results_csv)),
               gene12 = paste0(gene1, "_", gene2))

group_by(df_o, gene12) %>% 
  summarise(cor = cor(pseudotime, ouija_pseudotime))

unique_combs <- unique(df_o$gene12)

# uc <- unique_combs[1]

gsn <- fData(sce_marker)$gene_short_name

plots <- lapply(unique_combs, function(uc) {
  dfuc <- filter(df_o, gene12 == uc)
  
  # flip pseudotime if cor < 0
  if(cor(dfuc$pseudotime, dfuc$ouija_pseudotime) < 0) {
    dfuc$pseudotime <- 1 - dfuc$pseudotime
  }
  
  exprs_df <- as_data_frame(t(exprs(sce_marker)))
  names(exprs_df) <- gsn
  exprs_df <- mutate(exprs_df, cell = sampleNames(sce_marker)) %>% 
    gather(gene, expression, -cell)
  
  dfg <- mutate(dfuc, gene1_str = gsn[gene1], gene2_str = gsn[gene2]) %>% 
    select(-gene1, -gene2, -gene12) %>% 
    gather(gene_ind, gene, -cell, -pseudotime,-ouija_pseudotime) %>% 
    rename(`Leave-two-out pseudotime` = pseudotime, `All marker pseudotime` = ouija_pseudotime) %>% 
    gather(pseudotime_type, pseudotime, -cell, -gene_ind, -gene)
  
  dfg2 <- inner_join(dfg, exprs_df)
  
  ggplot(dfg2, aes(x = pseudotime, y = expression)) +
    facet_grid(gene ~ pseudotime_type, scales = "free") +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE, color = 'red') +
    labs(x = "Pseudotime", y = "Expression") +
    theme(strip.background = element_rect(fill = 'grey90')) +
    theme(axis.text = element_text(size = 9)) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8))
})

plt <- plot_grid(plotlist = plots, nrow = 5, labels = "AUTO")

ggsave("../../figs/jackknife/l2o_expression.png", plt, width = 9, height = 14)
  

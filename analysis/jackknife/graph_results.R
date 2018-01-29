library(scater)
library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(cowplot)
library(tidyr)

theme_set(theme_cowplot(font_size = 11))

results_dir <- "../../data/jackknife/l1o"

results_csv <- dir(results_dir, full.names = TRUE)

read_csv_and_annotate_gene <- function(r) {
  gene_lo <- strsplit(r, "_")[[1]][4]
  read_csv(r) %>% 
    mutate(gene_lo = gene_lo)
}

df <- map_df(results_csv, read_csv_and_annotate_gene)

sce_list <- readRDS("../../data/scesets_with_marker_pseudotime.rds")
sce <- sce_list[['trapnell']]
sce_marker <- sce[fData(sce)$is_marker, ]

oui_pst <- pData(sce_marker)[["ouija_pseudotime"]]
gsn <- fData(sce_marker)$gene_short_name
gsn_df <- gsn[as.numeric(df$gene_lo)]
df_o <- mutate(df, ouija_pseudotime = rep(oui_pst, 5)) %>% 
  mutate(gene_lo_str = paste(gsn_df, "left out"),
         gene = gsn_df)

ggplot(df_o, aes(x = ouija_pseudotime, y = pseudotime)) +
  geom_point() +
  facet_wrap(~ gene_lo_str, nrow = 1) +
  labs(x = "All marker pseudotime",
       y = "Leave-one-out\npseudotime")

ggsave("../../figs/jackknife/l1o.png", width = 8, height = 3)


# Expression plots --------------------------------------------------------

dfe <- as_data_frame(t(exprs(sce_marker)))
names(dfe) <- gsn
dfe <- mutate(dfe, cell = sampleNames(sce_marker)) %>% 
  gather(gene, expression, -cell)

df_all <- inner_join(df_o, dfe, by = c("cell", "gene")) %>% 
  rename(`Leave-one-out pseudotime` = pseudotime, `All marker pseudotime` = ouija_pseudotime) %>% 
  select(-gene_lo_str, -gene_lo) %>% 
  gather(pseudotime_type, pseudotime, -cell, -gene, -expression)


ggplot(df_all, aes(x = pseudotime, y = expression)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE, color = 'red') +
  facet_grid(gene ~ pseudotime_type, scales = "free") +
  labs(x = "Pseudotime", y = "Expression") +
  theme(strip.background = element_rect(fill = 'grey90'))

ggsave("../../figs/jackknife/l1o_expression.png", width = 5, height = 7)













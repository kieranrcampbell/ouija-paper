
#' Invoke with
#' Rscript fit_differential_expression.R [dataset] [rep] [nam_index]

library(scater)
library(switchde)
library(tidyverse)

load("../../data/r1_scesets_with_marker_pseudotime.Rdata")

sce_list_with_marker_pseudotime <- readRDS("data/scesets/scesets_with_marker_pseudotime.rds")


all_dataset_de <- lapply(1:3, function(dset) {
  sce <- sce_list_with_marker_pseudotime[[dset]]
  
  pseudotime_str <- c("monocle_pseudotime", "tscan_pseudotime",
                      "ouija_pseudotime", "pc1_pseudotime",
                      "dpt_pseudotime", "pc1_marker_pseudotime",
                      "tscan_marker_pseudotime", "monocle_marker_pseudotime",
                      "dpt_marker_pseudotime")
  
  de_fits <- lapply(pseudotime_str, function(pst_str) {
    pst <- sce[[pst_str]]
    sde <- switchde(sce, pst)
    sde <- select(sde, gene, pval) %>% 
      mutate(pst_str = pst_str)
    return(sde)
  })
  
  de_fits_tidy <- bind_rows(de_fits)
  de_fits_tidy <- mutate(de_fits_tidy,
                         dataset = names(sce_list_with_marker_pseudotime)[dset])
  return(de_fits_tidy)
})

all_dataset_de_df <- bind_rows(all_dataset_de)

write_csv(all_dataset_de_df, "data/marker-vs-txome/txwide_marker_de.csv")






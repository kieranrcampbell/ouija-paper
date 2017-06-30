library(scater)
# library(tidyverse)

source("scripts/pseudotime_fitting_functions.R")

scelist <- readRDS("data/scesets/scesets.rds")


fit_all_pseudotime <- function(sce) {
  sce$monocle_pseudotime <- fit_monocle_pseudotime(sce)
  
  sce$tscan_pseudotime <- tryCatch(fit_tscan_pseudotime(sce),
                                   error = function(e) {
                                     message("TSCAN failed")
                                     return(rep(NA, ncol(sce)))
                                   })
  
  
  sce$ouija_pseudotime <- fit_ouija_pseudotime(sce)
  sce$pc1_pseudotime <- fit_pc1_pseudotime(sce)
  sce$dpt_pseudotime <- fit_dpt_pseudotime(sce)
  return( sce )
}



sce_list_with_pseudotime <- lapply(sce_list, fit_all_pseudotime)

saveRDS(sce_list_with_pseudotime, file = "data/scesets/scesets_with_pseudotime.rds")





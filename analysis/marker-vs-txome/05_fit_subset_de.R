
#' Invoke with
#' Rscript fit_differential_expression.R [dataset] [rep] [nam_index]

library(scater)
library(switchde)
library(tidyverse)

sce_list_with_marker_pseudotime <- readRDS("data/scesets/scesets_with_pseudotime.rds")
cmats <- readRDS("data/marker-vs-txome/pseudotimes_across_marker_subset_array.rds")

args <- commandArgs(trailingOnly = TRUE)

dataset <- as.numeric(args[1]) # which dataset
rep <- as.numeric(args[2]) # which rep
nam_index <- as.numeric(args[3]) # number of additional markers


sce <- sce_list_with_marker_pseudotime[[dataset]]
rm(sce_list_with_marker_pseudotime)

ouija_filename <- paste("ouija_pseudotime", dataset, rep, nam_index, sep = "_")
ouija_file_path <- file.path("data", "marker-vs-txome", "ouija_fits", paste0(ouija_filename, ".csv"))

ouija_hmc <- read_csv(ouija_file_path)
ouija_pseudotime <- ouija_hmc$pseudotime


monocle_pseudotime <- cmats[[dataset]][rep, nam_index, 1, ]
tscan_pseudotime <- cmats[[dataset]][rep, nam_index, 2, ]
pc1_pseudotime <- cmats[[dataset]][rep, nam_index, 4, ]
dpt_pseudotime <- cmats[[dataset]][rep, nam_index, 5, ]

if(any(is.na(ouija_pseudotime))) {
  sde_ouija <- NULL
} else {
  sde_ouija <- switchde(sce, ouija_pseudotime) %>% 
    mutate(algorithm = "ouija")
}

if(any(is.na(tscan_pseudotime))) {
  sde_tscan <- NULL
} else {
  sde_tscan <- switchde(sce, tscan_pseudotime) %>% 
    mutate(algorithm = "tscan")
}

if(any(is.na(monocle_pseudotime))) {
  sde_monocle <- NULL
} else {
  sde_monocle <- switchde(sce, monocle_pseudotime) %>% 
    mutate(algorithm = "monocle")
}

if(any(is.na(pc1_pseudotime))) {
  sde_pc1 <- NULL
} else {
  sde_pc1 <- switchde(sce, pc1_pseudotime) %>% 
    mutate(algorithm = "pc1")
}

if(any(is.na(dpt_pseudotime))) {
  sde_dpt <- NULL
} else {
  sde_dpt <- switchde(sce, dpt_pseudotime) %>% 
    mutate(algorithm = "dpt")
}



sde <- bind_rows(sde_ouija, sde_tscan, sde_monocle, sde_pc1, sde_dpt) %>% 
  select(algorithm, gene, pval)

output_file <- paste0(paste("data/marker-vs-txome/subset_de_fits/sde_fit", dataset, rep, nam_index, sep = "_"), ".csv")

write_csv(sde, output_file)





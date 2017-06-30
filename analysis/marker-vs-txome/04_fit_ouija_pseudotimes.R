set.seed(1234L)

library(tidyverse)
library(scater)

## Ouija consistency using HMC

source("scripts/pseudotime_fitting_functions.R")

sce_list_with_pseudotime <- readRDS("data/scesets/scesets_with_pseudotime.rds")
markers_to_use <- readRDS("data/marker-vs-txome/markers_to_use.rds")

# number of additional markers
n_additional_markers <- c(1,5,10,20,50,100,500,1000)


## Find what we're fitting

args <- commandArgs(trailingOnly = TRUE)

i <- as.numeric(args[1]) # which dataset
rep <- as.numeric(args[2]) # which rep
nam_index <- as.numeric(args[3]) # number of additional markers

sce <- sce_list_with_pseudotime[[i]]
genes_to_use <- c(which(fData(sce)$is_marker),
                  markers_to_use[[i]][[rep]][[nam_index]])
sce_reduced <- sce[genes_to_use,]


## Need to set the ouija priors
is_marker_gene <- fData(sce)$is_marker
k_mean <- fData(sce)$k_mean[is_marker_gene]
k_mean <- c(k_mean, rep(0, n_additional_markers[nam_index]))

# ouija_pseudotime <- fit_ouija_pseudotime(sce_reduced, inference_type = "hmc",
#                                         k_mean = k_mean, marker_only = FALSE,
#                                         iter = 6000)

ouija_pseudotime <- fit_ouija_pseudotime(sce_reduced, inference_type = "vb",
                                        k_mean = k_mean, marker_only = FALSE)


oui_df <- data_frame(algorithm = "ouija", 
                     dataset = names(sce_list_with_pseudotime)[i],
                     rep = rep, n_markers = n_additional_markers[nam_index],
                     pseudotime = ouija_pseudotime)

output_filename <- paste("ouija_pseudotime", i, rep, nam_index, sep = "_")
output_filename <- file.path("data", "marker-vs-txome", "ouija_fits", paste0(output_filename, ".csv"))

write_csv(oui_df, path = output_filename)







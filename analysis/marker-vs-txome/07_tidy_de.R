library(tidyverse)
library(stringr)
library(magrittr)


sce_list_with_marker_pseudotime <- readRDS("data/scesets/scesets_with_marker_pseudotime.rds")

datasets <- names(sce_list_with_marker_pseudotime)
rm(sce_list_with_marker_pseudotime)

alpha <- 0.05


## Parse the txwide and marker results
txwide_marker <- read_csv("txwide_marker_de.csv")

pseudotime_str <- c(# "monocle_pseudotime", 
                    "tscan_pseudotime",
                    "ouija_pseudotime", "pc1_pseudotime",
                    "dpt_pseudotime", "pc1_marker_pseudotime",
                    "tscan_marker_pseudotime", # "monocle_marker_pseudotime",
                    "dpt_marker_pseudotime")

pseudotime_comparisons <- list(
  # monocle = c("monocle_pseudotime", "monocle_marker_pseudotime"),
  tscan = c("tscan_pseudotime", "tscan_marker_pseudotime"),
  ouija = c("ouija_pseudotime"),
  pc1 = c("pc1_pseudotime", "pc1_marker_pseudotime"),
  dpt = c("dpt_pseudotime", "dpt_marker_pseudotime")
)

extract_sig <- function(df, dset, pst) {
  qval <- filter(df, dataset == dset & pst_str == pst) %>% 
    extract2("pval") %>% p.adjust(method = "BH")
  qval < alpha
}

sig_calls <- lapply(datasets, function(dset) {
  sub_list <- lapply(pseudotime_str, function(pst_str) {
    extract_sig(txwide_marker, dset, pst_str)
  })
  names(sub_list) <- pseudotime_str
  sub_list
})
names(sig_calls) <- datasets


## Parse the individual files

sde_files <- dir("switchde_fits")
sde_files <- sapply(sde_files, function(x) gsub(".csv", "", x, fixed = TRUE))

n_files <- length(sde_files)

name_mat <- str_split_fixed(sde_files, "_", Inf)

datasets <- as.numeric(name_mat[,3])
reps <- as.numeric(name_mat[,4])
nam_indices <- as.numeric(name_mat[,5])





## We're going to do this one file at a time to avoid RAM overflows

dset_list <- lapply(seq_len(n_files), function(n) {
  dataset <- datasets[n]
  pvals_n <- read_csv(file.path("switchde_fits", paste0(sde_files[n], ".csv")))
  
  sdfs <- lapply(names(pseudotime_comparisons), function(alg) {
    p_vals <- filter(pvals_n, algorithm == alg) %>% 
      extract2("pval") %>% 
      p.adjust(method = "BH")
    is_sig <- p_vals < alpha
    n_sig <- sum(is_sig)  
    
    tpr <- sapply(pseudotime_comparisons[[alg]], function(pst_str) {
      overlap <- sum(is_sig & sig_calls[[dataset]][[pst_str]])
      overlap / n_sig
    })
    data_frame(algorithm = alg, tpr = tpr, pst_str = names(tpr))
  })
  bind_rows(sdfs) %>% 
    mutate(dataset = dataset, rep = reps[n],
           nam_index = nam_indices[n])
})

bind_rows(dset_list) %>% 
  write_csv("tidy_de.csv")



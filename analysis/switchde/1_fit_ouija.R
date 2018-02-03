
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(aargh))

devtools::load_all("~/oxford/bnlfa/oui/ouija")

set.seed(123L)

fit_ouija_traces <- function(output_csv = "output.csv",
                             traces_to_sample = 100) {
  sce_list <- readRDS("../../data/scesets_with_marker_pseudotime.rds")
  sce <- sce_list[['trapnell']]
  sce_marker <- sce[fData(sce)$is_marker, ]
  
  
  oui <- ouija(sce_marker, iter = 3000, 
               switch_strength_sd = rep(5, 5))
  
  pst <- rstan::extract(oui$fit, "t")$t
  n_iter <- nrow(pst)
  pst_sampled <- pst[sample(n_iter, traces_to_sample),]
  
  # Transpose and tidy
  trace_df <- as_data_frame(t(pst_sampled))
  names(trace_df) <- seq_len(traces_to_sample)
  trace_df <- mutate(trace_df, cell = sampleNames(sce_marker)) %>% 
    gather(sample, pseudotime, -cell)
  
  write_csv(trace_df, output_csv)
}

aargh(fit_ouija_traces)
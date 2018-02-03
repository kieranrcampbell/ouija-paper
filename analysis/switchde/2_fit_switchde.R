  

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(aargh))
suppressPackageStartupMessages(library(switchde))


fit_switchde <- function(input_csv = "input.csv", 
                         output_csv = "output.csv",
                         trace_sample = 1) {
  pst_df <- read_csv(input_csv)
  pst_df <- filter(pst_df, sample == trace_sample)
  
  
  sce_list <- readRDS("../../data/scesets_with_marker_pseudotime.rds")
  sce <- sce_list[['trapnell']]
  
  stopifnot(all.equal(pst_df$cell, sampleNames(sce)))
  
  pseudotime <- pst_df$pseudotime
  
  sde <- switchde(sce, pseudotime)
  
  sde <- mutate(sde, trace_sample = trace_sample,
                gsn = fData(sce)$gene_short_name)
  
  write_csv(sde, output_csv)

  
}

aargh(fit_switchde)



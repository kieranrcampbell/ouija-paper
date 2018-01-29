

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(aargh))

devtools::load_all("~/oxford/bnlfa/oui/ouija")


jk_ouija <- function(gene1 = 1, gene2 = 0, output_csv = "output.csv") {
  if(gene1 == gene2) {
    write_csv(data_frame(), output_csv)
    return(NULL)
  }
  genes_to_remove <- gene1
  if(gene2 > 0) {
    genes_to_remove <- c(genes_to_remove, gene2)
  }
  print(genes_to_remove)
  sce_list <- readRDS("../../data/scesets_with_marker_pseudotime.rds")
  sce <- sce_list[['trapnell']]
  sce_marker <- sce[fData(sce)$is_marker, ]
  
  sce_marker <- sce_marker[-genes_to_remove, ]
  
  oui <- ouija(sce_marker, iter = 3000, 
               switch_strength_sd = rep(5, 5 - length(genes_to_remove)))
  
  pst <- map_pseudotime(oui)
  
  df <- data_frame(cell = sampleNames(sce_marker), pseudotime = pst)
  write_csv(df, output_csv)
}

aargh(jk_ouija)

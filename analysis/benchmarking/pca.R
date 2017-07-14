## Go through the synthetic pseudotimes and see how PCA 
## corresponds to the truth
library(rhdf5)
library(dplyr)
library(VGAM)
library(readr)

pca_pseudotime <- function(d) {
  Y <- select(d, -pseudotime)
  tpca <- prcomp(Y)$x[,1]
  return( c(abscor = abs(cor(d$pseudotime, tpca)),
               kendall_tau = abs(kendall.tau(d$pseudotime, tpca)),
               sign = sign(cor(d$pseudotime, tpca))) )
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

h5_file <- input_file

struct <- h5ls(h5_file, recursive = FALSE)
Gs <- sapply(struct$name, function(x) as.numeric(gsub("G", "", x)))

pca <- lapply(Gs, function(G) {
  group_name <- paste0("G", G)
  dfs <- h5read(h5_file, group_name)
  dfs_d <- lapply(dfs, `[[`, "d")
  pca_res <- data.frame(t(sapply(dfs_d, pca_pseudotime)))
  pca_res$d <- rownames(pca_res)
  pca_res$G <- G
  pca_res
})

pca_df <- do.call("rbind", pca)
pca_df$G <- as.factor(pca_df$G)

saveRDS(pca_df, file = output_file)


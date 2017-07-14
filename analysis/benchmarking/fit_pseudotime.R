library(scater)

library(ggplot2)
library(readr)
library(dplyr)
library(rhdf5)
library(VGAM)


get_pseudotime <- function(sce, algorithm = c("dpt", "monocle", "tscan")) {
  algorithm <- match.arg(algorithm)
  if(algorithm == "dpt") {
    return(fit_dpt_pseudotime(sce))
  }
  if(algorithm == "monocle") {
    return(fit_monocle_pseudotime(sce))
  }
  if(algorithm == "tscan") {
    return(fit_tscan_pseudotime(sce))
  }
}

set.seed(123L)

args <- commandArgs(trailingOnly = TRUE)
algorithm <- args[1]
input_file <- args[2]
output_file <- args[3]

stopifnot(algorithm %in% c("dpt", "monocle", "tscan"))

h5str <- h5ls(input_file)

Gs <- c(6,9,12,15)
d <- 1:40

cors <- matrix(0, nrow = 0, ncol = 5)

for(g in Gs) {
  for(r in d) {
    Y <- h5read(input_file, paste0("G", g, "/", r))$d
    X <- Y[,1:g]
    sce <- newSCESet(exprsData = t(X))
    true_pst <- Y[,g+1]
    pst <- get_pseudotime(sce, algorithm)
    abscor <- abs(cor(true_pst, pst))
    sgn <- sign(cor(true_pst, pst))
    ktau <- abs(VGAM::kendall.tau(true_pst, pst))
    cors <- rbind(cors, c(abscor, ktau, sgn, r, g))
  }
}

colnames(cors) <- c("abscor", "kendall_tau", "sign", "d", "G")
cors <- data.frame(cors)
cors$condition <- algorithm

saveRDS(cors, file = output_file)

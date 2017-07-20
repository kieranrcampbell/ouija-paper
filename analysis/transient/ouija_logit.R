
library(dplyr)
library(aargh)
library(readr)

library(ouija)

bnlfa_pseudotime <- function(d) {
  Y <- as.matrix(select(d, -pseudotime))
  G <- ncol(Y)
    
  k_means <- rep(0, G)
  t0_means <- rep(0.5, G)
  strength_sd <- time_sd <- rep(5, G)

  bm <- ouija(Y, strengths = k_means, times = t0_means, 
              strength_sd = strength_sd, time_sd = time_sd)
  
  tmap <- map_pseudotime(bm)
  return( c(abscor = abs(cor(d$pseudotime, tmap)),
               kendall_tau = abs(cor(d$pseudotime, tmap, method = "kendall")),
               sign = sign(cor(d$pseudotime, tmap))) )
}


benchmark_ouija <- function(rep = 1,
                            G = 8,
                            prop_switch = 0.25,
                            input_file = "input_file",
                            output_file = "output_file") {

  set.seed(123)
  
  dfs <- read_csv(input_file)
  bnlfa_results <- bnlfa_pseudotime(dfs)
  bnlfa_res <- data.frame(t(bnlfa_results))
  bnlfa_res$G <- G
  bnlfa_res$prop_switch <- prop_switch
  bnlfa_res$condition <- condition
  bnlfa_res$rep <- rep

  write_csv(bnlfa_res, output_file)
}


aargh(benchmark_ouija)
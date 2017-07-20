
library(rhdf5)
library(dplyr)
library(aargh)
library(readr)

library(ouija)

bnlfa_pseudotime <- function(d, condition = "noninformative") {
  Y <- as.matrix(select(d, -pseudotime))
  
  if(condition == "true") {
    k_means <- k
    t0_means <- t0
    strength_sd <- time_sd <- rep(0.1, length(k))
  } else if(condition == "noninformative") {
    k_means <- rep(0, length(k))
    t0_means <- rep(0.5, length(k))
    strength_sd <- time_sd <- rep(5, length(k))
  } else if(condition == "t0_uncertainty") {
    k_means <- k
    t0_means <- t0 + rnorm(length(t0), 0, 0.1)
    strength_sd <- time_sd <- rep(0.1, length(k))
  } else if(condition == "t0_midpoint") {
    k_means <- k
    t0_means <- rep(0.5, length(k))
    strength_sd <- time_sd <- rep(0.1, length(k))
  }

  bm <- ouija(Y, strengths = k_means, times = t0_means, 
              strength_sd = strength_sd, time_sd = time_sd,
              normalise_expression = (condition == "noninformative"))
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
  bnlfa_results <- bnlfa_pseudotime(dfs$d, dfs$k, dfs$t0, condition)
  bnlfa_res <- data.frame(t(bnlfa_results))
  bnlfa_res$G <- G
  bnlfa_res$condition = condition
  bnlfa_res$rep = rep

  write_csv(bnlfa_res, output_file)
}


aargh(benchmark_ouija)
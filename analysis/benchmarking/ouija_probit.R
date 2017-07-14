
library(rhdf5)
library(dplyr)
library(aargh)
library(readr)

library(ouija)

bnlfa_pseudotime <- function(d, k_sgn, t0, condition) {
  Y <- as.matrix(select(d, -pseudotime))
  activation_guess <- 50 # guess of k for noninformative runs
  
  if(condition == "true") {
    k_means <- k_sgn * activation_guess
    t0_means <- t0
    strength_sd <- time_sd <- rep(0.1, length(k_sgn))
  } else if(condition == "noninformative") {
    k_means <- rep(0, length(k_sgn))
    t0_means <- rep(0.5, length(k_sgn))
    strength_sd <- time_sd <- rep(5, length(k_sgn))
  }

  bm <- ouija(Y, strengths = k_means, times = t0_means, 
              strength_sd = strength_sd, time_sd = time_sd,
              normalise_expression = (condition == "noninformative"))
  tmap <- map_pseudotime(bm)
  return( c(abscor = abs(cor(d$pseudotime, tmap)),
               kendall_tau = abs(cor(d$pseudotime, tmap, method = "kendall")),
               sign = sign(cor(d$pseudotime, tmap))) )
}


benchmark_ouija <- function(condition = "noninformative",
                            rep = 1,
                            G = 6,
                            input_file = "input_file",
                            output_file = "output_file") {
  stopifnot(condition %in% c("true", "noninformative"))

  
  
  set.seed(123)
  
  group_name <- paste0("G", G, "/", rep)
  dfs <- h5read(input_file, group_name)

  bnlfa_results <- bnlfa_pseudotime(dfs$d, dfs$k_sgn, dfs$t0, condition)
  bnlfa_res <- data.frame(t(bnlfa_results))
  bnlfa_res$G <- G
  bnlfa_res$condition = condition
  bnlfa_res$rep = rep
  
  write_csv(bnlfa_res, output_file)
}


aargh(benchmark_ouija)
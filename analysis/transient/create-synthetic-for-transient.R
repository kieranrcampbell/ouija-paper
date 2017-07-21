
library(ggplot2)
library(readr)


source("analysis/transient/common_simulation_functions.R")



#' For C cells and G genes
generate_dataset <- function(N, G, prop_switch) {
  G_switch <- floor(G * prop_switch)
  G_transient <- G - G_switch
  true_pst <- runif(N)
  
  t0 <- runif(G_switch, 0.1, 0.9)
  mu0 <- runif(G_switch, 3, 4)
  k <- runif(G_switch, 5, 20) * sample(c(-1, 1), G_switch, replace = TRUE)
  
  Y_switch <- mapply(rsinglecell_sigmoid, k, mu0, t0, MoreArgs = list(pst = true_pst))
  colnames(Y_switch) <- paste0("switch_gene_", seq_len(G_switch))
  
  mu0 <- runif(G_transient, 3, 4)
  p <- runif(G_transient, 0.3, 0.7)
  b <- runif(G_transient, 1, 5)
  Y_transient <- mapply(rsinglecell_transient, mu0, b, p, MoreArgs = list(pst = true_pst))
  colnames(Y_transient) <- paste0("transient_gene_", seq_len(G_transient))
  
  Y <- cbind(Y_switch, Y_transient)
  
  return( list(Y = data.frame(Y, pseudotime = true_pst), G = G, prop_switch = prop_switch) )
}


set.seed(123L)

N <- 100 # cells
Gs <- c(8, 12, 16, 24) # genes
prop_switches <- c(0.5, 0.75)
N_rep <- 100 # number of datasets at each 'condition

base_filename <- "data/transient/synthetic/transient_sim_"

for(G in Gs) {
  for(ps in prop_switches) {
    for(rep in seq_len(N_rep)) {
      fname <- paste0(base_filename, G, "_", ps, "_", rep, ".csv")
      d <- generate_dataset(N, G, ps)
      write_csv(d$Y, fname)
    }
  }
}


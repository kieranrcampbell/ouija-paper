
library(ggplot2)
library(rhdf5)


source("analysis/benchmarking/common_simulation_functions.R")

sigmoid <- function(pst, k, mu0, t0) {
  2 * mu0 / (1 + exp(-k * (pst - t0)))
}

transient <- function(pst, mu0, b, p) {
  2 * mu0 * exp(-10 * b(pst - p)^2)
}


#' Synthetic single-cells with mean and dispersion
rsinglecell_sigmoid <- function(k, mu0, t0, pst) {
  mean <- pstmean(pst, k, mu0, t0)
  x <- rcensnorm(mean, sqrt(gvar(mean)))
  
  pdrop <- pdropout(mean)
  is_dropout <- rbernoulli(pdrop)
  x[is_dropout] <- 0
  return( x )
}

rsinglecell_transient <- function(mu0, b, p, pst) {
  mean <- transient(pst, mu0, b, p)
  x <- rcensnorm(mean, sqrt(gvar(mean)))
  
  pdrop <- pdropout(mean)
  is_dropout <- rbernoulli(pdrop)
  x[is_dropout] <- 0
  return( x )
}

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
  Y_transient <- mapply(rsinglecell_transient, mu0, p, b, MoreArgs = list(pst = true_pst))
  colnames(Y_transient) <- paste0("switch_gene_", seq_len(G_transient))
  
  Y <- rbind(Y_switch, Y_transient)
  
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
    for(rep in N_rep) {
      fname <- paste0(base_filename, "_", G, "_", ps, "_", rep, ".csv")
      d <- generate_dataset(N, G, ps)
      write_csv(d$Y, fname)
    }
  }
}


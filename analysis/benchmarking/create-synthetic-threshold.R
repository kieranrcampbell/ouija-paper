
library(ggplot2)
library(rhdf5)

source("analysis/benchmarking/common_simulation_functions.R")

threshold <- function(pst, t0, mu0, k_sgn) {
  m <- runif(1, 1, 2) * k_sgn
  c <- -m * t0
  y_star <- rnorm(length(pst), m * pst + c, 0.1)
  y <- rep(0, length(pst))
  y[y_star > 0] <- 2 * mu0
  return(y)
}

#' Synthetic single-cells with mean and dispersion
rsinglecell <- function(k_sgn, mu0, t0, pst) {
  mean <- threshold(pst, t0, mu0, k_sgn)
  x <- rcensnorm(mean, sqrt(gvar(mean)))
  
  pdrop <- pdropout(mean)
  is_dropout <- rbernoulli(pdrop)
  x[is_dropout] <- 0
  return( x )
}

#' For C cells and G genes
generate_dataset <- function(C, G) {
  true_pst <- runif(C)
  
  t0 <- runif(G, 0.1, 0.9)
  mu0 <- runif(G, 3, 4)
  k_sgn <- sample(c(-1, 1), G, replace = TRUE)
  
  Y <- mapply(rsinglecell, k_sgn, mu0, t0, MoreArgs = list(pst = true_pst))
  
  return( list(d = data.frame(Y, pseudotime = true_pst), k_sgn = k_sgn, t0 = t0, mu0 = mu0) )
}



set.seed(123L)

N <- 100 # cells
Gs <- c(6, 9, 12, 15) # genes
N_rep <- 40 # number of datasets at each 'condition

h5_file <- "data/benchmarking/threshold_synthetic.h5"
if(file.exists(h5_file)) file.remove(h5_file)
h5createFile(h5_file)

for(G in Gs) {
  group_name <- paste0("G", G)
  h5createGroup(h5_file, group_name)
  for(i in seq_len(N_rep)) {
    obj_name <- paste0(group_name, "/", i)
    h5createGroup(h5_file, obj_name)
    d <- generate_dataset(N, G)
    h5write(d$d, h5_file, paste0(obj_name, "/d"))
    h5write(d$k_sgn, h5_file, paste0(obj_name, "/k_sgn"))
    h5write(d$t0, h5_file, paste0(obj_name, "/t0"))
    h5write(d$mu0, h5_file, paste0(obj_name, "/mu0"))
  }
}

G <- 12
d <- generate_dataset(N, G)
dfy <- d$d
names(dfy) <- c(paste0("Gene", 1:G), "Pseudotime")
dfym <- reshape2::melt(dfy, id.vars = "Pseudotime", 
                       variable.name = "Gene",
                       value.name = "Expression")


gene_plt <- ggplot(dfym, aes(x = Pseudotime, y = Expression)) +
  geom_point() + facet_wrap(~ Gene, nrow = 3, scales = "free_y") + 
  stat_smooth(colour = 'red', se = FALSE) + theme_bw()

ggsave(gene_plt, file = "figs/example_genes_threshold.png", width=8, height=4)

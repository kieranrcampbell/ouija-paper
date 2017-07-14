
library(ggplot2)
library(rhdf5)


source("analysis/benchmarking/common_simulation_functions.R")

pstmean <- function(pst, k, mu0, t0) {
  return( 2 * mu0 / (1 + exp(-k * (pst - t0))))
}

#' Synthetic single-cells with mean and dispersion
rsinglecell <- function(k, mu0, t0, pst) {
  mean <- pstmean(pst, k, mu0, t0)
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
  k <- runif(G, 5, 20) * sample(c(-1, 1), G, replace = TRUE)
  
  Y <- mapply(rsinglecell, k, mu0, t0, MoreArgs = list(pst = true_pst))
  
  return( list(d = data.frame(Y, pseudotime = true_pst), k = k, t0 = t0, mu0 = mu0) )
}


set.seed(123L)

N <- 100 # cells
Gs <- c(6, 9, 12, 15) # genes
N_rep <- 40 # number of datasets at each 'condition

h5_file <- "data/benchmarking/logit_synthetic.h5"
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
    h5write(d$k, h5_file, paste0(obj_name, "/k"))
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

ggsave(gene_plt, file = "figs/example_genes_logit.png", width=8, height=4)

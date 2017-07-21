beta_0_param <- 1.763
beta_1_param <- -1.156
disp_param <- 3.5

pdropout <- function(x, beta_0 = beta_0_param, beta_1 = beta_1_param) {
  mu <- beta_0 + beta_1 * x
  return(1 / (1 + exp(-mu)))
}

gvar <- function(mean_expr) {
  return (disp_param * mean_expr)
}

rcensnorm <- function(mean, sd, lower_bound = 0) {
  x <- rnorm(length(mean), mean = mean, sd = sd)
  x[x < lower_bound] <- lower_bound
  return( x )
}

rbernoulli <- function(pp){
  sapply(pp, function(p) sample(0:1, 1, prob = c(p, 1 - p)))
}


sigmoid <- function(pst, k, mu0, t0) {
  2 * mu0 / (1 + exp(-k * (pst - t0)))
}

transient <- function(pst, mu0, b, p) {
  2 * mu0 * exp(-10 * b * (pst - p)^2)
}


#' Synthetic single-cells with mean and dispersion
rsinglecell_sigmoid <- function(k, mu0, t0, pst) {
  mean <- sigmoid(pst, k, mu0, t0)
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
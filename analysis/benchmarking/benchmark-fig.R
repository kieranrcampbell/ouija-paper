
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
library(cowplot)
library(readr)

get_sim_df <- function(root_dir, sim_str = "Latent probit") {
  
  files <- dir(root_dir)
  
  pca_file <- files[grep("pca", files)]
  
  ## pca results
  pca_df <- readRDS(file.path(root_dir, pca_file)) %>% 
    tbl_df()
  pca_df$condition <- "PCA"
  pca_df <- rename(pca_df, rep = d)
  
  ## bnlfa results
  ouija_files <- files[grep("csv", files)]
  ouija_list <- lapply(ouija_files, function(f) {
    dat <- read_csv(file.path(root_dir, f))
  })
  df <- do.call("rbind", ouija_list)
  df <- rbind(df, pca_df)
  
  
  
  dpt_cors <- readRDS(file.path(root_dir, "dpt.rds")) %>% 
    tbl_df()
  dpt_cors <- rename(dpt_cors, rep = d)
  df <- rbind(df, dpt_cors)
  
  # ## Monocle results
  # mfiles <- load("data/monocle.Rdata")
  # 
  # df <- rbind(df, monocle_cors)
  
  condition_levels <- c("PCA",  "DPT", "noninformative", "true", "t0_uncertainty", "t0_midpoint")
  alg_names <- c("PCA", "DPT", "Ouija noninformative ", "Ouija informative", "Ouija switch uncertainty", "Ouija switch midpoint")
  df$G <- factor(as.numeric(df$G))
  
  df$condition <- factor(df$condition, levels = condition_levels, ordered = TRUE)
  df$sim_regime <- sim_str
  return(df)
}

sim_dfs <- Map(
  get_sim_df,
  paste0("../", c("cloglog-benchmarking", "logit-benchmarking", "probit-benchmarking", "threshold-benchmarking"), "/data/"),
  c("Complementary log-log", "Sigmoidal", "Probit", "Threshold")
)
sim_df_all <- bind_rows(sim_dfs)

sim_df <- filter(sim_df_all, condition != "t0_uncertainty", condition != "t0_midpoint")

sim_df$sim_regime <- fct_relevel(sim_df$sim_regime,
                                 c("Sigmoidal", "Complementary log-log", "Probit", "Threshold"))

condition_levels <- c("PCA",  "DPT", "noninformative", "true")
alg_names <- c("PCA", "DPT", "Ouija noninformative ", "Ouija informative")

ggplot(sim_df, aes(x = G, y = abscor, fill = condition, color = condition)) + 
  geom_boxplot(outlier.shape = NA) +
  xlab("Number of genes modelled") +
  scale_fill_brewer(palette = "Set1", name = "Algorithm", 
                    breaks = condition_levels,
                    labels = alg_names) +
  scale_color_brewer(palette = "Set1", name = "Algorithm", 
                     breaks = condition_levels,
                     labels = alg_names) +
  cowplot::theme_cowplot(font_size = 11) + 
  ylab(expression("Pearson" ~ rho)) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold")) +
    facet_wrap(~ sim_regime, nrow = 1) +
  ylim(0.4, NA)
boxplt <- last_plot()

# Mean function plots -----------------------------------------------------

cloglog_mean <- function(pst, k = 10, t0 = 0.5, mu0 = 5) {
  y <- 1 - exp(-exp(k * (pst - t0)))
  return(2 * mu0 * y)
}
threshold <- function(pst, t0, mu0, k_sgn) {
  m <- runif(1, 1, 2) * k_sgn
  c <- -m * t0
  y_star <- rnorm(length(pst), m * pst + c, 0.1)
  y <- rep(0, length(pst))
  y[y_star > 0] <- 2 * mu0
  return(y)
}
pstmean <- function(pst, k, mu0, t0) {
  return( 2 * mu0 / (1 + exp(-k * (pst - t0))))
}
probit <- function(pst, t0, mu0, k) {
  2 * mu0 * pnorm(k * (pst - t0))
}


k <- 10; mu0 <- 5; t0 <- 0.5

pst <- seq(from = 0, to = 1, length.out = 500)

mean_df <- data_frame(
  pst = pst,
  Sigmoidal = pstmean(pst, k, mu0, t0),
  Probit = probit(pst, t0, mu0, k),
  Threshold = threshold(pst, t0, mu0, sign(k)),
  "Complementary log-log" = cloglog_mean(pst, k, t0, mu0)
)

mean_df_tidy <- gather(mean_df, fn, value, -pst)
mean_df_tidy$fn <-  fct_relevel(mean_df_tidy$fn,
                                c("Sigmoidal", "Complementary log-log", "Probit", "Threshold"))

ggplot(mean_df_tidy, aes(x = pst, y = value)) +
  geom_point(alpha = 0.5, size = 1) + 
  facet_wrap(~ fn, nrow = 1) +
  theme_cowplot(font_size = 11) +
  xlab(paste("Pseudotime", sprintf('\u2192'))) + ylab("Mean expression") +
  theme(strip.background = element_blank(),
      strip.text = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(face = 'bold'),
      axis.line = element_blank(),
      axis.text.y = element_text(color = 'white'),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = "grey92")) +
  scale_x_continuous(position = "top")

mean_plot <- last_plot()

plot_grid(mean_plot, boxplt, ncol = 1, rel_heights = c(3,5))

ggsave("figs/fig4.png", width = 8, height = 6)



# ggplot(mean_df_tidy, aes(x = pst, y = value)) +
#   geom_point(alpha = 0.5) + 
#   facet_wrap(~ fn) +
#   theme_cowplot(font_size = 11) +
#   xlab(paste("Pseudotime", sprintf('\u2192'))) + ylab("Mean value") +
#   theme(strip.background = element_rect(fill = "grey90"),
#         strip.text = element_text(face = "bold"))

## Save summary statistics to file

sum_res <- group_by(sim_df, G, condition, sim_regime) %>% 
  summarise(mean_cor = mean(abscor)) %>% 
  arrange(sim_regime)

write_csv(sum_res, "data/benchmarking/summarised_results.csv")



# t0 uncertainty plot -----------------------------------------------------


condition_levels <- c("PCA",  "DPT", "noninformative", "t0_midpoint", "t0_uncertainty", "true")
alg_names <- c("PCA", "DPT", "Ouija noninformative",
               "Ouija switch midpoint", "Ouija switch uncertainty", "Ouija informative")

sim_df_logit <- filter(sim_df_all, sim_regime == "Sigmoidal")
sim_df_logit$condition <- fct_relevel(sim_df_logit$condition, condition_levels)

ggplot(sim_df_logit, aes(x = G, y = abscor, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Number of genes modelled") +
  scale_fill_brewer(palette = "Set1", name = "Algorithm",
                    breaks = condition_levels,
                    labels = alg_names) +
  scale_color_brewer(palette = "Set1", name = "Algorithm",
                     breaks = condition_levels,
                     labels = alg_names) +
  cowplot::theme_cowplot(font_size = 11) +
  ylab(expression("Pearson" ~ rho)) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold")) +
  # facet_wrap(~ sim_regime, nrow = 1) +
  ylim(0.5, NA)
boxplt2 <- last_plot()

ggsave("figs/logit-only-ouija.png", width = 6, height = 5)


library(scater)
library(tidyverse)
library(GGally)

sce_list_with_marker_pseudotime <- readRDS("data/scesets/scesets_with_marker_pseudotime.rds")
cmats <- readRDS("data/marker-vs-txome/pseudotimes_across_marker_subset_array.rds")

n_additional_markers <- c(1,5,10,20,50,100,500,1000)
n_reps <- 50

theme_set(theme_bw())

names(sce_list_with_marker_pseudotime) <- plyr::mapvalues(names(sce_list_with_marker_pseudotime),
                                  to = c("Trapnell", "Shin", "Zhou"),
                                  from = c("trapnell", "shin", "hsc"))

## Swap out ouija VB for ouija HMC

for(i in seq_along(sce_list_with_marker_pseudotime)) {
  for(rep in seq_len(n_reps)) {
    for(nam_index in seq_along(n_additional_markers)) {
      ouija_file <- paste("ouija_pseudotime", i, rep, nam_index, sep = "_")
      ouija_file_path <- file.path("data", "marker-vs-txome", "ouija_fits", paste0(ouija_file, ".csv"))
      ouija_hmc <- read_csv(ouija_file_path)
      cmats[[i]][rep, nam_index, 3,] <- ouija_hmc$pseudotime
    }
  }
}


## Make correlation matrix

cmat_list <- lapply(seq_along(cmats), function(i) {
  sce <- sce_list_with_marker_pseudotime[[i]]
  sdf <- data_frame(dataset = character(), algorithm = character(),
                    rep = numeric(), n_markers = numeric(), cor = numeric())
  class(cmats[[i]]) <- "numeric"
  for(rep in seq_len(n_reps)) {
    for(n_markers in seq_along(n_additional_markers)) {
      # monocle_correlation <- cor(sce$monocle_marker_pseudotime,
      #                            cmats[[i]][rep, n_markers, 1, ])
      tscan_correlation <- cor(sce$tscan_marker_pseudotime,
                                 cmats[[i]][rep, n_markers, 2, ], use = "na")
      ouija_correlation <- cor(sce$ouija_pseudotime,
                                 cmats[[i]][rep, n_markers, 3, ])
      pc1_correlation <- cor(sce$pc1_marker_pseudotime,
                                 cmats[[i]][rep, n_markers, 4, ])
      dpt_correlation <- cor(sce$dpt_marker_pseudotime,
                             cmats[[i]][rep, n_markers, 5, ])
      sdf_ind <- data_frame(dataset = names(sce_list_with_marker_pseudotime)[i],
                            algorithm = c("TSCAN", "Ouija", "PC1", "DPT"),
                            rep, n_markers = n_additional_markers[n_markers],
                            cor = abs(c(tscan_correlation,
                                    ouija_correlation, pc1_correlation, dpt_correlation)))
      sdf <- rbind(sdf, sdf_ind)
    }
  }
  return(sdf)
})

c_tidy <- bind_rows(cmat_list)
c_tidy$n_markers <- as.factor(c_tidy$n_markers)

c_tidy$algorithm <- factor(c_tidy$algorithm, levels = c("Ouija", "PC1",
                                                        "TSCAN", "DPT"))

ggplot(data = c_tidy, aes(x = n_markers, y = cor, color = algorithm, fill = algorithm)) +
  geom_boxplot() + facet_wrap(~ dataset) +
  xlab("Number of additional genes") +
  ylab("Correlation to marker pseudotime") +
  scale_fill_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm")

marker_plot <- last_plot()

# marker_plot_line <- group_by(c_tidy, dataset, algorithm, n_markers) %>% 
#   summarise(mean_cor = mean(cor), lower_5 = quantile(cor, 0.05), 
#             upper_95 = quantile(cor, 0.95)) %>% 
#   ggplot(aes(x = n_markers, group = algorithm, y = mean_cor, color = algorithm)) +
#   geom_line() + facet_wrap(~ dataset) +
#   xlab("Number of additional genes") +
#   ylab("Correlation to marker pseudotime") +
#   scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") 




# To global pseudotime ----------------------------------------------------


gmat_list <- lapply(seq_along(cmats), function(i) {
  sce <- sce_list_with_marker_pseudotime[[i]]
  sdf <- data_frame(dataset = character(), algorithm = character(),
                    rep = numeric(), n_markers = numeric(), cor = numeric())
  for(rep in seq_len(n_reps)) {
    for(n_markers in seq_along(n_additional_markers)) {
      # monocle_correlation <- cor(sce$monocle_pseudotime,
      #                            cmats[[i]][rep, n_markers, 1, ])
      
      if(all(is.na(cmats[[i]][rep, n_markers, 2, ])) | !is.numeric(cmats[[i]][rep, n_markers, 2, ])) {
        tscan_correlation <- NA
      } else {
        tscan_correlation <- cor(sce$tscan_pseudotime,
                                 cmats[[i]][rep, n_markers, 2, ], use = "na")
      }
      ouija_correlation <- cor(sce$ouija_pseudotime,
                               cmats[[i]][rep, n_markers, 3, ])
      pc1_correlation <- cor(sce$pc1_pseudotime,
                             cmats[[i]][rep, n_markers, 4, ])
      dpt_correlation <- cor(sce$dpt_pseudotime,
                             cmats[[i]][rep, n_markers, 5, ])
      sdf_ind <- data_frame(dataset = names(sce_list_with_marker_pseudotime)[i],
                            algorithm = c("TSCAN", "Ouija", "PC1", "DPT"),
                            rep, n_markers = n_additional_markers[n_markers],
                            cor = abs(c(tscan_correlation,
                                        ouija_correlation, pc1_correlation, dpt_correlation)))
      sdf <- rbind(sdf, sdf_ind)
    }
  }
  return(sdf)
})

g_tidy <- bind_rows(gmat_list)
g_tidy$n_markers <- as.factor(g_tidy$n_markers)

g_tidy$algorithm <- factor(g_tidy$algorithm, levels = c("Ouija", "PC1", "TSCAN", "DPT"))
g_tidy <- filter(g_tidy, algorithm != "Ouija")

ggplot(data = g_tidy, aes(x = n_markers, y = cor, color = algorithm, fill = algorithm)) +
  geom_boxplot() + facet_wrap(~ dataset) +
  xlab("Number of additional genes") +
  ylab("Correlation to global pseudotime") +
  scale_fill_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm")


global_plot <- last_plot()

# global_plot_line <- group_by(g_tidy, dataset, algorithm, n_markers) %>% 
#   summarise(mean_cor = mean(cor), lower_5 = quantile(cor, 0.05), 
#             upper_95 = quantile(cor, 0.95)) %>% 
#   ggplot(aes(x = n_markers, group = algorithm, y = mean_cor, color = algorithm)) +
#   geom_line() + facet_wrap(~ dataset) +
#   xlab("Number of additional genes") +
#   ylab("Correlation to global pseudotime") +
#   scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") 



# Within-algorithm correlation --------------------------------------------

wmat_list <- lapply(seq_along(cmats), function(i) {
  # sce <- sce_list_with_marker_pseudotime[[i]]
  sdf <- data_frame(algorithm = character(), correlation = numeric(),
                    n_markers = numeric(), dataset = character())

  for(n_markers in seq_along(n_additional_markers)) {
    # order always monocle-tscan-ouija-pc1
    matrices <- lapply(2:5, function(alg) { # Change back to 1:5 if using Monocle
      cor(t(cmats[[i]][, n_markers, alg, ]), use = 'na')
      })
  
    cors <- sapply(matrices, function(mat) abs(mat[lower.tri(mat)]))    
    cors <- as_data_frame(cors)
    # names(cors) <- c("Monocle", "TSCAN", "Ouija", "PC1", "DPT")
    names(cors) <- c("TSCAN", "Ouija", "PC1", "DPT")
    
    sdf_ind <- gather(cors, algorithm, correlation) %>% 
      mutate(n_markers = n_additional_markers[n_markers], 
             dataset = names(sce_list_with_marker_pseudotime)[i])
  
    sdf <- rbind(sdf, sdf_ind)
  }
  
  return(sdf)
})



w_tidy <- bind_rows(wmat_list)
w_tidy$n_markers <- as.factor(w_tidy$n_markers)

# w_tidy$algorithm <- factor(w_tidy$algorithm, levels = c("Ouija", "PC1",
#                                                         "Monocle", "TSCAN", "DPT"))
w_tidy$algorithm <- factor(w_tidy$algorithm, levels = c("Ouija", "PC1",
                                                        "TSCAN", "DPT"))

ggplot(data = w_tidy, aes(x = n_markers, y = correlation, color = algorithm, fill = algorithm)) +
  geom_boxplot() + facet_wrap(~ dataset) +
  xlab("Number of additional genes") +
  ylab("Within algorithm correlation") +
  scale_fill_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") +
  scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm")

within_plot <- last_plot()

# within_plot_line <- group_by(w_tidy, dataset, algorithm, n_markers) %>% 
#   summarise(mean_cor = mean(correlation)) %>% 
#   ggplot(aes(x = n_markers, group = algorithm, y = mean_cor, color = algorithm)) +
#   geom_line() + facet_wrap(~ dataset) +
#   xlab("Number of additional genes") +
#   ylab("Within algorithm correlation") +
#   scale_color_brewer(palette = "Set1", drop = FALSE, name = "Algorithm") 


# Across number of marker correlation -------------------------------------

amat_list <- lapply(seq_along(cmats), function(i) {
  sce <- sce_list_with_marker_pseudotime[[i]]
  sdf <- data_frame(algorithm = character(), correlation = numeric(),
                    dataset = character())
  
  # order always monocle-tscan-ouija-pc1
  cor_matrices <- sapply(2:5, function(alg) {
    marray <- cmats[[i]][, , alg, ]
    xc <- apply(marray, 1, function(m) {
      cor_mat <- cor(t(m), use = "na")
      abs(cor_mat[lower.tri(cor_mat)])
    })
    as.vector(xc)
  })
  
  cor_matrices <- as_data_frame(cor_matrices)
  # names(cor_matrices) <- c("monocle", "tscan", "ouija", "pc1", "dpt")
  names(cor_matrices) <- c("tscan", "ouija", "pc1", "dpt")
  
  sdf_ind <- gather(cor_matrices, algorithm, correlation) %>% 
    mutate(dataset = names(sce_list_with_marker_pseudotime)[i])
  
  sdf <- rbind(sdf, sdf_ind)
  
  return(sdf)
})

a_tidy <- bind_rows(amat_list)

# a_tidy$algorithm <- factor(a_tidy$algorithm, levels = c("ouija", "pc1",
#                                                         "monocle", "tscan", "dpt"))
a_tidy$algorithm <- factor(a_tidy$algorithm, levels = c("ouija", "pc1",
                                                        "tscan", "dpt"))

ggplot(data = a_tidy, aes(x = algorithm, y = correlation, fill = algorithm)) +
  geom_boxplot() + facet_wrap(~ dataset) +
  xlab("Algorithm") +
  ylab("Correlation across marker genes")

across_algs <- last_plot()

cowplot::plot_grid(marker_plot, global_plot, within_plot, ncol = 1, labels = "AUTO")

full_plot <- last_plot()

ggsave("figs/marker-vs-txome-pseudotime-subsets.png", width=12,height=12)


stop("Done")



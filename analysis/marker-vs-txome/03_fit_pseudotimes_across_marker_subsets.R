set.seed(123L)

source("scripts/pseudotime_fitting_functions.R")

sce_list_with_pseudotime <- readRDS("data/scesets/scesets_with_pseudotime.rds")
markers_to_use <- readRDS("data/marker-vs-txome/markers_to_use.rds")

n_additional_markers <- c(1,5,10,20,50,100,500,1000)

n_reps <- 50


# Fit marker pseudotimes first --------------------------------------------

sce_list_with_marker_pseudotime <- lapply(sce_list_with_pseudotime, function(sce) {
  sce_marker <- sce[fData(sce)$is_marker, ]
  
  sce$pc1_marker_pseudotime <- fit_pc1_pseudotime(sce_marker)
  sce$tscan_marker_pseudotime <- fit_tscan_pseudotime(sce_marker)
  sce$monocle_marker_pseudotime <- fit_monocle_pseudotime(sce_marker)
  sce$dpt_marker_pseudotime <- fit_dpt_pseudotime(sce_marker)
  
  return(sce)
})

saveRDS(sce_list_with_marker_pseudotime, 
        file = "data/scesets/scesets_with_marker_pseudotime.rds")



# Consistency calculations ------------------------------------------------


## Consistency matrix

cmats <- lapply(seq_along(sce_list_with_marker_pseudotime), function(i) {
  sce <- sce_list_with_marker_pseudotime[[i]]
  
  # We have a 
  # replication x number_of_markers x algorithm x cell array
  # to store the pseudotimes
  cmat <- array(dim = c(n_reps, length(n_additional_markers),
                        5, ncol(sce)))
  
  for(rep in 1:n_reps) {
    for(nam_index in seq_along(n_additional_markers)) {
      genes_to_use <- c(which(fData(sce)$is_marker),
                        markers_to_use[[i]][[rep]][[nam_index]])
      sce_reduced <- sce[genes_to_use,]
      
      cmat[rep, nam_index, 1, ] <- fit_monocle_pseudotime(sce_reduced)
      cmat[rep, nam_index, 2, ] <- fit_tscan_pseudotime(sce_reduced)
      
      ## Need to set the ouija priors
      is_marker_gene <- fData(sce)$is_marker
      k_mean <- fData(sce)$k_mean[is_marker_gene]
      k_mean <- c(k_mean, rep(0, n_additional_markers[nam_index]))
      
      # cmat[rep, nam_index, 3, ] <- fit_ouija_pseudotime(sce_reduced, inference_type = "vb",
                                                        # k_mean = k_mean, marker_only = FALSE)
      
      cmat[rep, nam_index, 4, ] <- fit_pc1_pseudotime(sce_reduced)
      cmat[rep, nam_index, 5, ] <- fit_dpt_pseudotime(sce_reduced)
    }
  }
  
  return(cmat)
  
})

saveRDS(cmats, file = "data/marker-vs-txome/pseudotimes_across_marker_subset_array.rds")









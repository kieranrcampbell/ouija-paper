set.seed(123L)

library(scater)

sce_list_with_pseudotime <- readRDS("data/scesets/scesets_with_pseudotime.rds")

n_additional_markers <- c(1,5,10,20,50,100,500,1000)

n_reps <- 50


## First choose the additional genes to use
## This list is formatted as
## markers_to_use$dataset[[replication]][[nmarkers]]

markers_to_use <- lapply(sce_list_with_pseudotime, function(sce) {
  non_marker_inds <- which(!fData(sce)$is_marker)
  replicate(n_reps, 
    lapply(n_additional_markers, function(n) {
      sample(non_marker_inds, n, replace = FALSE)
    }), simplify = FALSE
  )
})

save(markers_to_use, file = "data/marker-vs-txome/markers_to_use.rds")

suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(dpt))
# Pseudotime libraries

fit_monocle_pseudotime <- function(sce, ncenter = NULL) {
  fData(sce)$gene_short_name <- featureNames(sce)
  cds <- toCellDataSet(sce)
  cds <- monocle::setOrderingFilter(cds, ordering_genes = featureNames(sce))
  sizeFactors(cds) <- rep(1, ncol(cds))
  cds <- tryCatch(monocle::reduceDimension(cds, norm_method = "none", pseudo_expr = 0, ncenter = ncenter),
                  error = function(e) {
                    message("Monocle dimensionality reduction failed: ", e)
                    return(NA)
                  })
  if(!is(cds, "CellDataSet")) {
    if(is.na(cds)) {
      return(rep(NA, ncol(sce)))
    }
  }
  
  cds <- tryCatch(monocle::orderCells(cds, num_paths = 1),
                  error = function(e) {
                    message("Monocle ordering failed")
                    return(NA)
                  })
  if(!is(cds, "CellDataSet")) {
    if(is.na(cds)) {
      return(rep(NA, ncol(sce)))
    }
  }
  unloadNamespace("monocle")
  R.utils::gcDLLs()
  return(cds$Pseudotime)
}

fit_tscan_pseudotime <- function(sce, clusternum = 5) {
  tscan_pseudotime <- rep(NA, ncol(sce))
  
  cl_data <- tryCatch(TSCAN::exprmclust(exprs(sce), clusternum = clusternum),
                      error = function(e) {
                        message("TSCAN failed (clustering) ", e)
                        return(NA)
                      })
  
  if(!is.list(cl_data)) {
    if(is.na(cl_data)) {
      return(tscan_pseudotime)
    }
  }
  
  tscan_order <- tryCatch(TSCAN::TSCANorder(cl_data, orderonly = FALSE),
                          error = function(e) {
                            message("TSCAN failed (ordering) ", e)
                            return(rep(NA, ncol(sce)))
                          })
  

  if(any(is.na(tscan_order))) return( tscan_pseudotime )
  
  tscan_cell_inds <- match(tscan_order$sample_name, sampleNames(sce))
  tscan_pseudotime[tscan_cell_inds] <- tscan_order$Pseudotime
  unloadNamespace("TSCAN")
  R.utils::gcDLLs()
  return( tscan_pseudotime )
}

fit_pc1_pseudotime <- function(sce) {
  pc1 <- prcomp(t(exprs(sce)))$x[,1]
  return(pc1)
}


fit_ouija_pseudotime <- function(sce, marker_only = TRUE, k_mean = NULL,
                                 iter = 10000, inference_type = "hmc") {
  ouija_pseudotime <- NULL
  
  if(marker_only) {
    is_marker_gene <- fData(sce)$is_marker
    k_mean <- fData(sce)$k_mean[is_marker_gene]
    k_sd <- fData(sce)$k_sd[is_marker_gene]
    if(any(is.na(k_sd))) k_sd <- NULL
    oui <- tryCatch(ouija(sce[is_marker_gene, ], strengths = k_mean, strength_sd = k_sd,
                          iter = iter, inference_type = inference_type),
                    error = function(e) {
                      message(paste("Ouija error!"))
                      return(NA)
                    })
    if(is.na(oui)) {
      return(rep(NA, ncol(sce)))
    }
    ouija_pseudotime <- map_pseudotime(oui)
  } else {
    oui <- tryCatch(ouija(sce, strengths = k_mean,
                          iter = iter, inference_type = inference_type),
                    error = function(e) {
                      message(paste("Ouija error!"))
                      return(NA)
                    })
    if(is.na(oui)) {
      return(rep(NA, ncol(sce)))
    }    
    ouija_pseudotime <- map_pseudotime(oui)
  }
  return(ouija_pseudotime)
}


fit_dpt_pseudotime <- function(sce) {
  # we need to add a small amount of noise to DPT in order
  # to avoid the 'non-unique' error
  gene_sds <- matrixStats::rowSds(exprs(sce))
  gene_collapse_sd <- 1e-6 * gene_sds
  new_exprs <- apply(exprs(sce), 2, function(x) rnorm(length(x), x, gene_collapse_sd))
  
  ts <- dpt::Transitions(t(new_exprs))
  pt <- dpt::dpt(ts, branching = FALSE)  
  unloadNamespace("dpt")
  R.utils::gcDLLs()
  return(pt$DPT)
}





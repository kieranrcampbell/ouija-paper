
# Pseudotime libraries
library(dpt)
library(monocle)
library(TSCAN)
library(ouija)


fit_monocle_pseudotime <- function(sce) {
  cds <- toCellDataSet(sce)
  cds <- setOrderingFilter(cds, ordering_genes = featureNames(sce))
  cds <- tryCatch(reduceDimension(cds, norm_method = "none"),
                  error = function(e) {
                    message("Monocle dimensionality reduction failed")
                    return(NA)
                  })
  if(!is(cds, "CellDataSet")) {
    if(is.na(cds)) {
      return(rep(NA, ncol(sce)))
    }
  }
  
  cds <- tryCatch(orderCells(cds, num_paths = 1),
                  error = function(e) {
                    message("Monocle ordering failed")
                    return(NA)
                  })
  if(!is(cds, "CellDataSet")) {
    if(is.na(cds)) {
      return(rep(NA, ncol(sce)))
    }
  }
  return(cds$Pseudotime)
}

fit_tscan_pseudotime <- function(sce) {
  tscan_pseudotime <- rep(NA, ncol(sce))
  
  cl_data <- tryCatch(exprmclust(exprs(sce)),
                      error = function(e) {
                        message("TSCAN failed")
                        return(NA)
                      })
  
  if(!is.list(cl_data)) {
    if(is.na(cl_data)) {
      return(tscan_pseudotime)
    }
  }
  
  tscan_order <- tryCatch(TSCANorder(cl_data, orderonly = FALSE),
                          error = function(e) {
                            message("TSCAN failed")
                            return(rep(NA, ncol(sce)))
                          })
  

  if(any(is.na(tscan_order))) return( tscan_pseudotime )
  
  tscan_cell_inds <- match(tscan_order$sample_name, sampleNames(sce))
  tscan_pseudotime[tscan_cell_inds] <- tscan_order$Pseudotime
  return( tscan_pseudotime )
}

fit_pc1_pseudotime <- function(sce) {
  pc1 <- prcomp(t(exprs(sce)))$x[,1]
  return(pc1)
}







library(readxl)
library(stringr)
library(scater)
library(tidyverse)
library(mclust)

set.seed(123L)


add_marker_info <- function(sce, marker_inds, k_m, k_sds = NULL) {
  ngenes <- nrow(sce)
  is_marker <- rep(FALSE, ngenes)
  is_marker[marker_inds] <- TRUE
  k_mean <- k_sd <- rep(NA, ngenes)
  k_mean[marker_inds] <- k_m  
  if(!is.null(k_sds)) {
    k_sd[marker_inds] <- k_sds
  }
  fData(sce)$is_marker <- is_marker
  fData(sce)$k_mean <- k_mean
  fData(sce)$k_sd <- k_sd
  
  return(sce)
}

raw <- read_excel("data/raw/mmc5.xlsx")

feature_names <- raw$X__1
raw$X__1 <- NULL
exprs_mat <- as.matrix(raw)
rownames(exprs_mat) <- feature_names

sample_names <- names(raw) 

cell_type <- sapply(str_split(sample_names, "_"), `[`, 1)

pdata <- data.frame(cell_type = cell_type)
rownames(pdata) <- sample_names

sce <- newSCESet(exprsData = log2(exprs_mat + 1), 
                 phenoData = AnnotatedDataFrame(pdata))
sce <- calculateQCMetrics(sce)

sce <- sce[, sce$cell_type %in% c("qNSC", "aNSC", "NPC")]

# Cluster away ologodendrocytes
sce <- plotPCA(sce, colour_by = "cell_type", ncomponents =3, return_SCESet = TRUE)

mc <- Mclust(redDim(sce)[,1:2], G = 2)
sce$cluster <- mc$classification

# remove the cluster with fewest cells
to_keep <- which.max(tabulate(sce$cluster))
sce_pst <- sce[, sce$cluster == to_keep]

sce_pst <- plotPCA(sce_pst, colour_by = "cell_type", return_SCESet = TRUE)
sce_pst$PC1 <- redDim(sce_pst)[,1]


# "Id3", 
markers <- c("Id3", "Clu", "Rpl32", "Egfr",
             "Cdk4", "Cdk1", "Dlx2", "Dcx")

fData(sce_pst)$is_marker <- featureNames(sce_pst) %in% markers

saveRDS(sce_pst, "data/scesets/dulken-sce.rds")


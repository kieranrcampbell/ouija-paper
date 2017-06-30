library(scater)
library(tidyverse)
library(monocle)
library(HSMMSingleCell)
library(readxl)

#' Here we use genes that have a coefficient of variation 
#' greater than 0.5 and an 'exprs' greater than 1 in any cell:

filter_sceset <- function(sce) {
  variance <- matrixStats::rowVars(exprs(sce))
  any_1 <- apply(exprs(sce), 1, function(x) any(x > 1))
  
  to_use <- (variance > 1 & !is.na(variance) & any_1) | fData(sce)$is_marker
  return(sce[to_use, ])
}

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

# Trapnell ----------------------------------------------------------------

data(HSMM)
HSMM_trajectory <- HSMM[,HSMM$State != 3]

sce_trapnell <- fromCellDataSet(HSMM_trajectory, "fpkm")

## add marker info
trapnell_marker_genes <- c("CDK1", "ID1", "MYOG", "MEF2C", "MYH3")
trapnell_marker_gene_inds <- match(trapnell_marker_genes, 
                                   fData(sce_trapnell)$gene_short_name)
sign_bits <- c(-1, -1, 1, 1, 1)
k_means <- sign_bits * 10

sce_trapnell <- add_marker_info(sce_trapnell, trapnell_marker_gene_inds,
                                k_means)

sce_trapnell <- filter_sceset(sce_trapnell)


# HSCs --------------------------------------------------------------------

## this requires "scripts/read-hsc.R" to be run first

load("../data/hsc-sce.Rdata")
sce <- updateSCESet(sce)
sce_hsc <- sce[, which(!(sce$cell_type %in% c("Adult_HSC", "E11.0_T1CD201neg", "E11.0_T2_cd41_low")))]


ct <- plyr::mapvalues(as.character(sce_hsc$cell_type),
                      from = c("E11.0_EC", "E11.0_T1", "E11.0_T2_pre_hsc",
                               "E12.5_FL", "E14.5_FL"),
                      to = c("EC", "T1 pre-HSC", "T2 pre-HSC",
                             "E12 HSC", "E14 HSC"))
sce_hsc$CellType <- ct

## add marker info

hsc_marker_genes <- c("Nrp1", "Hey1", "Efnb2", "Ephb4", "Nrp2", "Nr2f2")
hsc_marker_inds <- match(hsc_marker_genes, featureNames(sce_hsc))
hsc_k_means <- 10 * rep(-1, 6)

sce_hsc <- add_marker_info(sce_hsc, hsc_marker_inds, hsc_k_means)


# Dataset specific outlier removal
# set.seed(123L)
# sce_hsc <- plotPCA(sce_hsc, colour_by = "CellType", 
#               ntop = nrow(sce_hsc), scale_features = T, return_SCESet = TRUE)

# to_keep <- redDim(sce_hsc)[,2] > -50 
# sce_hsc <- sce_hsc[, to_keep]

sce_hsc <- filter_sceset(sce_hsc)



# Shin --------------------------------------------------------------------

filename <- "data/raw/waterfall_data.xlsx"
if(!file.exists(filename)) {
  download.file("http://www.cell.com/cms/attachment/2038326541/2052521610/mmc7.xlsx",
                filename)
}

d <- read_excel(filename)

d <- as.data.frame(d)
cellnames <- as.character(d[1,-1])
genenames <- d[-c(1,2, nrow(d)),1]
w_pseudotimes <- as.numeric(d[2,-1])
tpm_data <- as.matrix(d[-c(1,2, nrow(d)),-1])
tpm_data <- t(apply(tpm_data, 1, as.numeric))

colnames(tpm_data) <- cellnames
rownames(tpm_data) <- genenames

pd <- new("AnnotatedDataFrame", data.frame(wpst = w_pseudotimes))
rownames(pd) <- cellnames

sce_shin <- newSCESet(tpmData = tpm_data, phenoData = pd, logExprsOffset = 1)
is_exprs(sce_shin) <- exprs(sce_shin) > 0
sce_shin <- calculateQCMetrics(sce_shin)

## add marker info
shin_marker_genes <- c("Sox11", "Eomes", "Stmn1", "Apoe", "Aldoc","Gfap")
shin_marker_inds <- match(shin_marker_genes, featureNames(sce_shin))

sign_bits <-c(1, 1, 1, -1, -1, -1)
k_means <- sign_bits * 5

sce_shin <- add_marker_info(sce_shin, shin_marker_inds, k_means)

sce_shin <- filter_sceset(sce_shin)



# Save all SCESets --------------------------------------------------------

sce_list <- list(
  trapnell = sce_trapnell,
  hsc = sce_hsc,
  shin = sce_shin
)


saveRDS(sce_list, file = "data/scesets/scesets.rds")




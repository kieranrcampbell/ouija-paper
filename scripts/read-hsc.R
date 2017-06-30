library(readr)
library(scater)

hsc_dir <- "../data/GSE67120_RAW/"
files <- dir(hsc_dir)

parse_name <- function(f) {
  token <- strsplit(f, ".genes", fixed = TRUE)[[1]][1]
  by_ <- strsplit(token, "_", fixed = TRUE)[[1]]
  name <- by_[1]
  cell_type <- paste(by_[2:3], collapse = "_")
  return(c(name, cell_type))
}

pd <- data.frame(t(sapply(files, parse_name)))
names(pd) <- c("name", "cell_type")
pd$cell_type <- as.character(pd$cell_type)
rownames(pd) <- pd$name

## Time for a bit of voodoo - according to cell numbers, the final 22 cells are T2 pre-HSC
## and the ones before that are T2 CD41-low

t2_inds <- grep("E11.0_T2", pd$cell_type, fixed = TRUE)
t2_cd41_low_inds <- t2_inds[1:26]
t2_prehsc_inds <- t2_inds[27:48]

pd$cell_type[t2_cd41_low_inds] <- paste0(pd$cell_type[t2_cd41_low_inds], "_cd41_low")
pd$cell_type[t2_prehsc_inds] <- paste0(pd$cell_type[t2_prehsc_inds], "_pre_hsc")

# probe first file to get layout
gex_table <- read_tsv(file.path(hsc_dir, files[1]))
gene_id <- gex_table$gene_id

## need a first pass to find genes common to all cells
for(i in seq_along(files)) {
  gti <- read_tsv(file.path(hsc_dir, files[i]))
  gene_id <- intersect(gene_id, gti$gene_id)
}

## now we can build our fpkm matrix
n_genes <- length(gene_id)
fpkm_data <- matrix(NA, ncol = nrow(pd), nrow = n_genes)
rownames(fpkm_data) <- gene_id
colnames(fpkm_data) <- pd$name

for(i in seq_along(files)) {
  gti <- read_tsv(file.path(hsc_dir, files[i]))
  where_gene <- match(gene_id, gti$gene_id)
  gti <- gti[where_gene, ]
  stopifnot(all.equal(gti$gene_id, gene_id))
  fpkm_data[,i] <- gti$FPKM
}

sce <- newSCESet(fpkmData = fpkm_data,
                 phenoData = new("AnnotatedDataFrame", pd))

sce <- calculateQCMetrics(sce, feature_controls = grep("ERCC", featureNames(sce)))

cell_stage <- as.character(sce$cell_type)
cell_stage[grep("E11", cell_stage)] <- "E11"

saveRDS(sce, file = "data/scesets/hsc-sce.rds")

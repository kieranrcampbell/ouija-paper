
library(scater)
library(tidyverse)

expr_df <- read_tsv("data/raw/GSE76157_Expression.txt")
expr_df <- rename(expr_df, gene = X1)
expr_df <- head(expr_df, -1) # last row is just NAs

expr_mat <- as.matrix(select(expr_df, -gene))
rownames(expr_mat) <- expr_df$gene


sce <- newSCESet(exprsData = expr_mat)
sce <- calculateQCMetrics(sce)

sce <- plotPCA(sce, ntop = 5e3, return_SCESet = TRUE)


# Remove “outlier” cells --------------------------------------------------

sce$is_outlier <- exprs(sce)["Col9a2",] < 9 | sce$pct_dropout > 80


sce <- sce[, !sce$is_outlier]

saveRDS(sce, "data/scesets/li-sce.rds")



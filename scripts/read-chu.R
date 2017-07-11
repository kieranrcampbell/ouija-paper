library(scater)
library(tidyverse)
library(stringr)
devtools::load_all("~/oxford/bnlfa/oui/ouija")

expr_df <- read_csv("data/raw/GSE75748_sc_time_course_ec.csv.gz")
expr_df <- rename(expr_df, gene = X1)

expr_mat <- as.matrix(select(expr_df, -gene))
rownames(expr_mat) <- expr_df$gene

sample_names <- colnames(expr_mat)

ss <- strsplit(sample_names, ".", fixed = T)
split2 <- strsplit(sapply(ss, `[`, 2), "h")
capture_time <- sapply(split2, `[`, 1)
split3 <- sapply(split2, `[`, 2)


cell_names <- paste0("cell_", seq_len(ncol(expr_mat)))
pdata <- data.frame(capture_time)
rownames(pdata) <- colnames(expr_mat) <- cell_names

sce <- newSCESet(countData = expr_mat,
                 phenoData = AnnotatedDataFrame(pdata))
sce <- calculateQCMetrics(sce)

sce <- normaliseExprs(sce, method = "TMM")

plotPCA(sce, colour_by = "capture_time", ntop = 1e4)

markers <- c("POU5F1", "NANOG", "SOX2", "EOMES", "CER1", "GATA4", "DKK4",
             "MYCT1", "PRDM1")
fData(sce)$is_marker <- featureNames(sce) %in% markers

saveRDS(sce, "data/scesets/chu-sce.rds")

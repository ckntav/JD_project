# setwd("/Users/chris/Desktop/JD_project")

library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)

# Load raw counts
raw_counts <- read_tsv("output/rna-pipeline-hg38/DGE/rawCountMatrix.csv")
colnames(raw_counts)

##### Get gene lengths
geneLengths <- read_tsv("input/ensembl/Homo_sapiens.GRCh38.Ensembl87.genes.length.tsv",
                        col_names = c("Gene", "geneLength"))

# make TPM matrix
df <- left_join(raw_counts, geneLengths, by = "Gene")
is.na(df$geneLength) %>% sum

raw_counts_mat <- raw_counts %>% dplyr::select(-Gene, -Symbol) %>% as.matrix
rownames(raw_counts_mat) <- raw_counts$Gene
geneLengths_vector <- df$geneLength
x <- raw_counts_mat/geneLengths_vector
TPM_raw <- t(t(x)*1e6/colSums(x))
TPM_df <- data.frame(Gene = raw_counts$Gene, Symbol = raw_counts$Symbol, TPM_raw)
colnames(TPM_df) <- colnames(raw_counts)

write_tsv(TPM_df, path = "output/rna-pipeline-hg38/DGE/TPM_matrix.tsv")

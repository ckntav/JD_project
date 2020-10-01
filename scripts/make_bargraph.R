# setwd("/Users/chris/Desktop/JD_project")

library(tidyverse)
library(knitr)
library(ComplexHeatmap)
library(viridis)

#
samples <- c("UMSCC1_NTsh", "UMSCC47_NTsh")
gene_symbols <- c("FUCA2", "NFYA", "BRD4", "NANOG", "NBPF8")

#
TPM_df <- read_tsv(file = "output/rna-pipeline-hg38/DGE/TPM_matrix.tsv") %>%
  dplyr::select(-Gene) %>% 
  dplyr::filter(Symbol %in% gene_symbols)

kable(TPM_df)

#
mat <- TPM_df %>% dplyr::select(-Symbol)
rownames(mat) <- TPM_df$Symbol

data <- data.frame(samples, t(mat)) %>% 
  gather(key = "Symbol", value = "TPM", -samples)
data$samples <- factor(data$samples, levels = samples)
data$Symbol <- factor(data$Symbol, levels = gene_symbols)

barplot <- 
  ggplot(data, aes(x = Symbol, y = TPM, fill = samples)) +
  geom_bar(stat = "identity", position = "dodge", size = 0.5, col = "black") +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Samples", fill = "Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"))

barplot
load("/home/dermot.kelly/Dermot_analysis/Other/data/DE_results_BRCA_miRNA.RData")

library(DESeq2)

vsd <- varianceStabilizingTransformation(dds_er, blind = F) # Blind = F, respects design and accounts for covariates
vst_matrix <- assay(vsd)

# Remove genes with zero variance
vst_matrix_filtered <- vst_matrix[apply(vst_matrix, 1, var) >0, ]

nrow(vst_matrix)
nrow(vst_matrix_filtered)

pca <- prcomp(t(vst_matrix_filtered), scale. = TRUE)  # transpose so samples are rows
pca_df <- as.data.frame(pca$x)
pca_df$ER_status_created <- colData(vsd)$ER_status_created


library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = ER_status_created)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA of VST-transformed miRNA expression",
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% variance)")
  )

View(pca_df)

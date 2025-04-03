# DeSeq2 Pipeline:

setwd("/home/ibab/HM_VM_Project/DGEA/star_htseq_counts")

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(biomaRt)
x = read.delim("htseq_counts_matrix.txt", row.names = "Gene_Name")
x$row.names = NULL
nrow(x)
head(x)

# Convert ensmbl ids to gene symbols:
# Connect to Ensembl BioMart
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Extract rownames (Ensembl IDs) from your dataframe
ensembl_ids = rownames(x)  # Replace 'my_df' with your dataframe name

# Convert Ensembl IDs to Gene Symbols
genes = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
               filters = "ensembl_gene_id", 
               values = ensembl_ids, 
               mart = ensembl)

# Merge with original dataframe
x$ensembl_gene_id = rownames(x)  # Add Ensembl IDs as a column
x = merge(x, genes, by = "ensembl_gene_id", all.x = TRUE)
x = na.omit(x)
x = x[!(is.na(x$hgnc_symbol) | x$hgnc_symbol == ""), ]


x2 = x[!duplicated(x$hgnc_symbol),]
rownames(x2) = x2$hgnc_symbol
x2$ensembl_gene_id = NULL
x2$hgnc_symbol = NULL

nrow(x2)
write.csv(x = x2, file = "final_counts_matrix_all.csv", row.names = T)

# # AT vs Rest
# at_vs_rest_metadata = read.csv("at_vs_rest_metadata.csv", row.names = 1)
# stopifnot(all(colnames(x2) == rownames(at_vs_rest_metadata)))
# 
# dds = DESeqDataSetFromMatrix(countData = x2,
#                               colData = at_vs_rest_metadata,
#                               design = ~ Condition)
# 
# 
# dds = dds[rowSums(counts(dds)) > 10, ]
# 
# # # Boxplot Before Normalization:
# # 
# # boxplot(log2(x2 + 1), 
# #         main = "Raw Counts Before Normalization", 
# #         las = 2, 
# #         col = rainbow(ncol(x2)), 
# #         outline = FALSE)
# 
# 
# vsd_raw <- vst(dds, blind = TRUE)  # Variance stabilizing transformation
# plotPCA(vsd_raw, intgroup = "Condition") + ggtitle("PCA Before Normalization")
# 
# 
# # Perform DESeq:
# dds = DESeq(dds)
# 
# # # Boxplot after normalization
# # normalized_counts = counts(dds, normalized = TRUE)
# # 
# # boxplot(log2(normalized_counts + 1), 
# #         main = "After Normalization", 
# #         las = 2, 
# #         col = rainbow(ncol(normalized_counts)), 
# #         outline = FALSE)
# 
# vsd <- vst(dds, blind = FALSE)
# plotPCA(vsd, intgroup = "Condition") + ggtitle("PCA After Normalization")
# 
# 
# res = results(dds, contrast = c("Condition", "AT", "Control"))
# res = res[order(res$padj, na.last = NA), ]
# res_df = as.data.frame(res)
# at_vs_rest_upgenes = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
# at_vs_rest_downgenes = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)
# 
# write.csv(at_vs_rest_upgenes, file = "at_vs_rest_upgenes.csv", row.names = T)
# write.csv(at_vs_rest_downgenes, file = "at_vs_rest_downgenes.csv", row.names = T)
# 
# 
# # write.csv(as.data.frame(res), file = "AT_vs_Rest_DEG_results.csv", row.names = TRUE)
# 
# # CE vs Rest
# 
# ce_vs_rest_metadata = read.csv("ce_vs_rest_metadata.csv", row.names = 1)
# stopifnot(all(colnames(x2) == rownames(ce_vs_rest_metadata)))
# 
# dds = DESeqDataSetFromMatrix(countData = x2,
#                              colData = ce_vs_rest_metadata,
#                              design = ~ Condition)
# 
# 
# dds = dds[rowSums(counts(dds)) > 10, ]
# 
# # Perform DESeq:
# dds = DESeq(dds)
# 
# res = results(dds, contrast = c("Condition", "Case", "Control"))
# res = res[order(res$padj, na.last = NA), ]
# 
# res_df = as.data.frame(res)
# ce_vs_rest_upgenes = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
# ce_vs_rest_downgenes = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)
# 
# write.csv(ce_vs_rest_upgenes, file = "ce_vs_rest_upgenes.csv", row.names = T)
# write.csv(ce_vs_rest_downgenes, file = "ce_vs_rest_downgenes.csv", row.names = T)
# # write.csv(as.data.frame(res), file = "CE_vs_Rest_DEG_results.csv", row.names = TRUE)
# 
# 
# # ESUS vs Rest
# 
# esus_vs_rest_metadata = read.csv("esus_vs_rest_metadata.csv", row.names = 1)
# stopifnot(all(colnames(x2) == rownames(esus_vs_rest_metadata)))
# 
# dds = DESeqDataSetFromMatrix(countData = x2,
#                              colData = esus_vs_rest_metadata,
#                              design = ~ Condition)
# 
# 
# dds = dds[rowSums(counts(dds)) > 10, ]
# 
# # Perform DESeq:
# dds = DESeq(dds)
# 
# res = results(dds, contrast = c("Condition", "Case", "Control"))
# res = res[order(res$padj, na.last = NA), ]
# 
# res = results(dds, contrast = c("Condition", "Case", "Control"))
# res = res[order(res$padj, na.last = NA), ]
# 
# res_df = as.data.frame(res)
# esus_vs_rest_upgenes = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
# esus_vs_rest_downgenes = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)
# 
# write.csv(esus_vs_rest_upgenes, file = "esus_vs_rest_upgenes.csv", row.names = T)
# write.csv(esus_vs_rest_downgenes, file = "esus_vs_rest_downgenes.csv", row.names = T)
# # write.csv(as.data.frame(res), file = "ESUS_vs_Rest_DEG_results.csv", row.names = TRUE)

x2 = read.csv(file = "final_counts_matrix_all.csv", row.names = 1)
# Pairwise: 
# AT vs CE


at_ce = x2[, 0:16]
at_ce <- at_ce[, !colnames(at_ce) %in% c("ce2", "ce4", "at6", "ce6", "ce5", "ce7", "at5", "at1")]
write.csv(at_ce, file = "at_vs_ce_final_counts_matrix.csv", sep = ",", row.names = T)

# at_ce = read.csv(file = "at_vs_ce_final_counts_matrix.csv")
# rownames(at_ce) = at_ce$Gene
# at_ce$Gene = NULL
at_vs_ce_metadata = read.csv("at_vs_ce_rem2.csv", row.names = 1)
dds = DESeqDataSetFromMatrix(countData = at_ce,
                             colData = at_vs_ce_metadata,
                             design = ~ Condition)
dds = dds[rowSums(counts(dds)) > 10, ]
nrow(dds)
vsd <- vst(dds, blind = TRUE)  # Variance stabilizing transformation

pca_data = plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# Extract percentage variance explained
percentVar = round(100 * attr(pca_data, "percentVar"))

# Add sample names
# Create PCA plot with sample names and variance info
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Removing Outliers") +
  theme_minimal()

# vsd_raw <- vst(dds, blind = TRUE)  # Variance stabilizing transformation
# plotPCA(vsd_raw, intgroup = "Condition") + ggtitle("PCA Before Normalization")
# dds$Condition <- relevel(dds$Condition, ref = "CE")  # Set CE as reference
dds = DESeq(dds)  
resultsNames(dds)  # Now should list "Condition_AT_vs_CE"
# res_shrink = lfcShrink(dds, coef="Condition_CE_vs_AT", type="apeglm")

# # dds = DESeq(dds)
# vsd = vst(dds, blind = F)
# # Compute PCA data
# pca_data = plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# # Extract percentage variance explained
# percentVar = round(100 * attr(pca_data, "percentVar"))
# 
# # Add sample names
# # Create PCA plot with sample names and variance info
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_data))) +
#   geom_point(size = 3) +
#   geom_text_repel() +  # Add sample names
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label PC1 with variance
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label PC2 with variance
#   ggtitle("PCA AT vs CE") +
#   theme_minimal()
# resultsNames(dds)

# Perform LFC shrinkage
# res_shrink = lfcShrink(dds, coef="Condition_CE_vs_AT", type="apeglm")

# # Convert to data frame
# res_shrink_df = as.data.frame(res_shrink)
# 
# # Subset significant genes with |log2FC| ≥ 1.5 and padj ≤ 0.01
# at_vs_ce_upgenes_1_5_shrink = subset(res_shrink_df, log2FoldChange >= 1.5 & padj <= 0.05)
# at_vs_ce_downgenes_1_5_shrink = subset(res_shrink_df, log2FoldChange <= -1.5 & padj <= 0.05)
# 
# at_vs_ce_upgenes_1_shrink = subset(res_shrink_df, log2FoldChange >= 1 & padj <= 0.05)
# at_vs_ce_downgenes_1_shrink = subset(res_shrink_df, log2FoldChange <= -1 & padj <= 0.05)
# 
# at_vs_ce_upgenes_0_5_shrink = subset(res_shrink_df, log2FoldChange >= 0.5 & padj <= 0.05)
# at_vs_ce_downgenes_0_5_shrink = subset(res_shrink_df, log2FoldChange <= -0.5 & padj <= 0.05)
# 
# 
# write.csv(x = at_vs_ce_upgenes_1_5_shrink, file = "at_vs_ce_upgenes_1_5_shrink.csv")
# write.csv(x = at_vs_ce_downgenes_1_5_shrink, file = "at_vs_ce_downgenes_1_5_shrink.csv")
# write.csv(x = at_vs_ce_upgenes_1_shrink, file = "at_vs_ce_upgenes_1_shrink.csv")
# write.csv(x = at_vs_ce_downgenes_1_shrink, file = "at_vs_ce_downgenes_1_shrink.csv")
# write.csv(x = at_vs_ce_upgenes_0_5_shrink, file = "at_vs_ce_upgenes_0_5_shrink.csv")
# write.csv(x = at_vs_ce_downgenes_0_5_shrink, file = "at_vs_ce_downgenes_0_5_shrink.csv")



# vsd = vst(dds, blind = FALSE, fitType = "mean")
# plotPCA(vsd, intgroup = "rownames(at_vs_ce_metadata)") + ggtitle("PCA After Normalization")

# # Extract PCA data with variance explained
# pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# pca_data$Sample <- rownames(at_vs_ce_metadata)  # Assign sample names from metadata
# 
# # Perform PCA manually to get variance explained
# pca <- prcomp(t(assay(vsd)))
# percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# # Add percent variance to axis labels
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
#   geom_point(size = 3) +
#   geom_label(aes(label = Sample), vjust = 1.5, size = 3) +
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#   ggtitle("PCA After Normalization: AT vs CE") +
#   theme_minimal()

res = results(dds, contrast = c("Condition", "AT", "CE"))
res = res[order(res$padj, na.last = NA), ]
res_df = as.data.frame(res)
res_p_val_cut = subset(res_df, res_df$padj <= 0.05)
at_vs_ce_upgenes_1_5 = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
at_vs_ce_downgenes_1_5 = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)

at_vs_ce_upgenes_1 = subset(res_df, res_df$log2FoldChange >= 1 & res_df$padj <= 0.01)
at_vs_ce_downgenes_1 = subset(res_df, res_df$log2FoldChange <= -1 & res_df$padj <= 0.01)

at_vs_ce_upgenes_0_5 = subset(res_df, res_df$log2FoldChange >= 0.5 & res_df$padj <= 0.01)
at_vs_ce_downgenes_0_5 = subset(res_df, res_df$log2FoldChange <= -0.5 & res_df$padj <= 0.01)

write.csv(at_vs_ce_upgenes_1_5, file = "at_vs_rest_upgenes_1_5.csv", row.names = T)
write.csv(at_vs_ce_downgenes_1_5, file = "at_vs_rest_downgenes_1_5.csv", row.names = T)

write.csv(at_vs_ce_upgenes_1, file = "at_vs_rest_upgenes_1.csv", row.names = T)
write.csv(at_vs_ce_downgenes_1, file = "at_vs_rest_downgenes_1.csv", row.names = T)

write.csv(at_vs_ce_upgenes_0_5, file = "at_vs_rest_upgenes_0_5.csv", row.names = T)
write.csv(at_vs_ce_downgenes_0_5, file = "at_vs_rest_downgenes_0_5.csv", row.names = T)

# ESUS vs CE
esus_ce = x2[, 7:21]
esus_ce <- esus_ce[, !colnames(esus_ce) %in% c("esus1", "ce10", "ce6", "ce4",
                                               "esus3", "ce1", "ce9", "ce8",
                                               "ce3")]
write.csv(esus_ce, file = "esus_ce_final_counts_matrix.csv")

# esus_ce = read.csv(file = "esus_ce_final_counts_matrix.csv")
# rownames(esus_ce) = esus_ce$X
# esus_ce$X = NULL

esus_vs_ce_metadata = read.csv("esus_vs_ce_rem2.csv", row.names = 1)
dds = DESeqDataSetFromMatrix(countData = esus_ce,
                             colData = esus_vs_ce_metadata,
                             design = ~ Condition)

dds = dds[rowSums(counts(dds)) > 10, ]
nrow(dds)
# vsd_raw <- vst(dds, blind = TRUE)  # Variance stabilizing transformation
# plotPCA(vsd_raw, intgroup = "Condition") + ggtitle("PCA Before Normalization")


vsd <- vst(dds, blind = TRUE)  # Variance stabilizing transformation

pca_data = plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# Extract percentage variance explained
percentVar = round(100 * attr(pca_data, "percentVar"))

# Add sample names
# Create PCA plot with sample names and variance info
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Removing Outliers") +
  theme_minimal()



dds = DESeq(dds)
vsd = vst(dds, blind = F)
pca_data = plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# Extract percentage variance explained
percentVar = round(100 * attr(pca_data, "percentVar"))

# Add sample names
# Create PCA plot with sample names and variance info
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA CE vs ESUS") +
  theme_minimal()

resultsNames(dds)

# Perform LFC shrinkage
res_shrink = lfcShrink(dds, contrast = c("Condition", "CE", "ESUS"), type="ashr")

# Convert to data frame
res_shrink_df = as.data.frame(res_shrink)

# Subset significant genes with |log2FC| ≥ 1.5 and padj ≤ 0.01
esus_vs_ce_upgenes_1_5_shrink = subset(res_shrink_df, log2FoldChange >= 1.5 & padj <= 0.05)
esus_vs_ce_downgenes_1_5_shrink = subset(res_shrink_df, log2FoldChange <= -1.5 & padj <= 0.05)

esus_vs_ce_upgenes_1_shrink = subset(res_shrink_df, log2FoldChange >= 1 & padj <= 0.05)
esus_vs_ce_downgenes_1_shrink = subset(res_shrink_df, log2FoldChange <= -1 & padj <= 0.05)

esus_vs_ce_upgenes_0_5_shrink = subset(res_shrink_df, log2FoldChange >= 0.5 & padj <= 0.05)
esus_vs_ce_downgenes_0_5_shrink = subset(res_shrink_df, log2FoldChange <= -0.5 & padj <= 0.05)


write.csv(x = esus_vs_ce_upgenes_1_5_shrink, file = "esus_vs_ce_upgenes_1_5_shrink.csv")
write.csv(x = esus_vs_ce_downgenes_1_5_shrink, file = "esus_vs_ce_downgenes_1_5_shrink.csv")
write.csv(x = esus_vs_ce_upgenes_1_shrink, file = "esus_vs_ce_upgenes_1_shrink.csv")
write.csv(x = esus_vs_ce_downgenes_1_shrink, file = "esus_vs_ce_downgenes_1_shrink.csv")
write.csv(x = esus_vs_ce_upgenes_0_5_shrink, file = "esus_vs_ce_upgenes_0_5_shrink.csv")
write.csv(x = esus_vs_ce_downgenes_0_5_shrink, file = "esus_vs_ce_downgenes_0_5_shrink.csv")

res = results(dds, contrast = c("Condition", "CE", "ESUS"))
res = res[order(res$padj, na.last = NA), ]
res_df = as.data.frame(res)
res_p_val_cut = subset(res_df, res_df$padj <= 0.05)
write.csv(res_p_val_cut, file = "ce_vs_esus_all_degs.csv")
esus_vs_ce_upgenes_1_5 = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
esus_vs_ce_downgenes_1_5 = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)

esus_vs_ce_upgenes_1 = subset(res_df, res_df$log2FoldChange >= 1 & res_df$padj <= 0.01)
esus_vs_ce_downgenes_1 = subset(res_df, res_df$log2FoldChange <= -1 & res_df$padj <= 0.01)

esus_vs_ce_upgenes_0_5 = subset(res_df, res_df$log2FoldChange >= 0.5 & res_df$padj <= 0.01)
esus_vs_ce_downgenes_0_5 = subset(res_df, res_df$log2FoldChange <= -0.5 & res_df$padj <= 0.01)

write.csv(esus_vs_ce_upgenes_1_5, file = "ce_vs_esus_upgenes_1_5.csv", row.names = T)
write.csv(esus_vs_ce_downgenes_1_5, file = "ce_vs_esus_downgenes_1_5.csv", row.names = T)

write.csv(esus_vs_ce_upgenes_1, file = "ce_vs_esus_upgenes_1.csv", row.names = T)
write.csv(esus_vs_ce_downgenes_1, file = "ce_vs_esus_downgenes_1.csv", row.names = T)

write.csv(esus_vs_ce_upgenes_0_5, file = "ce_vs_esus_upgenes_0_5.csv", row.names = T)
write.csv(esus_vs_ce_downgenes_0_5, file = "ce_vs_esus_downgenes_0_5.csv", row.names = T)

# vsd = vst(dds, blind = FALSE, fitType = 'mean')
# # plotPCA(vsd, intgroup = "Condition") + ggtitle("PCA After Normalization")
# 
# # Extract PCA data with variance explained
# pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# pca_data$Sample <- rownames(esus_vs_ce_metadata)  # Assign sample names from metadata
# 
# # Perform PCA manually to get variance explained
# pca <- prcomp(t(assay(vsd)))
# percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
# 
# # Add percent variance to axis labels
# ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
#   geom_point(size = 3) +
#   geom_label(aes(label = Sample), vjust = 1.5, size = 3) +
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#   ggtitle("PCA After Normalization: ESUS vs CE") +
#   theme_minimal()


# AT vs ESUS:
at_esus = x2[,c(0:6, 17:21)]
at_esus = at_esus[, !colnames(at_esus) %in% c("esus1", "at6", "at5", "at1")]
write.csv(at_esus, file = "at_esus_final_counts_matrix.csv")

# at_esus = read.csv(file = "at_esus_final_counts_matrix.csv")
# rownames(at_esus) = at_esus$X
# at_esus$X = NULL

at_vs_esus_metadata = read.csv("at_vs_esus.csv", row.names = 1)
dds = DESeqDataSetFromMatrix(countData = at_esus,
                             colData = at_vs_esus_metadata,
                             design = ~ Condition)

dds = dds[rowSums(counts(dds)) > 10, ]
nrow(dds)


vsd <- vst(dds, blind = TRUE)  # Variance stabilizing transformation

pca_data = plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# Extract percentage variance explained
percentVar = round(100 * attr(pca_data, "percentVar"))

# Add sample names
# Create PCA plot with sample names and variance info
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA After Normalization") +
  theme_minimal()














# vsd_raw <- vst(dds, blind = TRUE)  # Variance stabilizing transformation
# plotPCA(vsd_raw, intgroup = "Condition") + ggtitle("PCA Before Normalization")
dds$Condition <- relevel(dds$Condition, ref = "ESUS")  # Set CE as reference
dds = DESeq(dds)
vsd = vst(dds, blind = FALSE, fitType = "mean")
# Compute PCA data
pca_data = plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
# Extract percentage variance explained
percentVar = round(100 * attr(pca_data, "percentVar"))

# Add sample names
# Create PCA plot with sample names and variance info
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel() +  # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +  # Label PC1 with variance
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +  # Label PC2 with variance
  ggtitle("PCA AT vs ESUS") +
  theme_minimal()
resultsNames(dds)

# Perform LFC shrinkage
res_shrink = lfcShrink(dds, coef = "Condition_AT_vs_ESUS", type="apeglm")

# Convert to data frame
# res_shrink_df = as.data.frame(res_shrink)
# Subset significant genes with |log2FC| ≥ 1.5 and padj ≤ 0.01
# at_vs_esus_upgenes_1_5_shrink = subset(res_shrink_df, log2FoldChange >= 1.5 & padj <= 0.05)
# at_vs_esus_downgenes_1_5_shrink = subset(res_shrink_df, log2FoldChange <= -1.5 & padj <= 0.05)
# 
# at_vs_esus_upgenes_1_shrink = subset(res_shrink_df, log2FoldChange >= 1 & padj <= 0.05)
# at_vs_esus_downgenes_1_shrink = subset(res_shrink_df, log2FoldChange <= -1 & padj <= 0.05)
# 
# at_vs_esus_upgenes_0_5_shrink = subset(res_shrink_df, log2FoldChange >= 0.5 & padj <= 0.05)
# at_vs_esus_downgenes_0_5_shrink = subset(res_shrink_df, log2FoldChange <= -0.5 & padj <= 0.05)

res = results(dds, contrast = c("Condition", "ESUS", "AT"))
res = res[order(res$padj, na.last = NA), ]
res_df = as.data.frame(res)
res_p_val_cut = subset(res_df, res_df$padj <= 0.05)
write.csv(res_p_val_cut, file = "at_vs_esus_all_degs.csv")
esus_vs_at_upgenes_1_5 = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$padj <= 0.01)
esus_vs_at_downgenes_1_5 = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$padj <= 0.01)

esus_vs_at_upgenes_1 = subset(res_df, res_df$log2FoldChange >= 1 & res_df$padj <= 0.01)
esus_vs_at_downgenes_1 = subset(res_df, res_df$log2FoldChange <= -1 & res_df$padj <= 0.01)

esus_vs_at_upgenes_0_5 = subset(res_df, res_df$log2FoldChange >= 0.5 & res_df$padj <= 0.01)
esus_vs_at_downgenes_0_5 = subset(res_df, res_df$log2FoldChange <= -0.5 & res_df$padj <= 0.01)

write.csv(esus_vs_at_upgenes_1_5, file = "at_vs_esus_upgenes_1_5.csv", row.names = T)
write.csv(esus_vs_at_downgenes_1_5, file = "at_vs_esus_downgenes_1_5.csv", row.names = T)

write.csv(esus_vs_at_upgenes_1, file = "at_vs_esus_upgenes_1.csv", row.names = T)
write.csv(esus_vs_at_downgenes_1, file = "at_vs_esus_downgenes_1.csv", row.names = T)

write.csv(esus_vs_at_upgenes_0_5, file = "at_vs_esus_upgenes_0_5.csv", row.names = T)
write.csv(esus_vs_at_downgenes_0_5, file = "at_vs_esus_downgenes_0_5.csv", row.names = T)

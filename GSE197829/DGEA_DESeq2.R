# DeSeq2 Pipeline:

setwd("/home/ibab/HM_VM_Project/GSE197829/star_htseq_counts")

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

counts_matrix = read.csv(file = "GSE187829_final_counts_matrix_coding.csv", row.names = 1)

# HC vs LAA:
hc_vs_laa = counts_matrix %>% dplyr::select(1:6)
hc_vs_laa_metadata = read.csv(file = "hc_vs_laa_metadata.csv", row.names = 1)

dds = DESeqDataSetFromMatrix(countData = hc_vs_laa,
                             colData = hc_vs_laa_metadata,
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
  ggtitle("PCA After Normalization: Healthy Control vs Large Artery Atherosclerosis") +
  theme_minimal()
levels(dds$Condition)


dds = DESeq(dds)
resultsNames(dds)
# Extract normalized counts
normalized_counts = counts(dds, normalized = TRUE)

# Convert to data frame
normalized_counts_df = as.data.frame(normalized_counts)

res = results(dds, contrast = c("Condition", "LAA", "Control"))
res = res[order(res$pvalue, na.last = NA), ]
res_df = as.data.frame(res)
res_p_val_cut = subset(res_df, res_df$pvalue <= 0.05)
write.csv(res_p_val_cut, file = "laa_vs_ctrl_all_genes.csv", row.names = T)


laa_vs_ctrl_upgenes_1_5 = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$pvalue <= 0.05)
laa_vs_ctrl_downgenes_1_5 = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$pvalue <= 0.05)

laa_vs_ctrl_upgenes_1 = subset(res_df, res_df$log2FoldChange >= 1 & res_df$pvalue <= 0.05)
laa_vs_ctrl_downgenes_1 = subset(res_df, res_df$log2FoldChange <= -1 & res_df$pvalue <= 0.05)

laa_vs_ctrl_upgenes_0_5 = subset(res_df, res_df$log2FoldChange >= 0.5 & res_df$pvalue <= 0.05)
laa_vs_ctrl_downgenes_0_5 = subset(res_df, res_df$log2FoldChange <= -0.5 & res_df$pvalue <= 0.05)

write.csv(laa_vs_ctrl_upgenes_1_5, file = "laa_vs_ctrl_upgenes_1_5_deseq.csv", row.names = T)
write.csv(laa_vs_ctrl_downgenes_1_5, file = "laa_vs_ctrl_downgenes_1_5_deseq.csv", row.names = T)

write.csv(laa_vs_ctrl_upgenes_1, file = "laa_vs_ctrl_upgenes_1_deseq.csv", row.names = T)
write.csv(laa_vs_ctrl_downgenes_1, file = "laa_vs_ctrl_downgenes_1_deseq.csv", row.names = T)

write.csv(laa_vs_ctrl_upgenes_0_5, file = "laa_vs_ctrl_upgenes_0_5_deseq.csv", row.names = T)
write.csv(laa_vs_ctrl_downgenes_0_5, file = "laa_vs_ctrl_downgenes_0_5_deseq.csv", row.names = T)

# Get top 25 upregulated genes
top_25_up_laa_ctrl <- res_p_val_cut %>%
  filter(log2FoldChange >= 1.5 & pvalue <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_laa_ctrl <- res_p_val_cut %>%
  filter(log2FoldChange <= -1.5 & pvalue <= 0.05) %>%
  arrange(log2FoldChange) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts_df, rownames(normalized_counts_df) %in% 
                                   rownames(top_25_up_laa_ctrl))

normalized_counts_df_down = subset(normalized_counts_df, rownames(normalized_counts_df) %in% 
                                     rownames(top_25_down_laa_ctrl))
library(pheatmap)

# Convert data frame to matrix (if not already in matrix format)
heatmap_data_up <- as.matrix(normalized_counts_df_up)

# Create the heatmap
pheatmap(heatmap_data_up, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("HC1", "HC2", "HC3", "LAA1", "LAA2", "LAA3"),
         main = "Heatmap of Top 25 Upregulated Genes: Control vs Large Artery Atherosclerosis")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("HC1", "HC2", "HC3", "LAA1", "LAA2", "LAA3"),
         main = "Heatmap of Top 25 Downregulated Genes: Control vs Large Artery Atherosclerosis")


res_df$neg_logp_val = -log10(res_df$pvalue)
res_df$Expression = ifelse(res_df$log2FoldChange >= 0.5 & res_df$pvalue <= 0.05, "Upregulated",
                                 ifelse(res_df$log2FoldChange <= -0.5 & 
                                          res_df$pvalue <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(res_df, aes(x=log2FoldChange, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Large Artery Atherosclerosis")












# HC vs SAA:
hc_vs_saa = counts_matrix %>% dplyr::select(1:3, 7:9)
hc_vs_saa_metadata = read.csv(file = "hc_vs_saa_metadata.csv", row.names = 1)

dds = DESeqDataSetFromMatrix(countData = hc_vs_saa,
                             colData = hc_vs_saa_metadata,
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
  ggtitle("PCA After Normalization: Healthy Control vs Small Artery Atherosclerosis") +
  theme_minimal()
levels(dds$Condition)


dds = DESeq(dds)
resultsNames(dds)
# Extract normalized counts
normalized_counts = counts(dds, normalized = TRUE)

# Convert to data frame
normalized_counts_df = as.data.frame(normalized_counts)

res = results(dds, contrast = c("Condition", "SAA", "Control"))
res = res[order(res$pvalue, na.last = NA), ]
res_df = as.data.frame(res)
res_p_val_cut = subset(res_df, res_df$pvalue <= 0.05)
write.csv(res_p_val_cut, file = "saa_vs_ctrl_all_genes.csv", row.names = T)


saa_vs_ctrl_upgenes_1_5 = subset(res_df, res_df$log2FoldChange >= 1.5 & res_df$pvalue <= 0.05)
saa_vs_ctrl_downgenes_1_5 = subset(res_df, res_df$log2FoldChange <= -1.5 & res_df$pvalue <= 0.05)

saa_vs_ctrl_upgenes_1 = subset(res_df, res_df$log2FoldChange >= 1 & res_df$pvalue <= 0.05)
saa_vs_ctrl_downgenes_1 = subset(res_df, res_df$log2FoldChange <= -1 & res_df$pvalue <= 0.05)

saa_vs_ctrl_upgenes_0_5 = subset(res_df, res_df$log2FoldChange >= 0.5 & res_df$pvalue <= 0.05)
saa_vs_ctrl_downgenes_0_5 = subset(res_df, res_df$log2FoldChange <= -0.5 & res_df$pvalue <= 0.05)

write.csv(saa_vs_ctrl_upgenes_1_5, file = "saa_vs_ctrl_upgenes_1_5_deseq.csv", row.names = T)
write.csv(saa_vs_ctrl_downgenes_1_5, file = "saa_vs_ctrl_downgenes_1_5_deseq.csv", row.names = T)

write.csv(saa_vs_ctrl_upgenes_1, file = "saa_vs_ctrl_upgenes_1_deseq.csv", row.names = T)
write.csv(saa_vs_ctrl_downgenes_1, file = "saa_vs_ctrl_downgenes_1_deseq.csv", row.names = T)

write.csv(saa_vs_ctrl_upgenes_0_5, file = "saa_vs_ctrl_upgenes_0_5_deseq.csv", row.names = T)
write.csv(saa_vs_ctrl_downgenes_0_5, file = "saa_vs_ctrl_downgenes_0_5_deseq.csv", row.names = T)

# Get top 25 upregulated genes
top_25_up_saa_ctrl <- res_p_val_cut %>%
  filter(log2FoldChange >= 1.5 & pvalue <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_saa_ctrl <- res_p_val_cut %>%
  filter(log2FoldChange <= -1.5 & pvalue <= 0.05) %>%
  arrange(log2FoldChange) %>%
  head(25)

normalized_counts_df_up = subset(normalized_counts_df, rownames(normalized_counts_df) %in% 
                                   rownames(top_25_up_saa_ctrl))

normalized_counts_df_down = subset(normalized_counts_df, rownames(normalized_counts_df) %in% 
                                     rownames(top_25_down_saa_ctrl))
library(pheatmap)

# Convert data frame to matrix (if not already in matrix format)
heatmap_data_up <- as.matrix(normalized_counts_df_up)

# Create the heatmap
pheatmap(heatmap_data_up, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("HC1", "HC2", "HC3", "SAA1", "SAA2", "SAA3"),
         main = "Heatmap of Top 25 Upregulated Genes: Control vs Small Artery Atherosclerosis")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("HC1", "HC2", "HC3", "SAA1", "SAA2", "SAA3"),
         main = "Heatmap of Top 25 Downregulated Genes: Control vs Small Artery Atherosclerosis")


res_df$neg_logp_val = -log10(res_df$pvalue)
res_df$Expression = ifelse(res_df$log2FoldChange >= 0.5 & res_df$pvalue <= 0.05, "Upregulated",
                                  ifelse(res_df$log2FoldChange <= -0.5 & 
                                           res_df$pvalue <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(res_df, aes(x=log2FoldChange, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Small Artery Atherosclerosis")




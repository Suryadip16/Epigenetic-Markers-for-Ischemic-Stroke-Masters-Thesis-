library(dplyr)

setwd("/home/ibab/HM_VM_Project/array/GSE58294_(copy)/output_files/")

gene_lengths = read.csv(file = "/home/ibab/HM_VM_Project/hg38_main/hg38_ensGene_gene_lengths.csv")
counts_matrix = read.csv(file = "GSE187829_final_counts_matrix_coding.csv", row.names = 1)



# SAA vs Ctrl: 
saa_ctrl_counts = counts_matrix %>% dplyr::select(hc1:hc3, saa1:saa3)

# Convert gene_lengths dataframe to a named vector
gene_lengths_vector <- setNames(gene_lengths$Length, gene_lengths$hgnc_symbol)

# Ensure counts_matrix row names match gene IDs in gene_lengths_vector
matching_genes <- intersect(rownames(saa_ctrl_counts), names(gene_lengths_vector))
saa_ctrl_filtered <- saa_ctrl_counts[matching_genes, ]
gene_lengths_filtered <- gene_lengths_vector[matching_genes]

# Convert gene lengths to kilobases
gene_lengths_kb <- gene_lengths_filtered / 1000

# Compute RPK (Reads Per Kilobase)
rpk <- saa_ctrl_filtered / gene_lengths_kb

# Compute TPM using edgeR
library(edgeR)
saa_ctrl_tpm_matrix <- cpm(rpk, log = FALSE)

# View TPM-normalized matrix
head(saa_ctrl_tpm_matrix)


saa_ctrl_genes = read.csv(file = "saa_vs_ctrl_all_genes.csv", row.names = 1)
epifactors_list = read.csv(file = "/home/ibab/HM_VM_Project/Epigenes_epicypher_final.csv", row.names = 1)
saa_ctrl_intersection_genes = subset(saa_ctrl_genes, rownames(saa_ctrl_genes) %in% rownames(epifactors_list))
epifactors_list$rownames = rownames(epifactors_list)
saa_ctrl_intersection_genes$rownames = rownames(saa_ctrl_intersection_genes)
saa_ctrl_merged = merge(saa_ctrl_intersection_genes, epifactors_list, by = "rownames", all.x = T)
rownames(saa_ctrl_merged) = saa_ctrl_merged$rownames
saa_ctrl_merged$rownames = NULL

write.csv(saa_ctrl_merged, file = "saa_ctrl_epifactors.csv", row.names = T)

saa_ctrl_tpm_epi = subset(saa_ctrl_tpm_matrix, rownames(saa_ctrl_tpm_matrix) %in% rownames(epifactors_list))


# Get top 25 upregulated epifactors
top_25_up_saa_ctrl_epi <- saa_ctrl_merged %>%
  filter(log2FoldChange > 0.0 & pvalue <= 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  head(25)

# Get top 25 downregulated epifactors
top_25_down_saa_ctrl_epi <- saa_ctrl_merged %>%
  filter(log2FoldChange < 0.0 & pvalue <= 0.05) %>%
  arrange(log2FoldChange) %>%
  head(16)

normalized_counts_df_up = subset(saa_ctrl_tpm_epi, rownames(saa_ctrl_tpm_epi) %in% 
                                   rownames(top_25_up_saa_ctrl_epi))
# normalized_counts_df_up <- normalized_counts_df_up[!(rownames(normalized_counts_df_up) %in% "PHF13"), ]

normalized_counts_df_down = subset(saa_ctrl_tpm_epi, rownames(saa_ctrl_tpm_epi) %in% 
                                     rownames(top_25_down_saa_ctrl_epi))

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
         main = "Heatmap of Upregulated Epigenetic Factors\nSAA vs Control")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("HC1", "HC2", "HC3", "SAA1", "SAA2", "SAA3"),
         main = "Heatmap of Downregulated Epigenetic Factors\nSAA vs Control")




# CE vs Ctrl:
ce_ctrl_counts = read.csv(file = "ce_ctrl_expr_matrix_1.csv", row.names = 1) # Already RMA normalised

# Read DE results and epifactors list
ce_ctrl_genes = read.csv(file = "ce_vs_ctrl_all_genes_pval_cut.csv", row.names = 1)  # assuming this is DE results
epifactors_list = read.csv(file = "/home/ibab/HM_VM_Project/Epigenes_epicypher_final.csv", row.names = 1)

# Filter and merge with epigenetic factors
ce_ctrl_intersection_genes = subset(ce_ctrl_genes, rownames(ce_ctrl_genes) %in% rownames(epifactors_list))
epifactors_list$rownames = rownames(epifactors_list)
ce_ctrl_intersection_genes$rownames = rownames(ce_ctrl_intersection_genes)
ce_ctrl_merged = merge(ce_ctrl_intersection_genes, epifactors_list, by = "rownames", all.x = T)
rownames(ce_ctrl_merged) = ce_ctrl_merged$rownames
ce_ctrl_merged$rownames = NULL

write.csv(ce_ctrl_merged, file = "ce_ctrl_epifactors.csv", row.names = T)

# Subset normalized counts for epigenetic factors
ce_ctrl_norm_epi = subset(ce_ctrl_counts, rownames(ce_ctrl_counts) %in% rownames(epifactors_list))

# Get top 25 upregulated epifactors
library(dplyr)
top_25_up_ce_ctrl_epi <- ce_ctrl_merged %>%
  filter(logFC >= 0.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated epifactors
top_25_down_ce_ctrl_epi <- ce_ctrl_merged %>%
  filter(logFC <= -0.5 & adj.P.Val <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

# Subset normalized counts for top epifactors
normalized_counts_df_up = subset(ce_ctrl_norm_epi, rownames(ce_ctrl_norm_epi) %in% 
                                   rownames(top_25_up_ce_ctrl_epi))
normalized_counts_df_down = subset(ce_ctrl_norm_epi, rownames(ce_ctrl_norm_epi) %in% 
                                     rownames(top_25_down_ce_ctrl_epi))

# Convert data frames to matrices
heatmap_data_up <- as.matrix(normalized_counts_df_up)
heatmap_data_down <- as.matrix(normalized_counts_df_down)

# Create heatmaps
library(pheatmap)
group_df = data.frame(Groups=rep(c("HC", "CE"), c(23,69)))
rownames(group_df) <- colnames(heatmap_data_down)

ann_colors = list(
  Groups = c(Control="turquoise", Ischemic="darkred"))

# Plot heatmap with sample type annotations
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         labels_col = c("HC1", "HC2", "HC3", "HC4", "HC5", "HC6", "HC7", "HC8",
                        "HC9", "HC10", "HC11", "HC12", "HC13", "HC14", "HC15",
                        "HC16", "HC17", "HC18", "HC19", "HC20", "HC21", "HC22", "HC23", 
                        "CE1", "CE2", "CE3", "CE4", "CE5", "CE6", "CE7", "CE8", "CE9", 
                        "CE10", "CE11", "CE12", "CE13", "CE14", "CE15", "CE16", "CE17", 
                        "CE18", "CE19", "CE20", "CE21", "CE22", "CE23", "CE24", "CE25", 
                        "CE26", "CE27", "CE28", "CE29", "CE30", "CE31", "CE32", "CE33", 
                        "CE34", "CE35", "CE36", "CE37", "CE38", "CE39", "CE40", "CE41", 
                        "CE42", "CE43", "CE44", "CE45", "CE46", "CE47", "CE48", "CE49", 
                        "CE50", "CE51", "CE52", "CE53", "CE54", "CE55", "CE56", "CE57", 
                        "CE58", "CE59", "CE60", "CE61", "CE62", "CE63", "CE64", "CE65", 
                        "CE66", "CE67", "CE68", "CE69"),
         main = "Heatmap of Top Upregulated Epifactors",
         angle_col = 45,
         fontsize_col = 4,
         annotation_col = group_df,
         show_legend = TRUE)

group_df = data.frame(Groups=rep(c("HC", "CE"), c(23,69)))
rownames(group_df) <- colnames(heatmap_data_up)

# Plot heatmap with sample type annotations
pheatmap(heatmap_data_up, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         labels_col = c("HC1", "HC2", "HC3", "HC4", "HC5", "HC6", "HC7", "HC8",
                        "HC9", "HC10", "HC11", "HC12", "HC13", "HC14", "HC15",
                        "HC16", "HC17", "HC18", "HC19", "HC20", "HC21", "HC22", "HC23", 
                        "CE1", "CE2", "CE3", "CE4", "CE5", "CE6", "CE7", "CE8", "CE9", 
                        "CE10", "CE11", "CE12", "CE13", "CE14", "CE15", "CE16", "CE17", 
                        "CE18", "CE19", "CE20", "CE21", "CE22", "CE23", "CE24", "CE25", 
                        "CE26", "CE27", "CE28", "CE29", "CE30", "CE31", "CE32", "CE33", 
                        "CE34", "CE35", "CE36", "CE37", "CE38", "CE39", "CE40", "CE41", 
                        "CE42", "CE43", "CE44", "CE45", "CE46", "CE47", "CE48", "CE49", 
                        "CE50", "CE51", "CE52", "CE53", "CE54", "CE55", "CE56", "CE57", 
                        "CE58", "CE59", "CE60", "CE61", "CE62", "CE63", "CE64", "CE65", 
                        "CE66", "CE67", "CE68", "CE69"),
         main = "Heatmap of Top Downregulated Epifactors",
         angle_col = 45,
         fontsize_col = 4,
         annotation_col = group_df,
         show_legend = TRUE)


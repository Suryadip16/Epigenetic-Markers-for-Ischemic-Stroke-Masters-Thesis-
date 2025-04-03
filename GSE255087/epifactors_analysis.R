library(dplyr)

setwd("/home/ibab/HM_VM_Project/DGEA/edgeR_res/pea/edgeR_genes")

gene_lengths = read.csv(file = "/home/ibab/HM_VM_Project/hg38_main/hg38_ensGene_gene_lengths.csv")
counts_matrix = read.csv(file = "/home/ibab/HM_VM_Project/DGEA/star_htseq_counts/processed_final_counts_matrix_coding.csv", row.names = 1)

# AT vs CE: 
at_ce_counts = counts_matrix %>% dplyr::select(at1:at5, ce1, ce3, ce8, ce9)

# Convert gene_lengths dataframe to a named vector
gene_lengths_vector <- setNames(gene_lengths$Length, gene_lengths$hgnc_symbol)

# Ensure counts_matrix row names match gene IDs in gene_lengths_vector
matching_genes <- intersect(rownames(at_ce_counts), names(gene_lengths_vector))
at_ce_filtered <- at_ce_counts[matching_genes, ]
gene_lengths_filtered <- gene_lengths_vector[matching_genes]

# Convert gene lengths to kilobases
gene_lengths_kb <- gene_lengths_filtered / 1000

# Compute RPK (Reads Per Kilobase)
rpk <- at_ce_filtered / gene_lengths_kb

# Compute TPM using edgeR
library(edgeR)
at_ce_tpm_matrix <- cpm(rpk, log = FALSE)

# View TPM-normalized matrix
head(at_ce_tpm_matrix)


at_ce_genes = read.csv(file = "new_at_vs_ce_edgeR_all_genes.csv", row.names = 1)
epifactors_list = read.csv(file = "/home/ibab/HM_VM_Project/Epigenes_epicypher_final.csv", row.names = 1)
at_ce_intersection_genes = subset(at_ce_genes, rownames(at_ce_genes) %in% rownames(epifactors_list))
epifactors_list$rownames = rownames(epifactors_list)
at_ce_intersection_genes$rownames = rownames(at_ce_intersection_genes)
at_ce_merged = merge(at_ce_intersection_genes, epifactors_list, by = "rownames", all.x = T)
rownames(at_ce_merged) = at_ce_merged$rownames
at_ce_merged$rownames = NULL

write.csv(at_ce_merged, file = "at_ce_epifactors.csv", row.names = T)

at_ce_tpm_epi = subset(at_ce_tpm_matrix, rownames(at_ce_tpm_matrix) %in% rownames(epifactors_list))


# Get top 25 upregulated epifactors
top_25_up_at_ce_epi <- at_ce_merged %>%
  filter(logFC >= 1.0 & PValue <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(15)

# Get top 25 downregulated epifactors
top_25_down_at_ce_epi <- at_ce_merged %>%
  filter(logFC <= -1.0 & PValue <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(at_ce_tpm_epi, rownames(at_ce_tpm_epi) %in% 
                                   rownames(top_25_up_at_ce_epi))

normalized_counts_df_down = subset(at_ce_tpm_epi, rownames(at_ce_tpm_epi) %in% 
                                     rownames(top_25_down_at_ce_epi))

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
         labels_col = c("AT1", "AT2", "AT3", "AT4", "AT5", "CE1", "CE3", "CE8", "CE9"),
         main = "Heatmap of Top 15 Upregulated Epigenetic Factors:\nAtherothrombotic vs Cardioembolic Stroke")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("AT1", "AT2", "AT3", "AT4", "AT5", "CE1", "CE3", "CE8", "CE9"),
         main = "Heatmap of Top 15 Downregulated Epigenetic Factors:\nAtherothrombotic vs Cardioembolic Stroke")














# AT vs ESUS: 
at_esus_counts = counts_matrix %>% dplyr::select(at1:at2, at4, at5, esus4, esus5)

# Convert gene_lengths dataframe to a named vector
gene_lengths_vector <- setNames(gene_lengths$Length, gene_lengths$hgnc_symbol)

# Ensure counts_matrix row names match gene IDs in gene_lengths_vector
matching_genes <- intersect(rownames(at_esus_counts), names(gene_lengths_vector))
at_esus_filtered <- at_esus_counts[matching_genes, ]
gene_lengths_filtered <- gene_lengths_vector[matching_genes]

# Convert gene lengths to kilobases
gene_lengths_kb <- gene_lengths_filtered / 1000

# Compute RPK (Reads Per Kilobase)
rpk <- at_esus_filtered / gene_lengths_kb

# Compute TPM using edgeR
library(edgeR)
at_esus_tpm_matrix <- cpm(rpk, log = FALSE)

# View TPM-normalized matrix
head(at_esus_tpm_matrix)


at_esus_genes = read.csv(file = "new_at_vs_esus_edgeR_all_genes.csv", row.names = 1)
epifactors_list = read.csv(file = "/home/ibab/HM_VM_Project/Epigenes_epicypher_final.csv", row.names = 1)
at_esus_intersection_genes = subset(at_esus_genes, rownames(at_esus_genes) %in% rownames(epifactors_list))
epifactors_list$rownames = rownames(epifactors_list)
at_esus_intersection_genes$rownames = rownames(at_esus_intersection_genes)
at_esus_merged = merge(at_esus_intersection_genes, epifactors_list, by = "rownames", all.x = T)
rownames(at_esus_merged) = at_esus_merged$rownames
at_esus_merged$rownames = NULL

write.csv(at_esus_merged, file = "at_esus_epifactors.csv", row.names = T)

at_esus_tpm_epi = subset(at_esus_tpm_matrix, rownames(at_esus_tpm_matrix) %in% rownames(epifactors_list))


# Get top 25 upregulated epifactors
top_25_up_at_esus_epi <- at_esus_merged %>%
  filter(logFC >= 1.0 & PValue <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(15)

# Get top 25 downregulated epifactors
top_25_down_at_esus_epi <- at_esus_merged %>%
  filter(logFC <= -1.0 & PValue <= 0.05) %>%
  arrange(logFC) %>%
  head(15)

normalized_counts_df_up = subset(at_esus_tpm_epi, rownames(at_esus_tpm_epi) %in% 
                                   rownames(top_25_up_at_esus_epi))
normalized_counts_df_up <- normalized_counts_df_up[!(rownames(normalized_counts_df_up) %in% "PHF13"), ]

normalized_counts_df_down = subset(at_esus_tpm_epi, rownames(at_esus_tpm_epi) %in% 
                                     rownames(top_25_down_at_esus_epi))

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
         labels_col = c("AT1", "AT2", "AT4", "AT5", "ESUS4", "ESUS5"),
         main = "Heatmap of Top 15 Upregulated Epigenetic Factors:\nAtherothrombotic vs Embolic Stroke of Undetermined Source")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("AT1", "AT2", "AT4", "AT5", "ESUS4", "ESUS5"),
         main = "Heatmap of Top 15 Downregulated Epigenetic Factors:\nAtherothrombotic vs Embolic Stroke of Undetermined Source")












# CE vs ESUS: 
ce_esus_counts = counts_matrix %>% dplyr::select(ce2, ce4:ce7, esus2:esus5)

# Convert gene_lengths dataframe to a named vector
gene_lengths_vector <- setNames(gene_lengths$Length, gene_lengths$hgnc_symbol)

# Ensure counts_matrix row names match gene IDs in gene_lengths_vector
matching_genes <- intersect(rownames(ce_esus_counts), names(gene_lengths_vector))
ce_esus_filtered <- ce_esus_counts[matching_genes, ]
gene_lengths_filtered <- gene_lengths_vector[matching_genes]

# Convert gene lengths to kilobases
gene_lengths_kb <- gene_lengths_filtered / 1000

# Compute RPK (Reads Per Kilobase)
rpk <- ce_esus_filtered / gene_lengths_kb

# Compute TPM using edgeR
library(edgeR)
ce_esus_tpm_matrix <- cpm(rpk, log = FALSE)

# View TPM-normalized matrix
head(ce_esus_tpm_matrix)


ce_esus_genes = read.csv(file = "new_ce_vs_esus_edgeR_all_genes.csv", row.names = 1)
epifactors_list = read.csv(file = "/home/ibab/HM_VM_Project/Epigenes_epicypher_final.csv", row.names = 1)
ce_esus_intersection_genes = subset(ce_esus_genes, rownames(ce_esus_genes) %in% rownames(epifactors_list))
epifactors_list$rownames = rownames(epifactors_list)
ce_esus_intersection_genes$rownames = rownames(ce_esus_intersection_genes)
ce_esus_merged = merge(ce_esus_intersection_genes, epifactors_list, by = "rownames", all.x = T)
rownames(ce_esus_merged) = ce_esus_merged$rownames
ce_esus_merged$rownames = NULL

write.csv(ce_esus_merged, file = "ce_esus_epifactors.csv", row.names = T)

ce_esus_tpm_epi = subset(ce_esus_tpm_matrix, rownames(ce_esus_tpm_matrix) %in% rownames(epifactors_list))


# Get top 25 upregulated epifactors
top_25_up_ce_esus_epi <- ce_esus_merged %>%
  filter(logFC >= 1.0 & PValue <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated epifactors
top_25_down_ce_esus_epi <- ce_esus_merged %>%
  filter(logFC <= -1.0 & PValue <= 0.05) %>%
  arrange(logFC) %>%
  head(25)

normalized_counts_df_up = subset(ce_esus_tpm_epi, rownames(ce_esus_tpm_epi) %in% 
                                   rownames(top_25_up_ce_esus_epi))


normalized_counts_df_down = subset(ce_esus_tpm_epi, rownames(ce_esus_tpm_epi) %in% 
                                     rownames(top_25_down_ce_esus_epi))

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
         labels_col = c("CE2", "CE4", "CE5", "CE6", "CE7","ESUS2", "ESUS3", "ESUS4", "ESUS5"),
         main = "Heatmap of Top 25 Upregulated Epigenetic Factors:\nCardioembolic vs Embolic Stroke of Undetermined Source")


heatmap_data_down <- as.matrix(normalized_counts_df_down)
pheatmap(heatmap_data_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100), # Gradient colors
         cluster_rows = T,  # Cluster rows
         cluster_cols = T,  # Cluster columns
         scale = "row",        # Normalize by row (optional)
         show_rownames = TRUE,  # Show row labels
         show_colnames = TRUE,  # Show column labels
         labels_col = c("CE2", "CE4", "CE5", "CE6", "CE7","ESUS2", "ESUS3", "ESUS4", "ESUS5"),
         main = "Heatmap of Top 25 Upregulated Epigenetic Factors:\nCardioembolic vs Embolic Stroke of Undetermined Source")

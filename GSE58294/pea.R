# PEA using clusterProfiler

# Credits: https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/


# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations

# Set input path
in_path = "C:/Users/dipak/Desktop/Riku/VM2_Project/edgeR_genes/"
out_path = "C:/Users/dipak/Desktop/Riku/VM2_Project/GSEA/PEA/"
bg_path = "C:/Users/dipak/Desktop/Riku/VM2_Project/GSEA/Background_genes/"
setwd(in_path)


  


df = read.csv("ce_esus_common_genes_final.csv")


# Get the genes that are present in your dataframe
genes_in_data = df$X


# Convert gene symbols to Entrez IDs
entrez_ids = bitr(genes_in_data, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_ids = entrez_ids$ENTREZID

# Perform Enrichment

kegg_result = enrichKEGG(gene = gene_ids, 
                          organism = "hsa", 
                          pvalueCutoff = 0.05)
write.csv(as.data.frame(kegg_result@result), file = paste0(out_path, "ce_esus_KEGG_Enrichment_Results.csv"), row.names = FALSE)


go_result = enrichGO(gene = gene_ids, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "BP", 
                      pvalueCutoff = 0.05)
write.csv(as.data.frame(go_result@result), file = paste0(out_path, "ce_esus_GO_CC_Enrichment_Results.csv"), row.names = FALSE)

# Visualise GO Enrichment
go_res_df = go_result@result
significant_go = subset(go_res_df, go_res_df$pvalue <= 0.05)

top20_go = significant_go %>% 
  slice_max(order_by = FoldEnrichment, n = 20)

ggplot(top20_go, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +  # Adjust bar width
  labs(
    x = "Enrichment Score (ES)", 
    y = "Pathway", 
    title = "Enriched GO (CC): CE vs ESUS"
  ) +
  theme_minimal(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center title, bold
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.title = element_text(face = "bold"),  # Bold axis labels
    panel.grid.major.y = element_blank(),  # Remove horizontal gridlines for clarity
    panel.grid.minor = element_blank()
  ) +
  geom_text(aes(label = round(FoldEnrichment, 2)), hjust = -0.2, size = 4)


# Visualise Pathway Enrichment
kegg_res_df = kegg_result@result
significant_kegg = subset(kegg_res_df, kegg_res_df$pvalue <= 0.05)

ggplot(significant_kegg, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +  # Adjust bar width
  labs(
    x = "Enrichment Score (ES)", 
    y = "Pathway", 
    title = "Enriched Pathways (KEGG): CE vs ESUS"
  ) +
  theme_minimal(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Center title, bold
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.title = element_text(face = "bold"),  # Bold axis labels
    panel.grid.major.y = element_blank(),  # Remove horizontal gridlines for clarity
    panel.grid.minor = element_blank()
  ) +
  geom_text(aes(label = round(FoldEnrichment, 2)), hjust = -0.2, size = 4)




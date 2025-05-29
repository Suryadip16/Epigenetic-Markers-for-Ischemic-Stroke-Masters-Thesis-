setwd("/home/ibab/HM_VM_Project/array/GSE58294_(copy)/output_files/")

library(dplyr)

array_res = read.csv(file = "GSE58294_filtered_limma_output.csv")
array_res = na.omit(array_res)
array_res = array_res[!duplicated(array_res$SYMBOL), ]


# array_res$FoldChange = 2^array_res$logFC
array_res_pval_cut = subset(array_res, array_res$adj.P.Val <= 0.05)
array_res_pval_cut_clean = array_res_pval_cut %>% dplyr::select(SYMBOL, DESC, logFC, adj.P.Val, adj.P.Val)
array_res_pval_cut_clean = array_res_pval_cut_clean[!duplicated(array_res_pval_cut_clean$SYMBOL), ]
rownames(array_res_pval_cut_clean) = array_res_pval_cut_clean$SYMBOL
array_res_pval_cut_clean$SYMBOL = NULL
write.csv(array_res_pval_cut_clean, file = "is_vs_ctrl_all_genes_pval_cut_ml_data2.csv", row.names = T)
array_res_fc_filter = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC >= 0.5 | array_res_pval_cut_clean$logFC <= -0.5)
write.csv(array_res_fc_filter, file = "ce_ctrl_all_genes_fc_filter_pval_cut.csv", row.names = T)




ce_vs_ctrl_upgenes_1_5 = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC >= 1.5 & array_res_pval_cut_clean$adj.P.Val <= 0.05)
ce_vs_ctrl_downgenes_1_5 = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC <= -1.5 & array_res_pval_cut_clean$adj.P.Val <= 0.05)

ce_vs_ctrl_upgenes_1 = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC >= 1 & array_res_pval_cut_clean$adj.P.Val <= 0.05)
ce_vs_ctrl_downgenes_1 = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC <= -1 & array_res_pval_cut_clean$adj.P.Val <= 0.05)

ce_vs_ctrl_upgenes_0_5 = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC >= 0.5 & array_res_pval_cut_clean$adj.P.Val <= 0.05)
ce_vs_ctrl_downgenes_0_5 = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC <= -0.5 & array_res_pval_cut_clean$adj.P.Val <= 0.05)

write.csv(ce_vs_ctrl_upgenes_1_5, file = "ce_vs_ctrl_upgenes_1_5_limma3.csv", row.names = T)
write.csv(ce_vs_ctrl_downgenes_1_5, file = "ce_vs_ctrl_downgenes_1_5_limma3.csv", row.names = T)

write.csv(ce_vs_ctrl_upgenes_1, file = "ce_vs_ctrl_upgenes_1_limma3.csv", row.names = T)
write.csv(ce_vs_ctrl_downgenes_1, file = "ce_vs_ctrl_downgenes_1_limma3.csv", row.names = T)

write.csv(ce_vs_ctrl_upgenes_0_5, file = "ce_vs_ctrl_upgenes_0_5_limma3.csv", row.names = T)
write.csv(ce_vs_ctrl_downgenes_0_5, file = "ce_vs_ctrl_downgenes_0_5_limma3.csv", row.names = T)




# array_res_sig_lfc = subset(array_res_pval_cut_clean, array_res_pval_cut_clean$logFC >= 0.5 | array_res_pval_cut_clean$logFC <= -0.5)
# array_res_sig_lfc = array_res_sig_lfc[!duplicated(rownames(array_res_sig_lfc)), ]
# write.csv(array_res_sig_lfc, file = "ce_vs_ctrl_all_genes_fc_filter.csv", row.names = T)
# colnames(array_res_sig_lfc)







expr_matrix = array_res_pval_cut[ , c("SYMBOL", colnames(array_res_pval_cut)[13:ncol(array_res_pval_cut)]) ]
rownames(expr_matrix) = expr_matrix$SYMBOL
expr_matrix$SYMBOL = NULL
colnames(expr_matrix)

write.csv(expr_matrix, file = "is_ctrl_expr_matrix_ml_data2.csv", row.names = T)

# Manually compute PCA using prcomp for 3 PCs
pca_res = prcomp(t(as.matrix(expr_matrix)))
pca_var = pca_res$sdev^2
percentVar = round(100 * pca_var / sum(pca_var))[1:3]  # Get % variance for PC1–PC3
pc_df = as.data.frame(pca_res$x[, 1:3])
pc_df$Condition = factor(c(rep("HC", 23), rep("CE", 69)))
pc_df$Sample = rownames(pc_df)

library(plotly)
# Plot with plotly
plot_ly(pc_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Condition, type = 'scatter3d', mode = 'markers+text') %>%
  layout(title = '3D PCA (Normalized Counts): HC vs CE',
         scene = list(
           xaxis = list(title = paste0("PC1: ", percentVar[1], "%")),
           yaxis = list(title = paste0("PC2: ", percentVar[2], "%")),
           zaxis = list(title = paste0("PC3: ", percentVar[3], "%"))
         ))




# Get top 25 upregulated genes
top_25_up_ce_ctrl = array_res_pval_cut_clean %>%
  filter(logFC >= 0.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

# Get top 25 downregulated genes
top_25_down_ce_ctrl = array_res_pval_cut_clean %>%
  filter(logFC <= -0.5 & adj.P.Val <= 0.05) %>%
  arrange(desc(logFC)) %>%
  head(25)

normalized_counts_df_up = subset(expr_matrix, rownames(expr_matrix) %in% 
                                  rownames(top_25_up_ce_ctrl))

normalized_counts_df_up = normalized_counts_df_up[!grepl("LINC", rownames(normalized_counts_df_up)), ]


normalized_counts_df_down = subset(expr_matrix, rownames(expr_matrix) %in% 
                                     rownames(top_25_down_ce_ctrl))
library(pheatmap)

# Convert data frame to matrix (if not already in matrix format)
heatmap_data_up = as.matrix(normalized_counts_df_up)
group_df = data.frame(Groups=rep(c("HC", "CE"), c(23,69)))
rownames(group_df) = colnames(heatmap_data_up)

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
         main = "Heatmap of Top 25 Upregulated Genes",
         angle_col = 45,
         fontsize_col = 4,
         annotation_col = group_df,
         show_legend = TRUE)


heatmap_data_down = as.matrix(normalized_counts_df_down)
group_df = data.frame(Groups=rep(c("HC", "CE"), c(23,69)))
rownames(group_df) = colnames(heatmap_data_down)

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
         main = "Heatmap of Top 25 Downregulated Genes",
         angle_col = 45,
         fontsize_col = 4,
         annotation_col = group_df,
         show_legend = TRUE)





array_res$neg_logp_val = -log10(array_res$adj.P.Val)
array_res$Expression = ifelse(array_res$logFC >= 0.5 & array_res$adj.P.Val <= 0.05, "Upregulated",
                           ifelse(array_res$logFC <= -0.5 & 
                                    array_res$adj.P.Val <= 0.05, "Downregulated", "Not Significant"))

volcano_plot = ggplot(array_res, aes(x=logFC, y=neg_logp_val, colour = Expression)) +
  geom_point() +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed", color="darkgray") +  # Vertical fold-change cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgray") +  # Horizontal p-value cutoff line
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10(P-value)")

volcano_plot + ggtitle("Differential Gene Expression Analysis of\nControl vs Cardio Embolic Stroke")

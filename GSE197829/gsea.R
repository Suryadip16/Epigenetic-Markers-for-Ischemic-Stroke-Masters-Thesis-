# GSEA
# Credits: https://biostatsquid.com/fgsea-tutorial-gsea/


library(fgsea)
library(tidyverse)
library(RColorBrewer)

in_path = "C:/Users/dipak/Desktop/Riku/VM2_Project/edgeR_genes"
out_path = "C:/Users/dipak/Desktop/Riku/VM2_Project/GSEA"
setwd(in_path)

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list = function(pws){
  pws.l = list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] = rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
## Function: prepare_gmt --------------------------------------
prepare_gmt = function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file = gmt_files[1]
  #genes_in_data = df$gene_symbol
  
  # Read in gmt file
  gmt = gmtPathways(gmt_file)
  hidden = unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat = matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] = as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 = intersect(genes_in_data, hidden)
  mat = mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list = matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}
# Analysis ====================================================
## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df = read.csv("ce_esus_common_genes_final.csv")
## 2. Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes = df$X
# list.files(bg_path)

gmt_file = "C:/Users/dipak/Desktop/Riku/VM2_Project/GSEA/Background_genes/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
bg_genes = prepare_gmt(gmt_file, my_genes, savefile = FALSE)


# Ranking Gene List
rankings = sign(df$logFC)*(-log10(df$PValue)) # we will use the signed p values from spatial DGE as ranking
names(rankings) = df$X # genes as names#
rankings = sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Run GSEA

GSEAres = fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation


## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.05])


topPathwaysUp = GSEAres[ES > 0][head(order(padj), n = 10), pathway]
topPathwaysDown = GSEAres[ES < 0][head(order(padj), n = 10), pathway]
topPathways = c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()
# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)

GSEAres$leadingEdge = sapply(GSEAres$leadingEdge, function(x) paste(x, collapse = "; "))

write.csv(GSEAres, file = paste0(out_path, "/CE_vs_ESUS_gsea_res.csv"), row.names = FALSE)


significant_gs = subset(GSEAres, GSEAres$padj<0.05)


up_pathways = subset(significant_gs, significant_gs$ES>0)
down_pathways = subset(significant_gs, significant_gs$ES<0)

top20_up_pathways = up_pathways %>% 
  slice_max(order_by = ES, n = 20)

top20_up_pathways$pathway = gsub("^REACTOME_", "", top20_up_pathways$pathway)  # Remove "REACTOME_"
top20_up_pathways$pathway = gsub("_", " ", top20_up_pathways$pathway)  # Replace underscores with spaces



ggplot(top20_up_pathways, aes(x = ES, y = reorder(pathway, ES))) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +  # Adjust bar width
  labs(
    x = "Enrichment Score (ES)", 
    y = "Pathway", 
    title = "Upregulated Pathways: AT vs ESUS"
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
  geom_text(aes(label = round(ES, 2)), hjust = -0.2, size = 4)  # Add value labels to bars



top20_down_pathways = down_pathways %>% 
  slice_min(order_by = ES, n = 30)

top20_down_pathways$pathway = gsub("^REACTOME_", "", top20_down_pathways$pathway)  # Remove "REACTOME_"
top20_down_pathways$pathway = gsub("_", " ", top20_down_pathways$pathway)  # Replace underscores with spaces




ggplot(top20_down_pathways, aes(x = ES, y = reorder(pathway, ES))) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +  # Adjust bar width
  labs(
    x = "Enrichment Score (ES)", 
    y = "Pathway", 
    title = "Downregulated Pathways: CE vs ESUS"
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
  geom_text(aes(label = round(ES, 2)), hjust = -0.2, size = 4)  # Add value labels to bars


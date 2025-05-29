setwd("/home/ibab/HM_VM_Project/GSE197829/star_htseq_counts")

library(biomaRt)

# Read data
x = read.delim("GSE187829_htseq_counts_matrix.txt", row.names = "Gene_Name")
x$row.names = NULL
nrow(x)
head(x)

# Connect to Ensembl BioMart
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Retrieve gene biotypes and symbols in one query
annotations = getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(x),
  mart = ensembl
)

# Retain only protein-coding genes
annotations = annotations[annotations$gene_biotype == "protein_coding", ]

# Filter count matrix for protein-coding genes
counts_matrix_filtered = x[rownames(x) %in% annotations$ensembl_gene_id, ]

# Merge annotations (keeping only coding genes)
counts_matrix_filtered$ensembl_gene_id = rownames(counts_matrix_filtered)
counts_matrix_filtered = merge(counts_matrix_filtered, annotations[, c("ensembl_gene_id", "hgnc_symbol")], by = "ensembl_gene_id", all.x = TRUE)

# Remove rows without a gene symbol
counts_matrix_filtered = counts_matrix_filtered[!is.na(counts_matrix_filtered$hgnc_symbol) & counts_matrix_filtered$hgnc_symbol != "", ]

# Remove non-numeric columns before aggregation
counts_matrix_filtered_numeric = counts_matrix_filtered[, sapply(counts_matrix_filtered, is.numeric)]

# Add back the hgnc_symbol column for aggregation
counts_matrix_filtered_numeric$hgnc_symbol = counts_matrix_filtered$hgnc_symbol

# Aggregate (sum counts for duplicated gene symbols)
counts_matrix_final = aggregate(. ~ hgnc_symbol, data = counts_matrix_filtered_numeric, sum)

# Set gene symbols as row names
rownames(counts_matrix_final) = counts_matrix_final$hgnc_symbol

# Remove redundant hgnc_symbol column
counts_matrix_final$hgnc_symbol = NULL


# Save final matrix
write.csv(counts_matrix_final, file = "GSE187829_final_counts_matrix_coding.csv", row.names = TRUE)

nrow(counts_matrix_final)





library(GenomicFeatures)
library(biomaRt)

# Load the GTF file and create a TxDb object
txdb <- makeTxDbFromGFF("/home/ibab/HM_VM_Project/hg38_main/hg38.ensGene.gtf", format="gtf")

# Extract exon lengths per gene
exons_by_gene <- exonsBy(txdb, by="gene")
gene_lengths <- sum(width(reduce(exons_by_gene)))

# Convert to a data frame
gene_lengths_df <- data.frame(Ensembl_ID=names(gene_lengths), Length=as.numeric(gene_lengths))

# Connect to Ensembl BioMart
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Retrieve HGNC symbols for Ensembl IDs
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_lengths_df$Ensembl_ID,
  mart = ensembl
)

# Merge with original dataframe to retain gene lengths
gene_lengths_mapped <- gene_lengths_df %>%
  left_join(mapping, by = c("Ensembl_ID" = "ensembl_gene_id")) %>%
  dplyr::select(hgnc_symbol, Length) %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol))  # Remove missing symbols

# Write to CSV
write.csv(gene_lengths_mapped, "/home/ibab/HM_VM_Project/hg38_main/hg38_ensGene_gene_lengths.csv", row.names=FALSE)

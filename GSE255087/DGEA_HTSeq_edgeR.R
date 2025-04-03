# DGEA: HTSeq Counts

setwd("/home/ibab/HM_VM_Project/DGEA/ht_seq_counts")

library(edgeR)

x = read.delim("htseq_counts_matrix.txt", row.names = "Gene_Name")
x$row.names = NULL
nrow(x)
head(x)
x2 = x[!duplicated(x),]
nrow(x2)
x2 = x2[rowSums(x2) >= 10,  ]

write.csv(x2, file = "final_counts_matrix.csv", sep = ",")

plot(density(x2$at1), col = 'blue',
     xlim = c(0, 1200), ylim = c(0, 0.002), main = "Density Plot of Non Normalised Reads",
     xlab = "Value", ylab = "Density", lwd = 2)
lines(density(x2$at2), col = 'red', lwd = 2)
lines(density(x2$at3), col = 'black', lwd = 2)
lines(density(x2$at4), col = 'green', lwd = 2)
lines(density(x2$at5), col = 'yellow',lwd = 2 )
lines(density(x2$at6), col = 'darkred',lwd = 2 )
lines(density(x2$ce1), col = 'magenta',lwd = 2 )
lines(density(x2$ce2), col = 'cyan',lwd = 2 )
lines(density(x2$ce3), col = 'brown',lwd = 2 )
lines(density(x2$ce4), col = 'purple',lwd = 2 )
lines(density(x2$ce5), col = 'violet',lwd = 2 )
lines(density(x2$ce6), col = 'gray',lwd = 2 )
lines(density(x2$ce7), col = 'orange',lwd = 2 )
lines(density(x2$ce8), col = 'turquoise',lwd = 2 )
lines(density(x2$ce9), col = 'beige',lwd = 2 )
lines(density(x2$ce10), col = 'maroon',lwd = 2 )
lines(density(x2$esus1), col = 'darkblue',lwd = 2 )
lines(density(x2$esus2), col = 'gold',lwd = 2 )
lines(density(x2$esus3), col = 'darkgreen',lwd = 2 )
lines(density(x2$esus4), col = 'darkorange',lwd = 2 )
lines(density(x2$esus5), col = 'lightblue',lwd = 2 )

legend("topright", legend = c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6", 
                              "CE1", "CE2", "CE3", "CE4", "CE5", "CE6", "CE7", 
                              "CE8", "CE9", "CE10", "ESUS1", "ESUS2", "ESUS3", 
                              "ESUS4", "ESUS5"),
       col = c("blue", "red", "black", "green", "yellow", "darkred", "magenta", 
               "cyan", "brown", "purple", "violet", "gray", "orange", "turquoise", 
               "beige", "maroon", "darkblue", "gold", "darkgreen", "darkorange", 
               "lightblue"), lty = 1, lwd = 3, cex = 1.0, ncol = 3)

group = c(rep("AT", 6), rep("CE", 10), rep("ESUS", 5))
y = DGEList(counts = x2, group = group)

# Box plots:
boxplot(log2(cpm(y) + 1), main = "Before Normalization", col = rainbow(ncol(y)),
        outline = FALSE, ylim = c(0, 15))

# Calculate normalization factors to scale the library sizes. Uses TMM Normalization
y = calcNormFactors(y)

# Calculate TMM-normalized CPM values
normalized_counts <- cpm(y, normalized.lib.sizes = TRUE)

boxplot(log2(normalized_counts + 1), main = "After Normalization", 
        col = rainbow(ncol(y)), outline = F, ylim = c(0, 15))


# Plot density of normalized counts
plot(density(normalized_counts[, "at1"]), col = 'blue',
     xlim = c(0, 1200), ylim = c(0, 0.05),
     main = "Density Plot of TMM Normalized Reads",
     xlab = "Normalized Value", ylab = "Density", lwd = 2)

lines(density(normalized_counts[, "at2"]), col = 'red', lwd = 2)
lines(density(normalized_counts[, "at3"]), col = 'black', lwd = 2)
lines(density(normalized_counts[, "at4"]), col = 'green', lwd = 2)
lines(density(normalized_counts[, "at6"]), col = 'darkred', lwd = 2)
lines(density(normalized_counts[, "at5"]), col = 'yellow', lwd = 2)
lines(density(normalized_counts[, "ce1"]), col = 'magenta', lwd = 2)
lines(density(normalized_counts[, "ce2"]), col = 'cyan', lwd = 2)
lines(density(normalized_counts[, "ce3"]), col = 'brown', lwd = 2)
lines(density(normalized_counts[, "ce4"]), col = 'purple', lwd = 2)
lines(density(normalized_counts[, "ce5"]), col = 'violet', lwd = 2)
lines(density(normalized_counts[, "ce6"]), col = 'gray', lwd = 2)
lines(density(normalized_counts[, "ce7"]), col = 'orange', lwd = 2)
lines(density(normalized_counts[, "ce8"]), col = 'turquoise', lwd = 2)
lines(density(normalized_counts[, "ce9"]), col = 'beige', lwd = 2)
lines(density(normalized_counts[, "ce10"]), col = 'maroon', lwd = 2)
lines(density(normalized_counts[, "esus1"]), col = 'darkblue', lwd = 2)
lines(density(normalized_counts[, "esus2"]), col = 'gold', lwd = 2)
lines(density(normalized_counts[, "esus3"]), col = 'darkgreen', lwd = 2)
lines(density(normalized_counts[, "esus4"]), col = 'darkorange', lwd = 2)
lines(density(normalized_counts[, "esus5"]), col = 'lightblue', lwd = 2)

legend("topright", legend = c("AT1", "AT2", "AT3", "AT4", "AT5", "AT6", 
                              "CE1", "CE2", "CE3", "CE4", "CE5", "CE6", "CE7", 
                              "CE8", "CE9", "CE10", "ESUS1", "ESUS2", "ESUS3", 
                              "ESUS4", "ESUS5"),
       col = c("blue", "red", "black", "green", "yellow", "darkred", "magenta", 
               "cyan", "brown", "purple", "violet", "gray", "orange", "turquoise", 
               "beige", "maroon", "darkblue", "gold", "darkgreen", "darkorange", 
               "lightblue"), lty = 1, lwd = 3, cex = 1.0, ncol = 3)


# # Convert ENSMBL ids to gene ids
# 
# library(biomaRt)
# 
# # Load BioMart
# mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# ensmbl_ids = rownames(x)
# genes = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
#                filters = "ensembl_gene_id", 
#                values = ensmbl_ids, 
#                mart = mart)


# Building design matrix
design = model.matrix(~group)

# Calculating Tagwise dispersion for all the tags
y = estimateDisp(y, design)
# This function estimates the dispersion of the count data in the DGEList object
# y. Dispersion is a measure of variability that is specific to RNA-seq and 
# other count-based data. It quantifies the degree of over-dispersion (greater 
# variability than expected from a Poisson model) in the data.
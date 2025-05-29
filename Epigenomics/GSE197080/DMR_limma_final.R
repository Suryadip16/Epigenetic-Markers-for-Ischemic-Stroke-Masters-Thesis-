# DMR from scratch

# Load necessary libraries
library(minfi)
library(DMRcate)  # For DMR analysis
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

setwd("/home/ibab/HM_VM_Project/methylome/Methylome")


# Read the unnormalized β-values CSV
meth_int = read.csv(file = "GSE197080_meth_ints.csv", row.names = 'TargetID')
meth_int$X = NULL

unmeth_int = read.csv(file = "GSE197080_unmeth_ints.csv", row.names = 'TargetID')
unmeth_int$X = NULL


# # Step 1: Filter out probes on chrX and chrY (since they often show gender-specific patterns)
# annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# keep = !(row.names(meth_int) %in% annEPIC$Name[annEPIC$chr %in% c("chrX", "chrY")])
# meth_int = meth_int[keep, ]
# rm(keep)  # Clean-up after filtering
# 
# # Step 2: Filter out probes with SNP overlap (removing probes with high SNP frequencies)
# snp.probe = annEPIC[!is.na(annEPIC$Probe_rs), ]
# snp5.probe = snp.probe$Name[snp.probe$Probe_maf <= 0.05]  # MAF threshold <= 0.05
# no.snp.probe = annEPIC$Name[is.na(annEPIC$Probe_rs)]  # Probes without SNPs
# meth_int = meth_int[row.names(meth_int) %in% c(no.snp.probe, snp5.probe), ]
# rm(snp.probe, snp5.probe, no.snp.probe)  # Clean-up after filtering
# 
# # Step 1: Filter out probes on chrX and chrY (since they often show gender-specific patterns)
# annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# keep = !(row.names(unmeth_int) %in% annEPIC$Name[annEPIC$chr %in% c("chrX", "chrY")])
# unmeth_int = unmeth_int[keep, ]
# rm(keep)  # Clean-up after filtering
# 
# # Step 2: Filter out probes with SNP overlap (removing probes with high SNP frequencies)
# snp.probe = annEPIC[!is.na(annEPIC$Probe_rs), ]
# snp5.probe = snp.probe$Name[snp.probe$Probe_maf <= 0.05]  # MAF threshold <= 0.05
# no.snp.probe = annEPIC$Name[is.na(annEPIC$Probe_rs)]  # Probes without SNPs
# unmeth_int = unmeth_int[row.names(unmeth_int) %in% c(no.snp.probe, snp5.probe), ]
# rm(snp.probe, snp5.probe, no.snp.probe)  # Clean-up after filtering



# # Create MethylSet
# mset = MethylSet(Meth = as.matrix(meth_int), Unmeth = as.matrix(unmeth_int),
#                   annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
# 
# # Quantile normalization
# mset_normalized = preprocessQuantile(mset)

# Load required package
library(limma)

# Assume you already have these:
# methylated_matrix: rows = CpG probes, cols = samples (raw methylated signals)
# unmethylated_matrix: same structure (raw unmethylated signals)

# Step 1: Quantile normalization
methylated_norm = normalizeQuantiles(as.matrix(meth_int))
unmethylated_norm = normalizeQuantiles(as.matrix(unmeth_int))

# # Step 2: Compute Beta values
# # Add a pseudocount to avoid division by zero
# offset = 100  # common offset used in Illumina protocols
# beta_values = methylated_norm / (methylated_norm + unmethylated_norm + offset)
# 
# # Optional: set row and column names to match original
# rownames(beta_values) = rownames(meth_int)
# colnames(beta_values) = colnames(meth_int)
# 
# # Sample column names must match colnames(beta_values)
# group = factor(c(rep("Control", 3), rep("Case", 3)))
# design = model.matrix(~ group)







beta_values = read.csv(file = "GSE197080_beta_vals_calculated.csv", row.names = "TargetID")
beta_values$X = NULL

group = factor(c(rep("HC", 3), rep("IS", 3)))
design = model.matrix(~ group)


# Transform beta values to M-values (better for statistical testing)
M_values = log2(beta_values / (1 - beta_values))

M_values_laa_saa = M_values[ , 6:ncol(M_values)]

# Step 1: Load 450k annotation and filter probes on chrX and chrY
ann450k = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep = !(row.names(M_values) %in% ann450k$Name[ann450k$chr %in% c("chrX", "chrY")])
M_values = M_values[keep, ]
rm(keep)

# Step 2: Filter out probes with SNPs (as available in 450k annotation)
snp.probe = ann450k[!is.na(ann450k$Probe_rs), ]
snp5.probe = snp.probe$Name[snp.probe$Probe_maf <= 0.05]
no.snp.probe = ann450k$Name[is.na(ann450k$Probe_rs)]
M_values = M_values[row.names(M_values) %in% c(no.snp.probe, snp5.probe), ]
rm(snp.probe, snp5.probe, no.snp.probe)



# 1. Define your condition vector
conds = factor(c(rep("HC", 3), rep("IS", 3)))

# 2. Subset the numeric matrix (drop non-numeric or metadata columns if any)
#    Assuming M_val is already all numeric, and you just want all 6 columns:
mat = as.matrix(M_values)

# 3. Run PCA (center and scale are almost always recommended)
pca = prcomp(mat, center = TRUE, scale. = TRUE)

# 4. Extract scores (PC coordinates) into a data frame
scores = as.data.frame(pca$x[, 1:2])   # take PC1 & PC2
scores$Condition = conds
scores$Sample = rownames(scores)

# 5. Plot with ggplot2
library(ggplot2)
ggplot(scores, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, show.legend = FALSE) +
  xlab(paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  ggtitle("PCA of M_val by Condition")







# Fit linear model
fit = lmFit(M_values, design)
fit = eBayes(fit)

M_values_mat = as.matrix(M_values)


# 1. Annotate CpGs for DMR detection
myAnnotation = cpg.annotate(
  datatype = "array",
  object = M_values_mat,
  what = "M",
  analysis.type = "differential",
  design = design,
  coef = 2,
  arraytype = "EPICv1"
)

# 2. Identify DMRs
dmrcoutput = dmrcate(myAnnotation, lambda = 1000, C = 2)

# 3. Extract DMR ranges
dmr_ranges = extractRanges(dmrcoutput, genome = "hg19")

# 4. Save as CSV
dmr_df = as.data.frame(dmr_ranges)
write.csv(dmr_df, "GSE69138_DMRcate_DMRs.csv", row.names = FALSE)

# --------- DMRcate analysis ends here ---------












# Get DMP results
results = topTable(fit, coef = 2, number = Inf, adjust = "fdr")

# Compute delta beta
beta_case = beta_values[, group == "SAA"]
beta_control = beta_values[, group == "LAA"]
delta_beta = rowMeans(beta_case) - rowMeans(beta_control)

results$delta_beta = delta_beta[rownames(results)]
results = na.omit(results)

# Apply FDR filter
sig_dmr = results[results$P.Value <= 0.05, ]
sig_dmr = na.omit(sig_dmr)

# Classify based on delta beta
sig_dmr$direction = ifelse(sig_dmr$delta_beta >= 0.17, "Hypermethylated",
                            ifelse(sig_dmr$delta_beta <= -0.17, "Hypomethylated", "Not significant"))
sig_dmr = sig_dmr[sig_dmr$direction != "Not significant", ]

# Annotate significant DMRs
probeIDs = rownames(sig_dmr)
selAnno.hg19 = ann450k[ann450k$Name %in% probeIDs, ]
annotated.hg19 = merge(
  selAnno.hg19,
  sig_dmr,
  by.x = "Name", by.y = "row.names",
  all.x = TRUE
)

# Export
write.csv(annotated.hg19, "GSE69138_annotated_450k_DMRs_norm.csv", row.names = FALSE)

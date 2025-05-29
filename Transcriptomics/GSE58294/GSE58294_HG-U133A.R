setwd("/home/chg/Suryadip Microarray/GSE58294")
#---------------------------------
#	SET PARAMETERS
#---------------------------------
# Affymetrix chip used
Affy_chip = "hgu133a2"
#BiocManager::install("hgu133a2.db")
#BiocManager::install("pd.hg.u133a.2")
#BiocManager::install("RColorBrewer")
# GEO number for experiment
gse_current = "GSE58294"
# covdesc file name

#covdesc_file = paste("covdesc",gse_current,sep="_") 
covdesc = read.delim("covdesc.txt", sep="\t", header=T)

# input and output file directories.
infile_dir = "data_files"

outfile_dir = "output_files"

test_between = c("control","ischemic")
## Decide which files are control and which are treatment.
control_cols = which(covdesc$treatment == test_between[1])

treatment_cols = which(covdesc$treatment == test_between[2])
#---------------------------------
#	LOAD PACKAGES
#---------------------------------
library(oligo) # normalization
library(limma) # DE analysis
library(multtest) # estimating False Discovery Rate (FDR)	

# load chip library for annotation
if(Affy_chip == "hgu133a2")
{
  library(hgu133a2.db)
  library(pd.hg.u133a.2)
}

if(Affy_chip == "hgu95av2")
{
  #library(hgu95av2cdf) # not required for oligo
  library(hgu95av2.db)
  library(pd.hg.u95av2)
}

if(Affy_chip == "hgu133plus2")
{
  library(hgu133plus2.db)
  #library(hgu133plus2cdf) ## not required for oligo
  library(pd.hg.u133.plus.2)
  
}

if(Affy_chip == "hgu95a")
{
  #library(hgu95acdf) # not required for oligo
  library(hgu95a.db)
  library(pd.hg.u95a)
}

if(Affy_chip == "hgu133a")
{
  library(hgu133a.db)
  # library(hgu133acdf) not required for oligo
  library(pd.hg.u133a)
}

if(Affy_chip == "hgu133b")
{
  library(hgu133b.db)
  #library(hgu133bcdf) # not required for oligo
  library(pd.hg.u133b)
}

# for hugene chips, oligo will automatically detect and load the required libraries. however, when not installed previously in the session, this show an error. 
if(Affy_chip == "hugene11st")
{
  library(pd.hugene.1.1.st.v1)
  library(hugene11stprobeset.db)
  library(hugene11sttranscriptcluster.db) # hugene platforms can be analysed at transcript level or probe level
}

#---------------------------------
#	NORMALIZE
#---------------------------------
# read the cel files
celFiles = list.celfiles("data_files",full.names=T)
rawData = read.celfiles(celFiles)

# rawData is now a GeneFeatureSet object. 
# this object does not contain information on phenotypes, which has to be loaded separately from the covdesc file
# pData(rawData) retrieves phenodata from the GenefeatureSet object as a dataframe.
pData(rawData)$gsm = rownames(pData(rawData))
pData(rawData)$source = covdesc$treatment


# Extract raw expression values (before normalization)
raw_exprs = log2(pm(rawData))  # Log2 transformation of probe intensities

# Perform RMA normalization
normalizedData = rma(rawData)
#png("volcano_plot.png", width = 1200, height = 800, res = 150, type = "cairo")
X11()
# Extract normalized expression values
norm_exprs = exprs(normalizedData)
# Define group labels dynamically
group_labels = covdesc$treatment

# Save plots to PNG file
#png("boxplot_normalization_comparison.png", width = 1200, height = 600, res = 150, type = "cairo")

X11()
# Set 1-row, 2-column layout for side-by-side plots
par(mfrow = c(1, 2))  

# Boxplot: Before Normalization
boxplot(raw_exprs, 
        main = "Before Normalization",
        col = brewer.pal(8, "Paired")[1:length(unique(group_labels))], 
        las = 2,
        names = group_labels,  
        ylab = "Log2 Intensity",
        outline = FALSE)

# Boxplot: After Normalization
boxplot(norm_exprs, 
        main = "After Normalization",
        col = brewer.pal(8, "Paired")[1:length(unique(group_labels))],  
        las = 2,  
        names = group_labels,  
        ylab = "Log2 Expression",
        outline = FALSE)

# Ensure the file is correctly saved
dev.off()

# Double-check that the file is a PNG
if (file.exists("boxplot_normalization_comparison.png")) {
  print("PNG file successfully saved.")
} else {
  print("Error: PNG file was not created.")
}

# generate expression values matrix
expval = exprs(normalizedData)
prbnames = rownames(expval)
clnames = colnames(expval)

## get rid of ".cel" extension names from clnames and from colnames(expval)
nme = strsplit(colnames(expval), "\\.")
for(i in 1:ncol(expval)){
  clnames[i] = strsplit(colnames(expval), "\\.")[[i]][1]
}
colnames(expval) = clnames

X11()
# Principal component analysis of the normalized data
data.PC = prcomp(t(expval),scale=T)

color = c("olivedrab","goldenrod")

plot(data.PC$x,col=color,type="p", pch=19)

## DE analysis using limma


# create design matrix for anova: group samples into lean and obese

groups = pData(rawData)$source
f = factor(groups,levels=test_between)
design = model.matrix(~ 0 + f)
colnames(design) = test_between


# lmFit() calculates mean expression levels in each group
data.fit = lmFit(expval,design)

# data.fit is an MArrayLM onject defined in limma
# first two columns contain mean log expression values of lean and obese

# use makeContrasts() to create a contrast matrix that show whic is the control and which is the test
contrast.matrix = makeContrasts(control-ischemic,levels=design)
# contrasts.fit() fits the contrast.matrix to the data.fit object with the log mean expression values
data.fit.con = contrasts.fit(data.fit,contrast.matrix)

# perform moderated t-test using limma (eBayes())
data.fit.eb = eBayes(data.fit.con)

#generate Volcano plot
#volcanoplot(data.fit.eb,coef=1,col="cornflowerblue")

# Generate a table of results from limma
tab = topTable(data.fit.eb, coef=1, number=Inf, adjust="BH")

#write.csv(tab)

# Add a new column to categorize the genes based on logFC and adjusted p-value
tab$Category = "Insignificant"
tab$Category[tab$logFC > 1 & tab$adj.P.Val < 0.05] = "Overexpressed"
tab$Category[tab$logFC < -1 & tab$adj.P.Val < 0.05] = "Underexpressed"

# Define colors based on the category
tab$Color = "grey" # Default color for insignificant genes
tab$Color[tab$Category == "Overexpressed"] = "darkorange"
tab$Color[tab$Category == "Underexpressed"] = "darkgreen"


X11()
# Create the volcano plot
plot(
  tab$logFC, 
  -log10(tab$adj.P.Val),
  main = "Volcano plot for GSE58294",
  col = tab$Color, 
  pch = 19, 
  xlab = "Log2 Fold Change", 
  ylab = "-Log10 Adjusted P-Value", 
  cex = 0.6
)
# Add legend
#legend("topright", legend=c("Overexpressed", "Underexpressed", "Insignificant"),
#      col=c("darkorange", "darkgreen", "grey"), pch=20, bty="n")


# Add horizontal and vertical lines for thresholds
abline(h = -log10(0.05), col = "black", lty = 2)
abline(v = c(-1, 1), col = "black", lty = 2)

#multiple hypothesis correction
tab = topTable(data.fit.eb,coef=1, number = nrow(normalizedData),adjust.method="BH",sort.by = "P") 

## In the above table, probes are sorted by p-value.
##  We will restore the order to same order of probes in our data so that we can combine with other data.

limma_output = tab[order(match(rownames(tab), rownames(fData(normalizedData)))), , drop = FALSE]

# The above data frame "limma_output" has same probe id order as normalizedData.

## =---------------Fold change computation--------------------

# Arithmetic Mean and Fold Change
# compute Fold change = obese/lean

control_samples = length(control_cols)
treatment_samples = length(treatment_cols)



# creating empty vectors of length equals number of rows in expval.
control_data =  rep(0, nrow(expval))
treatment_data   =  rep(0,nrow(expval))

# Summing all columns of lean
for(i in control_cols)
{
  control_data = control_data + 2^expval[,i]
}

# summing all columns of obese
for(i in treatment_cols)
{
  treatment_data = treatment_data + 2^expval[,i]
}
control_linear_mean = (control_data/control_samples)
treatment_linear_mean = (treatment_data/treatment_samples)

# Fold change is oa column divided by ra column. We get fold change for all genes.
FC_linear = as.vector(treatment_linear_mean/control_linear_mean)


FC_log2 = log2(FC_linear)

FC_linear = round(FC_linear, 4)
FC_log2 = round(FC_log2,4)

X11()

## Histogram of FC_log2
hist(FC_log2, breaks=60, xlim=c(-3,3))

# Geometric Mean and Fold Change

# compute Fold change = obese/lean

# creating empty vectors of length equals number of rows in expval.
obese_data =  rep(0, nrow(expval))
lean_data   =  rep(0,nrow(expval))

## geometric mean = exp(mean(log(x))). But here, expval is in log2 scale. So, geometric mean = 2^(mean(log(x)))

## control.g and treatment.g are mean(log(x))
control.g = apply(expval[, control_cols], 1, mean) # passing 1 means the mean will be calculated along rows)
treatment.g = apply(expval[, treatment_cols], 1, mean) # passing 1 means the mean will be calculated along rows)

## geometric mean is 2^( mean(log(x)) )
control.g.mean = 2^control.g
treatment.g.mean = 2^treatment.g

# Fold change of geometric mean is treatment column divided by control column. We get fold change for all genes.
FC.g_linear = as.vector(treatment.g.mean/control.g.mean)

FC.g_log2 = log2(FC.g_linear)

FC.g_linear = round(FC.g_linear, 4)
FC.g_log2 = round(FC.g_log2,4)

X11()
plot(FC_linear, FC.g_linear, xlim=c(0,10), ylim=c(0,10))

X11()

#Generate Volcano plot for uncorrected p-value
#plot( FC_log2, -log10(pvals.uncorrected.genes.present), col="darkblue", cex = 0.3, xlab="log2(Fold change)", ylab="log10(p-value)", 
#       main="Volcano plot for FDR p-values" )

#X11()

## Generate Volcano plot for fdr p-value
#plot( FC_log2, log10(pvals.fdr.genes.present), col="red", cex = 0.3, xlab="log2(Fold change)", ylab="log10(p-value)", 
#       main="Volcano plot for uncorrected p-values" )


#BiocManager::install("hgu219.db", force=TRUE)
#BiocManager::install("AnnotationDbi", force=TRUE)
#---------------------------------
#	ANNOTATE
#---------------------------------

# Load the appropriate annotation package for HGU129
library(AnnotationDbi)
library(hgu133a2.db)

annot = data.frame(
  ENTREZID = sapply(as.list(hgu133a2ENTREZID), paste, collapse = ", "),
  ACCNUM = sapply(as.list(hgu133a2ACCNUM), paste, collapse = ", "),
  SYMBOL = sapply(as.list(hgu133a2SYMBOL), paste, collapse = ", "),
  DESC = sapply(as.list(hgu133a2GENENAME), paste, collapse = ", "),
  CHROMOSOME = sapply(as.list(hgu133a2CHR), paste, collapse = ", ")
)

# Match limma output with normalized data and annotations
limma_output = tab[order(match(rownames(tab), rownames(fData(normalizedData)))), , drop = FALSE]

fData(normalizedData) = annot[match(rownames(fData(normalizedData)), rownames(annot)), ]
head(fData(normalizedData))

# Combine all data into a single data frame
df = data.frame(
  ProbeNames = prbnames,
  fData(normalizedData),
  limma_output,
  expval
)

library(dplyr)

# Filter out mismatched probes, probes without a name, and apply logFC and adjusted p-value criteria
filtered_df = df %>%
  # Remove rows where SYMBOL is NA or empty (probes without gene names)
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  # Remove rows with mismatched probes
  filter(!is.na(ProbeNames)) %>%
  # Filter rows based on logFC > 2 or logFC < -2 and adj.P.Val < 0.05
  filter( P.Value < 0.05)

filtered_df=na.omit(filtered_df)

# Define output file name
outfile_name = paste("output_files/", gse_current, "_filtered_limma_output.csv", sep = "")

# Write the filtered table to a CSV file
write.csv(filtered_df, file = outfile_name, row.names = FALSE)

#---------------------------------
#	END
#---------------------------------

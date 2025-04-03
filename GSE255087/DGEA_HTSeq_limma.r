setwd("/home/ibab/HM_VM_Project/harsha_di/limma3")
data <- read.csv("at_esus_final_counts_matrix.csv")
head(data)
rownames(data)<- data$X
data$X <- NULL
d0 <- DGEList(data)
#######Calculate normalization factors
d0 <- calcNormFactors(d0)
d0
#######Note: calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream.
#Filter low-expressed genes

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
##########
snames <- colnames(data) # Sample names
snames
group<-as.factor(c(rep("AT", 3), rep("ESUS", 4)))
head(group)
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

design = model.matrix(~0 + group)
design
colnames(design) = gsub("group", "", colnames(design))
group
contr.matrix = makeContrasts(ATvsESUS = AT-ESUS,
                             levels = colnames(design))
contr.matrix
#######Estimate contrast for each gene
# head(fit)
fit <- lmFit(y, design)
head(coef(fit))
tmp <- contrasts.fit(fit, contr.matrix)
head(tmp)
######Empirical Bayes smoothing of standard errors 
tmp <- eBayes(tmp)
head(tmp)
limma_res <- topTable(tmp, sort.by = "P", n = Inf)
nrow(limma_res)
head(limma_res)
pval <- limma_res[limma_res$P.Value<=0.05,]
nrow(pval)
write.csv(pval,"at_vs_esus_Cutoff_10_DGE_pval_0.05.csv")
upreg <- pval[pval$logFC>=1.5,]
nrow(upreg)
head(upreg)
#####Upreg genes_271
write.csv(upreg,"at_vs_esus_Cutoff_10_Upregulated_Genes_1.5.csv")
downreg <- pval[pval$logFC<=-1.5,]
nrow(downreg)
write.csv(downreg,"at_vs_esus_Cutoff_10_Downreg_genes_1.5.csv")
########Downreg 1210
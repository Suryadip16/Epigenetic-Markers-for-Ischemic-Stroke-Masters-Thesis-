library(dplyr)
setwd("/home/ibab/HM_VM_Project/methylome/Methylome")

meth = read.csv(file = "GSE197080_matrix_signal.txt", sep = "\t")
# Create an empty list to store the results

# List of sample names (column names that include 'Signal_A' and 'Signal_B')
sample_names = c("ZXJ", "YXG", "ZSS", "RWH", "HM", "ZM")

# Initialize result dataframe with TargetID
beta_df = data.frame(TargetID = meth$TargetID)

# Calculate beta values and add each as a column
for (sample in sample_names) {
  signal_A = meth[[paste0(sample, ".Signal_A")]]
  signal_B = meth[[paste0(sample, ".Signal_B")]]
  beta = signal_A / (signal_A + signal_B + 100)  # Adding α = 100
  beta_df[[sample]] = beta
}

colnames(beta_df)[2:4] = c("hc1", "hc2", "hc3")
colnames(beta_df)[5:7] = c("is1", "is2", "is3")
write.csv(beta_df, file = "GSE197080_beta_vals_calculated.csv")

# Process DMR res

setwd("/home/ibab/HM_VM_Project/methylome/Methylome")
laa_vs_saa_dmr_res = read.csv(file = "annotated_EPIC_DMRs_norm_1.csv", row.names = 1)
laa_vs_saa_dmr_res_cleaned = laa_vs_saa_dmr_res[!(is.na(laa_vs_saa_dmr_res$UCSC_RefGene_Name) | laa_vs_saa_dmr_res$UCSC_RefGene_Name == ""), ]
laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name
write.csv(laa_vs_saa_dmr_res_cleaned, file = "GSE197080_annotated_clean_850k_DMRs_norm.csv", row.names = T)
laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name <- sub(";.*", "", laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name)
laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name


dmrs = read.csv(file = "GSE69138_dmrs_res_refseq_symb.csv", row.names = 1)
dmrs = dmrs[!(is.na(dmrs$hgnc_symb) | dmrs$hgnc_symb == ""), ]

# intersect DMR and DEG res

laa_vs_ctrl_degs = read.csv("/home/ibab/HM_VM_Project/GSE197829/star_htseq_counts/laa_vs_ctrl_all_genes.csv", row.names = 1)

laa_vs_ctrl_degs_dmr_intersect = subset(laa_vs_ctrl_degs, rownames(laa_vs_ctrl_degs) %in% laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name)
laa_vs_ctrl_degs_dmr_intersect$UCSC_RefGene_Name = rownames(laa_vs_ctrl_degs_dmr_intersect)

laa_vs_ctrl_dmr_deg = inner_join(laa_vs_ctrl_degs_dmr_intersect, laa_vs_saa_dmr_res_cleaned, by = "UCSC_RefGene_Name")
laa_vs_ctrl_dmr_deg = laa_vs_ctrl_dmr_deg[!duplicated(laa_vs_ctrl_dmr_deg$UCSC_RefGene_Name), ]



rownames(laa_vs_ctrl_dmr_deg) = laa_vs_ctrl_dmr_deg$UCSC_RefGene_Name
laa_vs_ctrl_dmr_deg$UCSC_RefGene_Name = NULL

write.csv(laa_vs_ctrl_dmr_deg, file = "laa_vs_ctrl_dmr_deg_res.csv", row.names = T)

laa_vs_ctrl_dmr_deg_hyper = subset(laa_vs_ctrl_dmr_deg, laa_vs_ctrl_dmr_deg$log2FoldChange < 0 & laa_vs_ctrl_dmr_deg$direction == "Hypermethylated")
rownames(laa_vs_ctrl_dmr_deg_hyper)
write.csv(laa_vs_ctrl_dmr_deg_hyper, file = "laa_vs_ctrl_dmr_deg_hyper_res.csv", row.names = T)

laa_vs_ctrl_dmr_deg_hypo = subset(laa_vs_ctrl_dmr_deg, laa_vs_ctrl_dmr_deg$log2FoldChange > 0 & laa_vs_ctrl_dmr_deg$direction == "Hypomethylated")
rownames(laa_vs_ctrl_dmr_deg_hypo)
write.csv(laa_vs_ctrl_dmr_deg_hypo, file = "laa_vs_ctrl_dmr_deg_hypo_res.csv", row.names = T)







# intersect DMR and DEG res

saa_vs_ctrl_degs = read.csv("/home/ibab/HM_VM_Project/GSE197829/star_htseq_counts/saa_vs_ctrl_all_genes.csv", row.names = 1)

saa_vs_ctrl_degs_dmr_intersect = subset(saa_vs_ctrl_degs, rownames(saa_vs_ctrl_degs) %in% laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name)
saa_vs_ctrl_degs_dmr_intersect$UCSC_RefGene_Name = rownames(saa_vs_ctrl_degs_dmr_intersect)

saa_vs_ctrl_dmr_deg = inner_join(saa_vs_ctrl_degs_dmr_intersect, laa_vs_saa_dmr_res_cleaned, by = "UCSC_RefGene_Name")
saa_vs_ctrl_dmr_deg = saa_vs_ctrl_dmr_deg[!duplicated(saa_vs_ctrl_dmr_deg$UCSC_RefGene_Name), ]



rownames(saa_vs_ctrl_dmr_deg) = saa_vs_ctrl_dmr_deg$UCSC_RefGene_Name
saa_vs_ctrl_dmr_deg$UCSC_RefGene_Name = NULL

write.csv(saa_vs_ctrl_dmr_deg, file = "saa_vs_ctrl_dmr_deg_res.csv", row.names = T)

saa_vs_ctrl_dmr_deg_hyper = subset(saa_vs_ctrl_dmr_deg, saa_vs_ctrl_dmr_deg$log2FoldChange < 0 & saa_vs_ctrl_dmr_deg$direction == "Hypermethylated")
rownames(saa_vs_ctrl_dmr_deg_hyper)
write.csv(saa_vs_ctrl_dmr_deg_hyper, file = "saa_vs_ctrl_dmr_deg_hyper_res.csv", row.names = T)

saa_vs_ctrl_dmr_deg_hypo = subset(saa_vs_ctrl_dmr_deg, saa_vs_ctrl_dmr_deg$log2FoldChange > 0 & saa_vs_ctrl_dmr_deg$direction == "Hypomethylated")
rownames(saa_vs_ctrl_dmr_deg_hypo)
write.csv(saa_vs_ctrl_dmr_deg_hypo, file = "saa_vs_ctrl_dmr_deg_hypo_res.csv", row.names = T)








# intersect DMR and DEG res

ce_vs_ctrl_degs = read.csv("/home/ibab/HM_VM_Project/array/GSE58294_(copy)/output_files/ce_vs_ctrl_all_genes_pval_cut.csv", row.names = 1)

ce_vs_ctrl_degs_dmr_intersect = subset(ce_vs_ctrl_degs, rownames(ce_vs_ctrl_degs) %in% laa_vs_saa_dmr_res_cleaned$UCSC_RefGene_Name)
ce_vs_ctrl_degs_dmr_intersect$UCSC_RefGene_Name = rownames(ce_vs_ctrl_degs_dmr_intersect)

ce_vs_ctrl_dmr_deg = inner_join(ce_vs_ctrl_degs_dmr_intersect, laa_vs_saa_dmr_res_cleaned, by = "UCSC_RefGene_Name")
ce_vs_ctrl_dmr_deg = ce_vs_ctrl_dmr_deg[!duplicated(ce_vs_ctrl_dmr_deg$UCSC_RefGene_Name), ]



rownames(ce_vs_ctrl_dmr_deg) = ce_vs_ctrl_dmr_deg$UCSC_RefGene_Name
ce_vs_ctrl_dmr_deg$UCSC_RefGene_Name = NULL

write.csv(ce_vs_ctrl_dmr_deg, file = "ce_vs_ctrl_dmr_deg_res.csv", row.names = T)

ce_vs_ctrl_dmr_deg_hyper = subset(ce_vs_ctrl_dmr_deg, ce_vs_ctrl_dmr_deg$logFC.x > 0 & ce_vs_ctrl_dmr_deg$direction == "Hypermethylated")
rownames(ce_vs_ctrl_dmr_deg_hyper)
write.csv(ce_vs_ctrl_dmr_deg_hyper, file = "ce_vs_ctrl_dmr_deg_hyper_res.csv", row.names = T)

ce_vs_ctrl_dmr_deg_hypo = subset(ce_vs_ctrl_dmr_deg, ce_vs_ctrl_dmr_deg$logFC.x < 0 & ce_vs_ctrl_dmr_deg$direction == "Hypomethylated")
rownames(ce_vs_ctrl_dmr_deg_hypo)
write.csv(ce_vs_ctrl_dmr_deg_hypo, file = "ce_vs_ctrl_dmr_deg_hypo_res.csv", row.names = T)

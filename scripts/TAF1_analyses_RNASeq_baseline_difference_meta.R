library(ggplot2)
library(stringr)
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/lib/R/goplot.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/glmm.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/plot.R")


setwd("~/projects/xdp_aso/nsc/output/")

key_name = "ENSG00000147133"
res_xv = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPvanilla_vs_CON.txt", header = T, sep = "\t", check.names = F)
rownames(res_xv) = as.vector(res_xv$name)
res_xn = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPnoEdit_vs_CON.txt", header = T, sep = "\t", check.names = F)
rownames(res_xn) = as.vector(res_xn$name)
res_xc = read.table("res_FDR_LFC_glm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)
rownames(res_xc) = as.vector(res_xc$name)
res_cp0 = read.table("cellpaper_FDR_LFC_glmm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)
rownames(res_cp0) = as.vector(res_cp0$name)
res_cp1 = read.table("res_FDR_LFC_glmm_Exp_Genotype_NSC_XDP_vs_CON_cellpaper.txt", header = T, sep = "\t", check.names = F)
rownames(res_cp1) = as.vector(res_cp1$name)
res_cp2 = read.table("res_FDR_LFC_glmm_Exp_Genotype_XDP_vs_CON_mergedData.txt", header = T, sep = "\t", check.names = F)
rownames(res_cp2) = as.vector(res_cp2$name)

############################################
#Meta analysis with original Cell data
############################################
query = "cellpaper"
tt1 = meta_p(list(res_xc, res_cp0), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")
tt2 = meta_p(list(res_xn, res_cp0), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")
tt3 = meta_p(list(res_xv, res_cp0), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")

write.table(tt1, file = paste0("res_metaAnalysis_DEG_XDP_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tt2, file = paste0("res_metaAnalysis_DEG_XDPnoEdit_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tt3, file = paste0("res_metaAnalysis_DEG_XDPvanilla_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)

t1_up = rownames(tt1[which(tt1$meta_fdr<0.1 & tt1$direction == "up"), ])
t2_up = rownames(tt2[which(tt2$meta_fdr<0.1 & tt2$direction == "up"), ])
t3_up = rownames(tt3[which(tt3$meta_fdr<0.1 & tt3$direction == "up"), ])
t1_dn = rownames(tt1[which(tt1$meta_fdr<0.1 & tt1$direction == "down"), ])
t2_dn = rownames(tt2[which(tt2$meta_fdr<0.1 & tt2$direction == "down"), ])
t3_dn = rownames(tt3[which(tt3$meta_fdr<0.1 & tt3$direction == "down"), ])

t0_up = Reduce(intersect, list(t1_up, t2_up, t3_up))
t0_dn = Reduce(intersect, list(t1_dn, t2_dn, t3_dn))
t0 = unique(c(t0_up, t0_dn))

write.table(tt1[t0, ], file = paste0("res_metaAnalysis_with_", query,"_3overlaps.txt"), col.names = T, row.names = T, sep = "\t", quote = F)

ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")
deg_common_gs = ref[t0, 1]
shared_x = Reduce(intersect, list(rownames(res_xc), rownames(res_xn), rownames(res_xv), rownames(res_cp0)))
bg_gs = ref[shared_x, 1]
go_deg_common = gotest(deg_common_gs, bg_gs)
goplot(go_deg_common[1:5,], cutoff = 0.1, method = "p.value") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))
ggsave("pdf/GO_DEG_common.jpg")


############################################
#Meta analysis with reprocessed Cell data
############################################
query = "cellpaperReprocessd"
tt1 = meta_p(list(res_xc, res_cp1), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")
tt2 = meta_p(list(res_xn, res_cp1), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")
tt3 = meta_p(list(res_xv, res_cp1), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")

write.table(tt1, file = paste0("res_metaAnalysis_DEG_XDP_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tt2, file = paste0("res_metaAnalysis_DEG_XDPnoEdit_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tt3, file = paste0("res_metaAnalysis_DEG_XDPvanilla_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)

t1_up = rownames(tt1[which(tt1$meta_fdr<0.1 & tt1$direction == "up"), ])
t2_up = rownames(tt2[which(tt2$meta_fdr<0.1 & tt2$direction == "up"), ])
t3_up = rownames(tt3[which(tt3$meta_fdr<0.1 & tt3$direction == "up"), ])
t1_dn = rownames(tt1[which(tt1$meta_fdr<0.1 & tt1$direction == "down"), ])
t2_dn = rownames(tt2[which(tt2$meta_fdr<0.1 & tt2$direction == "down"), ])
t3_dn = rownames(tt3[which(tt3$meta_fdr<0.1 & tt3$direction == "down"), ])

t0_up = Reduce(intersect, list(t1_up, t2_up, t3_up))
t0_dn = Reduce(intersect, list(t1_dn, t2_dn, t3_dn))
t0 = unique(c(t0_up, t0_dn))

write.table(tt1[t0, ], file = paste0("res_metaAnalysis_with_", query,"_3overlaps_sigOnly.txt"), col.names = T, row.names = T, sep = "\t", quote = F)

ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")
deg_common_gs = ref[t0, 1]
bg_gs = ref[rownames(tt1), 1]
go_deg_common = gotest(deg_common_gs, bg_gs)
goplot(go_deg_common[1:5,], cutoff = 0.1, method = "p.value") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))
ggsave("pdf/GO_DEG_common.jpg")


######################################################
#Meta analysis with merged (Cell data + Current data)
######################################################
query = "cellpaperMergedData"
tt1 = meta_p(list(res_xc, res_cp2), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")
tt2 = meta_p(list(res_xn, res_cp2), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")
tt3 = meta_p(list(res_xv, res_cp2), col_feature = 4, col_pval = 3, col_fc = 2, padj_method = "bonferroni")

write.table(tt1, file = paste0("res_metaAnalysis_DEG_XDP_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tt2, file = paste0("res_metaAnalysis_DEG_XDPnoEdit_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)
write.table(tt3, file = paste0("res_metaAnalysis_DEG_XDPvanilla_with_", query,".txt"), col.names = T, row.names = T, sep = "\t", quote = F)

t1_up = rownames(tt1[which(tt1$meta_fdr<0.1 & tt1$direction == "up"), ])
t2_up = rownames(tt2[which(tt2$meta_fdr<0.1 & tt2$direction == "up"), ])
t3_up = rownames(tt3[which(tt3$meta_fdr<0.1 & tt3$direction == "up"), ])
t1_dn = rownames(tt1[which(tt1$meta_fdr<0.1 & tt1$direction == "down"), ])
t2_dn = rownames(tt2[which(tt2$meta_fdr<0.1 & tt2$direction == "down"), ])
t3_dn = rownames(tt3[which(tt3$meta_fdr<0.1 & tt3$direction == "down"), ])

t0_up = Reduce(intersect, list(t1_up, t2_up, t3_up))
t0_dn = Reduce(intersect, list(t1_dn, t2_dn, t3_dn))
t0 = unique(c(t0_up, t0_dn))

write.table(tt1[t0, ], file = paste0("res_metaAnalysis_with_", query,"_3overlaps_sigOnly.txt"), col.names = T, row.names = T, sep = "\t", quote = F)

ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")
deg_common_gs = ref[t0, 1]
bg_gs = ref[rownames(tt1), 1]
go_deg_common = gotest(deg_common_gs, bg_gs)
goplot(go_deg_common[1:5,], cutoff = 0.1, method = "p.value") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))
ggsave("pdf/GO_DEG_common.jpg")




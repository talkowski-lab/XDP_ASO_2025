library(ggplot2)
library(stringr)
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/lib/R/goplot.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/glmm.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/go.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/plot.R")


setwd("~/projects/xdp_aso/nsc/output/")


####################################################################################################################################################################################
#Untreated Overlaps
####################################################################################################################################################################################
####################################################################################################################################################################################
#This study
#Each DEG list
#############################################
#Baseline: XDP vs CON 
key_name = "ENSG00000147133"
grp = read.table(paste0("model_params_colData_glm_Exp_Genotype1_XDP_vs_CON.txt"), header = T, row.names = 1, sep = "\t", check.names = F)
grp$MIN = factor(grp$MIN)
grp[which(grp$Genotype1 == "XDP_noEdit"), "Genotype1"] = "XDP_unedited"
grp[which(grp$Genotype1 == "XDP_vanilla"), "Genotype1"] = "XDP_unexposed" 

mat1 = read.table("model_params_deSV_glm_Exp_Genotype1_XDP_vs_CON.txt", header = T, row.names = 1, sep = "\t")
rnaPCA(mat1, grp, c("Genotype1"), ntop = 500)
ggsave("pdf/PCA_rnaseq_Exp_Genotype1_CON_vs_XDP_corrected_new.jpg")

new1 = cbind(Exp = unlist(mat1[key_name, ]), grp)
ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_rnaseq_Exp_Genotype1_XDP_vs_CON_corrected_new.jpg")

mat2 = read.table("model_params_deSV_glm_Exp_Genotype_CON_vs_XDP.txt", header = T, row.names = 1, sep = "\t")
rnaPCA(mat2, grp, c("Genotype1"), ntop = 500)
ggsave("pdf/PCA_rnaseq_Exp_Genotype_CON_vs_XDP_corrected_new.jpg")
rnaPCA(mat2, grp, c("Genotype"), ntop = 500)
ggsave("pdf/PCA_rnaseq_Exp_Genotype_CON_vs_XDP_corrected_new2.jpg")

new2 = cbind(Exp = unlist(mat2[key_name, ]), grp)
ggplot(new2, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_rnaseq_Exp_Genotype_XDP_vs_CON_corrected_new.jpg")


res_xv = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPvanilla_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_xn = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPnoEdit_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_xc = read.table("res_FDR_LFC_glm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)

rownames(res_xv) = as.vector(res_xv$name)
rownames(res_xn) = as.vector(res_xn$name)
rownames(res_xc) = as.vector(res_xc$name)

deg_xv = rownames(res_xv)[which(res_xv$FDR.Genotype1XDP_vanilla < 0.1)]
deg_xv_up = rownames(res_xv)[which((res_xv$FDR.Genotype1XDP_vanilla < 0.1) & (res_xv$LFC.Genotype1XDP_vanilla > 0))]
deg_xv_dn = rownames(res_xv)[which((res_xv$FDR.Genotype1XDP_vanilla < 0.1) & (res_xv$LFC.Genotype1XDP_vanilla < 0))]

deg_xv_p = rownames(res_xv)[which(res_xv$pval.Genotype1XDP_vanilla < 0.05)]
deg_xv_up_p = rownames(res_xv)[which((res_xv$pval.Genotype1XDP_vanilla < 0.05) & (res_xv$LFC.Genotype1XDP_vanilla > 0))]
deg_xv_dn_p = rownames(res_xv)[which((res_xv$pval.Genotype1XDP_vanilla < 0.05) & (res_xv$LFC.Genotype1XDP_vanilla < 0))]

deg_xn = rownames(res_xn)[which(res_xn$FDR.Genotype1XDP_noEdit < 0.1)]
deg_xn_up = rownames(res_xn)[which((res_xn$FDR.Genotype1XDP_noEdit < 0.1) & (res_xn$LFC.Genotype1XDP_noEdit > 0))]
deg_xn_dn = rownames(res_xn)[which((res_xn$FDR.Genotype1XDP_noEdit < 0.1) & (res_xn$LFC.Genotype1XDP_noEdit < 0))]

deg_xn_p = rownames(res_xn)[which(res_xn$pval.Genotype1XDP_noEdit < 0.05)]
deg_xn_up_p = rownames(res_xn)[which((res_xn$pval.Genotype1XDP_noEdit < 0.05) & (res_xn$LFC.Genotype1XDP_noEdit > 0))]
deg_xn_dn_p = rownames(res_xn)[which((res_xn$pval.Genotype1XDP_noEdit < 0.05) & (res_xn$LFC.Genotype1XDP_noEdit < 0))]

deg_xc = rownames(res_xc)[which(res_xc$FDR.GenotypeXDP < 0.1)]
deg_xc_up = rownames(res_xc)[which((res_xc$FDR.GenotypeXDP < 0.1) & (res_xc$LFC.GenotypeXDP > 0))]
deg_xc_dn = rownames(res_xc)[which((res_xc$FDR.GenotypeXDP < 0.1) & (res_xc$LFC.GenotypeXDP < 0))]

deg_xc_p = rownames(res_xc)[which(res_xc$pval.GenotypeXDP < 0.05)]
deg_xc_up_p = rownames(res_xc)[which((res_xc$pval.GenotypeXDP < 0.05) & (res_xc$LFC.GenotypeXDP > 0))]
deg_xc_dn_p = rownames(res_xc)[which((res_xc$pval.GenotypeXDP < 0.05) & (res_xc$LFC.GenotypeXDP < 0))]

#"Replication Set"
deg_ov2 = Reduce(intersect, list(deg_xv, deg_xn))
deg_ov2_up = Reduce(intersect, list(deg_xv_up, deg_xn_up))
deg_ov2_dn = Reduce(intersect, list(deg_xv_dn, deg_xn_dn))
deg_ov2_dir = unique(c(deg_ov2_up, deg_ov2_dn))
deg_ov2_p = Reduce(intersect, list(deg_xv_p, deg_xn_p))
deg_ov2_up_p = Reduce(intersect, list(deg_xv_up_p, deg_xn_up_p))
deg_ov2_dn_p = Reduce(intersect, list(deg_xv_dn_p, deg_xn_dn_p))
deg_ov2_dir_p = unique(c(deg_ov2_up_p, deg_ov2_dn_p))
bg_ov2 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn)))

deg_ov3 = Reduce(intersect, list(deg_xv, deg_xn, deg_xc))
deg_ov3_up = Reduce(intersect, list(deg_xv_up, deg_xn_up, deg_xc_up))
deg_ov3_dn = Reduce(intersect, list(deg_xv_dn, deg_xn_dn, deg_xc_dn))
deg_ov3_dir = unique(c(deg_ov3_up, deg_ov3_dn))
bg_ov3 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc)))
shared_this = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc)))
####################################################################################################################################################################################


####################################################################################################################################################################################
#This study
#Comparisons of DEGs
#############################################
#XDP_noEdit DEGs vs XDP_vanilla DEGs
length(intersect(deg_xn_up_p, deg_xv_up_p))
length(intersect(deg_xn_dn_p, deg_xv_dn_p))
length(intersect(deg_xn_dn_p, deg_xv_up_p))
length(intersect(deg_xn_up_p, deg_xv_dn_p))
length(deg_xv_p)
length(deg_xn_p)
overlap_hgtest(deg_xn_p, deg_xv_p, shared_this)
overlap_hgtest(deg_xn_p, deg_xv_p, shared_this, overlaps = c(intersect(deg_xn_up_p, deg_xv_up_p), intersect(deg_xn_dn_p, deg_xv_dn_p)))
mat1 = res_xv
mat2 = res_xn
val1 = "LFC.Genotype1XDP_vanilla"
val2 = "LFC.Genotype1XDP_noEdit"
filter1 = "pval.Genotype1XDP_vanilla"
filter2 = "pval.Genotype1XDP_noEdit"
name1 = "XDP_unexposed"
name2 = "XDP_unedited"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2, cutoff1 = 0.05, cutoff2 = 0.05)
ggsave(paste0("pdf/corplot_pval_DEG_FC_", name1, "_vs_", name2, ".pdf"),)

length(intersect(deg_xn_up, deg_xv_up))
length(intersect(deg_xn_dn, deg_xv_dn))
length(intersect(deg_xn_dn, deg_xv_up))
length(intersect(deg_xn_up, deg_xv_dn))
length(deg_xv)
length(deg_xn)
overlap_hgtest(deg_xn, deg_xv, shared_this)
overlap_hgtest(deg_xn, deg_xv, shared_this, overlaps = c(intersect(deg_xn_up, deg_xv_up), intersect(deg_xn_dn, deg_xv_dn)))
mat1 = res_xv
mat2 = res_xn
val1 = "LFC.Genotype1XDP_vanilla"
val2 = "LFC.Genotype1XDP_noEdit"
filter1 = "FDR.Genotype1XDP_vanilla"
filter2 = "FDR.Genotype1XDP_noEdit"
name1 = "XDP_unexposed"
name2 = "XDP_unedited"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2)
ggsave(paste0("pdf/corplot_FDR_DEG_FC_", name1, "_vs_", name2, ".pdf"))


#XDP_noEdit DEGs vs XDP DEGs
length(intersect(deg_xn_up_p, deg_xc_up_p))
length(intersect(deg_xn_dn_p, deg_xc_dn_p))
length(intersect(deg_xn_dn_p, deg_xc_up_p))
length(intersect(deg_xn_up_p, deg_xc_dn_p))
length(deg_xc_p)
length(deg_xn_p)
overlap_hgtest(deg_xn_p, deg_xc_p, shared_this)
overlap_hgtest(deg_xn_p, deg_xc_p, shared_this, overlaps = c(intersect(deg_xn_up_p, deg_xc_up_p), intersect(deg_xn_dn_p, deg_xc_dn_p)))
mat1 = res_xc
mat2 = res_xn
val1 = "LFC.GenotypeXDP"
val2 = "LFC.Genotype1XDP_noEdit"
filter1 = "pval.GenotypeXDP"
filter2 = "pval.Genotype1XDP_noEdit"
name1 = "XDP"
name2 = "XDP_unedited"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2, cutoff1 = 0.05, cutoff2 = 0.05)
ggsave(paste0("pdf/corplot_pval_DEG_FC_", name1, "_vs_", name2, ".pdf"))

length(intersect(deg_xn_up, deg_xc_up))
length(intersect(deg_xn_dn, deg_xc_dn))
length(intersect(deg_xn_dn, deg_xc_up))
length(intersect(deg_xn_up, deg_xc_dn))
length(deg_xc)
length(deg_xn)
overlap_hgtest(deg_xn, deg_xc, shared_this)
overlap_hgtest(deg_xn, deg_xc, shared_this, overlaps = c(intersect(deg_xn_up, deg_xc_up), intersect(deg_xn_dn, deg_xc_dn)))
mat1 = res_xc
mat2 = res_xn
val1 = "LFC.GenotypeXDP"
val2 = "LFC.Genotype1XDP_noEdit"
filter1 = "FDR.GenotypeXDP"
filter2 = "FDR.Genotype1XDP_noEdit"
name1 = "XDP"
name2 = "XDP_unedited"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2)
ggsave(paste0("pdf/corplot_FDR_DEG_FC_", name1, "_vs_", name2, ".pdf"))


#XDP_vanilla DEGs vs XDP DEGs
length(intersect(deg_xv_up_p, deg_xc_up_p))
length(intersect(deg_xv_dn_p, deg_xc_dn_p))
length(intersect(deg_xv_dn_p, deg_xc_up_p))
length(intersect(deg_xv_up_p, deg_xc_dn_p))
length(deg_xc_p)
length(deg_xv_p)
overlap_hgtest(deg_xv_p, deg_xc_p, shared_this)
overlap_hgtest(deg_xv_p, deg_xc_p, shared_this, overlaps = c(intersect(deg_xv_up_p, deg_xc_up_p), intersect(deg_xv_dn_p, deg_xc_dn_p)))
mat1 = res_xc
mat2 = res_xv
val1 = "LFC.GenotypeXDP"
val2 = "LFC.Genotype1XDP_vanilla"
filter1 = "pval.GenotypeXDP"
filter2 = "pval.Genotype1XDP_vanilla"
name1 = "XDP"
name2 = "XDP_unexposed"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2, cutoff1 = 0.05, cutoff2 = 0.05)
ggsave(paste0("pdf/corplot_pval_DEG_FC_", name1, "_vs_", name2, ".pdf"))

length(intersect(deg_xv_up, deg_xc_up))
length(intersect(deg_xv_dn, deg_xc_dn))
length(intersect(deg_xv_dn, deg_xc_up))
length(intersect(deg_xv_up, deg_xc_dn))
length(deg_xc)
length(deg_xv)
overlap_hgtest(deg_xv, deg_xc, shared_this)
overlap_hgtest(deg_xv, deg_xc, shared_this, overlaps = c(intersect(deg_xv_up, deg_xc_up), intersect(deg_xv_dn, deg_xc_dn)))
mat1 = res_xc
mat2 = res_xv
val1 = "LFC.GenotypeXDP"
val2 = "LFC.Genotype1XDP_vanilla"
filter1 = "FDR.GenotypeXDP"
filter2 = "FDR.Genotype1XDP_vanilla"
name1 = "XDP"
name2 = "XDP_unexposed"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2)
ggsave(paste0("pdf/corplot_FDR_DEG_FC_", name1, "_vs_", name2, ".pdf"))
####################################################################################################################################################################################


####################################################################################################################################################################################
#This study vs Cell Paper
#(vs Original Cell paper)
#############################################
#Compare to 2018 Cell paper
res_cp = read.table("cellpaper_FDR_LFC_glmm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)
colnames(res_cp)[2] = "LFC.GenotypeXDP"
res_cp$LFC.GenotypeXDP = log(2^res_cp$LFC.GenotypeXDP)
rownames(res_cp) = as.vector(res_cp$name)
shared_all = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc), rownames(res_cp)))
res_cp = res_cp[shared_all, ]

deg_cp = rownames(res_cp)[which(res_cp$FDR.GenotypeXDP < 0.1)]
deg_cp_up = rownames(res_cp)[which((res_cp$FDR.GenotypeXDP < 0.1) & (res_cp$LFC.GenotypeXDP > 0))]
deg_cp_dn = rownames(res_cp)[which((res_cp$FDR.GenotypeXDP < 0.1) & (res_cp$LFC.GenotypeXDP < 0))]

deg_cp_p = rownames(res_cp)[which(res_cp$pval.GenotypeXDP < 0.05)]
deg_cp_up_p = rownames(res_cp)[which((res_cp$pval.GenotypeXDP < 0.05) & (res_cp$LFC.GenotypeXDP > 0))]
deg_cp_dn_p = rownames(res_cp)[which((res_cp$pval.GenotypeXDP < 0.05) & (res_cp$LFC.GenotypeXDP < 0))]

deg_ov_all = Reduce(intersect, list(deg_xv, deg_xn, deg_xc, deg_cp))


#XDP_vanilla DEGs vs Cell Paper DEGs
length(intersect(deg_xv_up, deg_cp_up))
length(intersect(deg_xv_dn, deg_cp_dn))
length(intersect(deg_xv_dn, deg_cp_up))
length(intersect(deg_xv_up, deg_cp_dn))
length(deg_cp)
length(deg_xv)
overlap_hgtest(deg_xv, deg_cp, shared_all)
overlap_hgtest(deg_xv, deg_cp, shared_all, overlaps = c(intersect(deg_xv_up, deg_cp_up), intersect(deg_xv_dn, deg_cp_dn)))


#XDP_noEdit DEGs vs Cell Paper DEGs
length(intersect(deg_xn_up, deg_cp_up))
length(intersect(deg_xn_dn, deg_cp_dn))
length(intersect(deg_xn_dn, deg_cp_up))
length(intersect(deg_xn_up, deg_cp_dn))
length(deg_cp)
length(deg_xn)
overlap_hgtest(deg_xn, deg_cp, shared_all)
overlap_hgtest(deg_xn, deg_cp, shared_all, overlaps = c(intersect(deg_xn_up, deg_cp_up), intersect(deg_xn_dn, deg_cp_dn)))


#XDP DEGs vs Cell Paper DEGs
length(intersect(deg_xc_up, deg_cp_up))
length(intersect(deg_xc_dn, deg_cp_dn))
length(intersect(deg_xc_dn, deg_cp_up))
length(intersect(deg_xc_up, deg_cp_dn))
length(deg_cp)
length(deg_xc)
overlap_hgtest(deg_xc, deg_cp, shared_all)
overlap_hgtest(deg_xc, deg_cp, shared_all, overlaps = c(intersect(deg_xc_up, deg_cp_up), intersect(deg_xc_dn, deg_cp_dn)))

deg_common = c(intersect(deg_xc_up, deg_cp_up), intersect(deg_xc_dn, deg_cp_dn))
ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")
deg_common_gs = ref[deg_common, 1]
bg_gs = ref[shared_this, 1]
go_deg_common = gotest(deg_common_gs, bg_gs)
goplot(go_deg_common[1:5,], cutoff = 0.1, method = "p.value") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))
ggsave("pdf/GO_DEG_common.jpg")

mod_info = read.table("cellpaper_WGCNA_modules.txt", header = T, row.names = 1, sep = "\t", check.names = F)
nsc_mod2 = rownames(mod_info[which(mod_info$NSC_module=="Module_2"), ])
ins_mod5 = rownames(mod_info[which(mod_info$iN_module=="Module_5"), ])
ov_modgene = intersect(nsc_mod2, ins_mod5)
ov_modgene = ov_modgene[ov_modgene %in% shared_this]

mat1 = res_xc
mat2 = res_cp
val1 = "LFC.GenotypeXDP"
val2 = "LFC.GenotypeXDP"
filter1 = "FDR.GenotypeXDP"
filter2 = "FDR.GenotypeXDP"
name1 = "XDP"
name2 = "CellPaper"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2)
ggsave(paste0("pdf/corplot_FDR_DEG_FC_", name1, "_vs_", name2, ".jpg"))


####################################################################################################################################################################################
#This study vs Cell Paper
#(vs Reprocessed Cell paper)
#############################################
#Compare to 2018 Cell paper (re-processed GLMM with newly supported SVA-correction)
res_cp1 = read.table("res_FDR_LFC_glmm_Exp_Genotype_NSC_XDP_vs_CON_cellpaper.txt", header = T, sep = "\t", check.names = F)
rownames(res_cp1) = as.vector(res_cp1$name)
shared_all1 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc), rownames(res_cp1)))
res_cp1 = res_cp1[shared_all1, ]

deg_cp1 = rownames(res_cp1)[which(res_cp1$FDR.GenotypeXDP < 0.1)]
deg_cp1_up = rownames(res_cp1)[which((res_cp1$FDR.GenotypeXDP < 0.1) & (res_cp1$LFC.GenotypeXDP > 0))]
deg_cp1_dn = rownames(res_cp1)[which((res_cp1$FDR.GenotypeXDP < 0.1) & (res_cp1$LFC.GenotypeXDP < 0))]

deg_cp1_p = rownames(res_cp1)[which(res_cp1$pval.GenotypeXDP < 0.05)]
deg_cp1_up_p = rownames(res_cp1)[which((res_cp1$pval.GenotypeXDP < 0.05) & (res_cp1$LFC.GenotypeXDP > 0))]
deg_cp1_dn_p = rownames(res_cp1)[which((res_cp1$pval.GenotypeXDP < 0.05) & (res_cp1$LFC.GenotypeXDP < 0))]

deg_ov_all = Reduce(intersect, list(deg_xv, deg_xn, deg_xc, deg_cp1))


#XDP_vanilla DEGs vs Cell Paper DEGs
length(intersect(deg_xv_up, deg_cp1_up))
length(intersect(deg_xv_dn, deg_cp1_dn))
length(intersect(deg_xv_dn, deg_cp1_up))
length(intersect(deg_xv_up, deg_cp1_dn))
length(deg_cp1)
length(deg_xv)
overlap_hgtest(deg_xv, deg_cp1, shared_all1)
overlap_hgtest(deg_xv, deg_cp1, shared_all1, overlaps = c(intersect(deg_xv_up, deg_cp1_up), intersect(deg_xv_dn, deg_cp1_dn)))


#XDP_noEdit DEGs vs Cell Paper DEGs
length(intersect(deg_xn_up, deg_cp1_up))
length(intersect(deg_xn_dn, deg_cp1_dn))
length(intersect(deg_xn_dn, deg_cp1_up))
length(intersect(deg_xn_up, deg_cp1_dn))
length(deg_cp1)
length(deg_xn)
overlap_hgtest(deg_xn, deg_cp1, shared_all1)
overlap_hgtest(deg_xn, deg_cp1, shared_all1, overlaps = c(intersect(deg_xn_up, deg_cp1_up), intersect(deg_xn_dn, deg_cp1_dn)))


#XDP DEGs vs Cell Paper DEGs
length(intersect(deg_xc_up, deg_cp1_up))
length(intersect(deg_xc_dn, deg_cp1_dn))
length(intersect(deg_xc_dn, deg_cp1_up))
length(intersect(deg_xc_up, deg_cp1_dn))
length(deg_cp1)
length(deg_xc)
overlap_hgtest(deg_xc, deg_cp1, shared_all1)
overlap_hgtest(deg_xc, deg_cp1, shared_all1, overlaps = c(intersect(deg_xc_up, deg_cp1_up), intersect(deg_xc_dn, deg_cp1_dn)))

deg_common1 = c(intersect(deg_xc_up, deg_cp1_up), intersect(deg_xc_dn, deg_cp1_dn))
ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")
deg_common1_gs = ref[deg_common1, 1]
bg_gs1 = ref[shared_all1, 1]
go_deg_common1 = gotest(deg_common1_gs, bg_gs1)
goplot(go_deg_common1[1:5,], cutoff = 0.1, method = "p.value") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))
ggsave("pdf/GO_DEG_common1.jpg", width = 5, height = 5)

mat1 = res_xc
mat2 = res_cp1
val1 = "LFC.GenotypeXDP"
val2 = "LFC.GenotypeXDP"
filter1 = "FDR.GenotypeXDP"
filter2 = "FDR.GenotypeXDP"
name1 = "XDP"
name2 = "CellPaper"
stat_corplot(mat1 = mat1, mat2 = mat2, val1 = val1, val2 = val2, filter1 = filter1, filter2 = filter2, name1 = name1, name2 = name2)
ggsave(paste0("pdf/corplot_FDR_DEG_FC_", name1, "_vs_", name2, "_new.jpg"), width = 5, height = 5)
####################################################################################################################################################################################


####################################################################################################################################################################################
#This study
#Determine useful DEG lists
#############################################
res_xv = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPvanilla_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_xn = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPnoEdit_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_xc = read.table("res_FDR_LFC_glm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)

rownames(res_xv) = as.vector(res_xv$name)
rownames(res_xn) = as.vector(res_xn$name)
rownames(res_xc) = as.vector(res_xc$name)

deg_xv = rownames(res_xv)[which(res_xv$FDR.Genotype1XDP_vanilla < 0.1)]
deg_xv_up = rownames(res_xv)[which((res_xv$FDR.Genotype1XDP_vanilla < 0.1) & (res_xv$LFC.Genotype1XDP_vanilla > 0))]
deg_xv_dn = rownames(res_xv)[which((res_xv$FDR.Genotype1XDP_vanilla < 0.1) & (res_xv$LFC.Genotype1XDP_vanilla < 0))]

deg_xv_p = rownames(res_xv)[which(res_xv$pval.Genotype1XDP_vanilla < 0.05)]
deg_xv_up_p = rownames(res_xv)[which((res_xv$pval.Genotype1XDP_vanilla < 0.05) & (res_xv$LFC.Genotype1XDP_vanilla > 0))]
deg_xv_dn_p = rownames(res_xv)[which((res_xv$pval.Genotype1XDP_vanilla < 0.05) & (res_xv$LFC.Genotype1XDP_vanilla < 0))]

deg_xn = rownames(res_xn)[which(res_xn$FDR.Genotype1XDP_noEdit < 0.1)]
deg_xn_up = rownames(res_xn)[which((res_xn$FDR.Genotype1XDP_noEdit < 0.1) & (res_xn$LFC.Genotype1XDP_noEdit > 0))]
deg_xn_dn = rownames(res_xn)[which((res_xn$FDR.Genotype1XDP_noEdit < 0.1) & (res_xn$LFC.Genotype1XDP_noEdit < 0))]

deg_xn_p = rownames(res_xn)[which(res_xn$pval.Genotype1XDP_noEdit < 0.05)]
deg_xn_up_p = rownames(res_xn)[which((res_xn$pval.Genotype1XDP_noEdit < 0.05) & (res_xn$LFC.Genotype1XDP_noEdit > 0))]
deg_xn_dn_p = rownames(res_xn)[which((res_xn$pval.Genotype1XDP_noEdit < 0.05) & (res_xn$LFC.Genotype1XDP_noEdit < 0))]

deg_xc = rownames(res_xc)[which(res_xc$FDR.GenotypeXDP < 0.1)]
deg_xc_up = rownames(res_xc)[which((res_xc$FDR.GenotypeXDP < 0.1) & (res_xc$LFC.GenotypeXDP > 0))]
deg_xc_dn = rownames(res_xc)[which((res_xc$FDR.GenotypeXDP < 0.1) & (res_xc$LFC.GenotypeXDP < 0))]

deg_xc_p = rownames(res_xc)[which(res_xc$pval.GenotypeXDP < 0.05)]
deg_xc_up_p = rownames(res_xc)[which((res_xc$pval.GenotypeXDP < 0.05) & (res_xc$LFC.GenotypeXDP > 0))]
deg_xc_dn_p = rownames(res_xc)[which((res_xc$pval.GenotypeXDP < 0.05) & (res_xc$LFC.GenotypeXDP < 0))]

deg_ov2 = Reduce(intersect, list(deg_xv, deg_xn))
deg_ov2_up = Reduce(intersect, list(deg_xv_up, deg_xn_up))
deg_ov2_dn = Reduce(intersect, list(deg_xv_dn, deg_xn_dn))
deg_ov2_dir = unique(c(deg_ov2_up, deg_ov2_dn))

deg_ov2_p = Reduce(intersect, list(deg_xv_p, deg_xn_p))
deg_ov2_up_p = Reduce(intersect, list(deg_xv_up_p, deg_xn_up_p))
deg_ov2_dn_p = Reduce(intersect, list(deg_xv_dn_p, deg_xn_dn_p))
deg_ov2_dir_p = unique(c(deg_ov2_up_p, deg_ov2_dn_p))
bg_ov2 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn)))

deg_ov3 = Reduce(intersect, list(deg_xv, deg_xn, deg_xc))
deg_ov3_up = Reduce(intersect, list(deg_xv_up, deg_xn_up, deg_xc_up))
deg_ov3_dn = Reduce(intersect, list(deg_xv_dn, deg_xn_dn, deg_xc_dn))
deg_ov3_dir = unique(c(deg_ov3_up, deg_ov3_dn))
bg_ov3 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc)))
shared_this = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc)))

#Permissive Set: XDP_vanilla vs CON at p-value significance
perm = deg_xc_p
perm_up = deg_xc_up_p
perm_dn = deg_xc_dn_p
perm_bg = rownames(res_xc)

#Replication Set: (DEGs of XDP_vanilla vs CON at p-value significance) overlap with (DEGs of XDP_noEdit vs CON at p-value significance)
repl = deg_ov2_dir_p
repl_up = deg_ov2_up_p
repl_dn = deg_ov2_dn_p
repl_bg = bg_ov2

#Stringent Set: (DEGs of XDP_vanilla vs CON at FDR significance) overlap with (DEGs of XDP_noEdit vs CON at FDR significance)
stri = deg_ov2_dir
stri_up = deg_ov2_up
stri_dn = deg_ov2_dn
stri_bg = bg_ov2


####################################################################################################################################################################################
#This study
#XDP vs CON Pathway analyses
#############################################
godb_list = list(
  GO_BP = "~/ref/gsea/gsea_go_bp.rds",
  GO_CC = "~/ref/gsea/gsea_go_cc.rds",
  GO_MF = "~/ref/gsea/gsea_go_mf.rds",
  PATHWAYS = "~/ref/gsea/gsea_pathways.rds",
  REG = "~/ref/gsea/gsea_reg.rds",
  HPO = "~/ref/gsea/gsea_hpo.rds"
)
ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")

perm_gs = ref[perm, 1]
perm_up_gs = ref[perm_up, 1]
perm_dn_gs = ref[perm_dn, 1]
perm_bg_gs = ref[perm_bg, 1]

repl_gs = ref[repl, 1]
repl_up_gs = ref[repl_up, 1]
repl_dn_gs = ref[repl_dn, 1]
repl_bg_gs = ref[repl_bg, 1]

stri_gs = ref[stri, 1]
stri_up_gs = ref[stri_up, 1]
stri_dn_gs = ref[stri_dn, 1]
stri_bg_gs = ref[stri_bg, 1]

#Permissive Set
perm_go = go_test(perm_gs, perm_bg_gs, godb_list)
perm_up_go = go_test(perm_up_gs, perm_bg_gs, godb_list)
perm_dn_go = go_test(perm_dn_gs, perm_bg_gs, godb_list)
write.table(perm_go, file = paste0("GSEA_GO_XDPDEG_","permissive_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(perm_up_go, file = paste0("GSEA_GO_XDPDEG_","permissive_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(perm_dn_go, file = paste0("GSEA_GO_XDPDEG_","permissive_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#Replication Set
repl_go = go_test(repl_gs, repl_bg_gs, godb_list)
repl_up_go = go_test(repl_up_gs, repl_bg_gs, godb_list)
repl_dn_go = go_test(repl_dn_gs, repl_bg_gs, godb_list)
write.table(repl_go, file = paste0("GSEA_GO_XDPDEG_","replication_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(repl_up_go, file = paste0("GSEA_GO_XDPDEG_","replication_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(repl_dn_go, file = paste0("GSEA_GO_XDPDEG_","replication_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#Stringent Set
stri_go = go_test(stri_gs, stri_bg_gs, godb_list)
stri_up_go = go_test(stri_up_gs, stri_bg_gs, godb_list)
stri_dn_go = go_test(stri_dn_gs, stri_bg_gs, godb_list)
write.table(stri_go, file = paste0("GSEA_GO_XDPDEG_","stringent_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(stri_up_go, file = paste0("GSEA_GO_XDPDEG_","stringent_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(stri_dn_go, file = paste0("GSEA_GO_XDPDEG_","stringent_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#Read in GO results
perm_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_all",".txt"), header = T, row.names = 1, sep = "\t")
perm_up_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_up",".txt"), header = T, row.names = 1, sep = "\t")
perm_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_dn",".txt"), header = T, row.names = 1, sep = "\t")
repl_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_all",".txt"), header = T, row.names = 1, sep = "\t")
repl_up_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_up",".txt"), header = T, row.names = 1, sep = "\t")
repl_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_dn",".txt"), header = T, row.names = 1, sep = "\t")
stri_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_all",".txt"), header = T, row.names = 1, sep = "\t")
stri_up_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_up",".txt"), header = T, row.names = 1, sep = "\t")
stri_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_dn",".txt"), header = T, row.names = 1, sep = "\t")

#force the same order
perm_go = perm_go[rownames(perm_go), ]
perm_up_go = perm_up_go[rownames(perm_go), ]
perm_dn_go = perm_dn_go[rownames(perm_go), ]
repl_go = repl_go[rownames(perm_go), ]
repl_up_go = repl_up_go[rownames(perm_go), ]
repl_dn_go = repl_dn_go[rownames(perm_go), ]
stri_go = stri_go[rownames(perm_go), ]
stri_up_go = stri_up_go[rownames(perm_go), ]
stri_dn_go = stri_dn_go[rownames(perm_go), ]

#extract terms of interest
aloy_vec = scan("go_terms_roi_v2.txt", what = "character")
keyword_vec = c(aloy_vec, "NEURO", "SYNAP")
perm_go_sel = go_extract_keywords(perm_go[which(perm_go$pvalue < 0.05), ], "term", keyword_vec)
repl_go_sel = go_extract_keywords(repl_go[which(repl_go$pvalue < 0.05), ], "term", keyword_vec)
stri_go_sel = go_extract_keywords(stri_go[which(stri_go$pvalue < 0.05), ], "term", keyword_vec)
perm_dn_go_sel = go_extract_keywords(perm_dn_go[which(perm_dn_go$pvalue < 0.05), ], "term", keyword_vec)
repl_dn_go_sel = go_extract_keywords(repl_dn_go[which(repl_dn_go$pvalue < 0.05), ], "term", keyword_vec)
stri_dn_go_sel = go_extract_keywords(stri_dn_go[which(stri_dn_go$pvalue < 0.05), ], "term", keyword_vec)
perm_up_go_sel = go_extract_keywords(perm_up_go[which(perm_up_go$pvalue < 0.05), ], "term", keyword_vec)
repl_up_go_sel = go_extract_keywords(repl_up_go[which(repl_up_go$pvalue < 0.05), ], "term", keyword_vec)
stri_up_go_sel = go_extract_keywords(stri_up_go[which(stri_up_go$pvalue < 0.05), ], "term", keyword_vec)

perm_go_sel0 = go_extract_keywords(perm_go, "term", keyword_vec)
repl_go_sel0 = go_extract_keywords(repl_go, "term", keyword_vec)
stri_go_sel0 = go_extract_keywords(stri_go, "term", keyword_vec)
perm_dn_go_sel0 = go_extract_keywords(perm_dn_go, "term", keyword_vec)
repl_dn_go_sel0 = go_extract_keywords(repl_dn_go, "term", keyword_vec)
stri_dn_go_sel0 = go_extract_keywords(stri_dn_go, "term", keyword_vec)
perm_up_go_sel0 = go_extract_keywords(perm_up_go, "term", keyword_vec)
repl_up_go_sel0 = go_extract_keywords(repl_up_go, "term", keyword_vec)
stri_up_go_sel0 = go_extract_keywords(stri_up_go, "term", keyword_vec)

perm_go_sig = perm_go[which(perm_go$pvalue < 0.05), ]
repl_go_sig = repl_go[which(repl_go$pvalue < 0.05), ]
stri_go_sig = stri_go[which(stri_go$pvalue < 0.05), ]
perm_dn_go_sig = perm_dn_go[which(perm_dn_go$pvalue < 0.05), ]
repl_dn_go_sig = repl_dn_go[which(repl_dn_go$pvalue < 0.05), ]
stri_dn_go_sig = stri_dn_go[which(stri_dn_go$pvalue < 0.05), ]
perm_up_go_sig = perm_up_go[which(perm_up_go$pvalue < 0.05), ]
repl_up_go_sig = repl_up_go[which(repl_up_go$pvalue < 0.05), ]
stri_up_go_sig = stri_up_go[which(stri_up_go$pvalue < 0.05), ]

overlap_hgtest(rownames(perm_go_sig), rownames(perm_go_sel0), rownames(perm_go))
overlap_hgtest(rownames(repl_go_sig), rownames(repl_go_sel0), rownames(repl_go))
overlap_hgtest(rownames(stri_go_sig), rownames(stri_go_sel0), rownames(stri_go))

overlap_hgtest(rownames(perm_up_go_sig), rownames(perm_up_go_sel0), rownames(perm_up_go))
overlap_hgtest(rownames(repl_up_go_sig), rownames(repl_up_go_sel0), rownames(repl_up_go))
overlap_hgtest(rownames(stri_up_go_sig), rownames(stri_up_go_sel0), rownames(stri_up_go))

overlap_hgtest(rownames(perm_dn_go_sig), rownames(perm_dn_go_sel0), rownames(perm_dn_go))
overlap_hgtest(rownames(repl_dn_go_sig), rownames(repl_dn_go_sel0), rownames(repl_dn_go))
overlap_hgtest(rownames(stri_dn_go_sig), rownames(stri_dn_go_sel0), rownames(stri_dn_go))




#non-directional DEGs for GO
go_plot(perm_go_sel[order(perm_go_sel$pvalue), ][1:min(nrow(perm_go_sel), 20), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","permissive_all",".jpg"), width = 6.1, height = 4.24)
go_plot(repl_go_sel[order(repl_go_sel$pvalue), ][1:min(nrow(repl_go_sel), 20), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","replication_all",".jpg"), width = 6.1, height = 4.24)
go_plot(stri_go_sel[order(stri_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","stringent_all",".jpg"), width = 6.1, height = 4.24)
intersect(perm_go_sel[order(perm_go_sel$pvalue), ][1:min(nrow(perm_go_sel), 20), "term"], repl_go_sel[order(repl_go_sel$pvalue), ][1:min(nrow(repl_go_sel), 20), "term"])

#Downregulated DEGs for GO
go_plot(perm_dn_go_sel[order(perm_dn_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","permissive_dn",".jpg"), width = 6.1, height = 4.24)
go_plot(repl_dn_go_sel[order(repl_dn_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","replication_dn",".jpg"), width = 6.1, height = 4.24)
go_plot(stri_dn_go_sel[order(stri_dn_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","stringent_dn",".jpg"), width = 6.1, height = 4.24)
intersect(perm_dn_go_sel[order(perm_dn_go_sel$pvalue), ][, "term"], repl_dn_go_sel[order(repl_dn_go_sel$pvalue), ][, "term"])
intersect(stri_dn_go_sel[order(stri_dn_go_sel$pvalue), ][, "term"], repl_dn_go_sel[order(repl_dn_go_sel$pvalue), ][, "term"])


#Upregulated DEGs for GO
go_plot(perm_up_go_sel[order(perm_up_go_sel$pvalue), ][1:min(nrow(perm_up_go_sel), 18), ], "pvalue", 0.05, text_width = 30, text_size = 5.75, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","permissive_up",".jpg"), width = 6.1, height = 4.24)
go_plot(repl_up_go_sel[order(repl_up_go_sel$pvalue), ][1:min(nrow(repl_up_go_sel), 18), ], "pvalue", 0.05, text_width = 30, text_size = 5.75, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","replication_up",".jpg"), width = 6.1, height = 4.24)
go_plot(stri_up_go_sel[order(stri_up_go_sel$pvalue), ][1:min(nrow(stri_up_go_sel), 18), ], "pvalue", 0.05, text_width = 30, text_size = 5.75, show_hit = F)
ggsave(paste0("pdf/GOplot_XDPDEG_","stringent_up",".jpg"), width = 6.1, height = 4.24)
intersect(perm_up_go_sel[order(perm_up_go_sel$pvalue), ][1:min(nrow(perm_up_go_sel), 20), "term"], repl_up_go_sel[order(repl_up_go_sel$pvalue), ][1:min(nrow(repl_up_go_sel), 20), "term"])
Reduce(intersect, list(perm_up_go_sel[order(perm_up_go_sel$pvalue), ][1:min(nrow(perm_up_go_sel), 20), "term"], repl_up_go_sel[order(repl_up_go_sel$pvalue), ][1:min(nrow(repl_up_go_sel), 20), "term"], stri_up_go_sel[order(stri_up_go_sel$pvalue), ][1:min(nrow(stri_up_go_sel), 20), "term"]))










#############################################
#Treated
#############################################
key_name = "ENSG00000147133"

for (aso in  c("880", "131", "307", "877", "876", "881", "879")){
  grp = read.table(paste0("model_params_colData_glmm_Exp_Treatment_",aso,".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
  grp$MIN = factor(grp$MIN)
  mat = read.table(paste0("model_params_deSV_glmm_Exp_Treatment_",aso,".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
  rnaPCA(as.matrix(mat), grp, c("Treatment"), textLabel = "ParentLine", ntop = 500, PC1 = 3, PC2 = 4)
  ggsave(paste0("pdf/PCA_rnaseq_Exp_untreated_XDP_vs_ASO", aso, "_new.jpg"))
  
  bl1 = read.table(paste0("res_FDR_LFC_glm_Exp_Genotype_CON_vs_XDP.txt"), header = T, sep = "\t")
  rownames(bl1) = as.vector(bl1$name)
  bl2 = read.table(paste0("res_FDR_LFC_glmm_Exp_Treatment_",aso,".txt"), header = T, sep = "\t")
  rownames(bl2) = as.vector(bl2$name)
  bl3 = read.table(paste0("res_FDR_LFC_crossModel_Exp_Treatment",aso,"_vs_GenotypeCON.txt"), header = T, sep = "\t")
  rownames(bl3) = as.vector(bl3$name)
  shared_genes = intersect(intersect(rownames(bl1), rownames(bl2)), rownames(bl3))
  new = cbind(bl1[shared_genes, 1:2], bl2[shared_genes, 1:2], bl3[shared_genes, 1:2])
  new$status = "Other"
  new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.Treatment < 0.1) & (sign(new$LFC.GenotypeCON) == sign(new$LFC.Treatment)), "status"] = "Rescued"
  new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.Treatment < 0.1) & (sign(new$LFC.GenotypeCON) != sign(new$LFC.Treatment)), "status"] = "Un-rescued"
  new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.Treatment >= 0.1), "status"] = "Un-rescued"
  new[(new[[paste0("FDR.Treatment",aso)]] < 0.1) & (new$status != "Rescued"), "status"] = "Off-target"
  table(new$status)
  write.table(as.data.frame(rownames(new)[new$status == "Rescued"]), file = paste0("ASO",aso,"_Rescued_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(rownames(new)[new$status == "Un-rescued"]), file = paste0("ASO",aso,"_Unrescued_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(rownames(new)[new$status == "Off-target"]), file = paste0("ASO",aso,"_Offtarget_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(rownames(new)), file = paste0("ASO",aso,"_background_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  
  ggplot(tt, aes(x = val0, y = val1, color = color)) + geom_point() + 
    xlim(c(-60,60)) + ylim(c(-60,60)) + 
    geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = c("red","black","red")) + 
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed", color = c("red","black","red")) + 
    xlab(paste0("Significance of IR Changes in ASO", aso, " (directional)")) + 
    ylab("Significance of IR Changes in dSVA (directional)") +
    geom_color(c("blue","red")) + geom_noBG()
  ggsave(paste0("pdf/scatter_Exp_change_treated_ASO", aso, "_new.jpg"))
}


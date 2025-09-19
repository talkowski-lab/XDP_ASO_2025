library(DESeq2)
library(data.table)
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/countsFilter.R")
source("~/lib/R/geom_noBG.R")
source("~/projects/xdp_aso/nsc/scripts/glmm.R")
source("~/projects/xdp_aso/nsc/scripts/go.R")


setwd("~/projects/xdp_aso/nsc/output/")

key_name = "ENSG00000147133"

meta = read.table("ASO_RNAseq_new.569_namecorrected.metadata.tab", header = T, check.names = F, sep = "\t", row.names = 1)
raw = read.table("ASO_RNAseq_new.exp.txt", header = T, check.names = F)
libsz = read.table("ASO_RNAseq_new.libSize.txt", header = F, check.names = F, row.names = 1)
grp = meta[colnames(raw), ]
grp$libsz = as.vector(libsz[rownames(grp), "V2"])

grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$geno %in% c("XDP_non_edit", "XDP", "CON")), ]
grp$Genotype = as.vector(grp$geno)
grp$Genotype[which(grp$geno == "XDP" | grp$geno == "XDP_non_edit")] = "XDP"
#grp$Genotype[which(grp$geno == "XDP_dSVA")] = "dSVA"
grp$geno = factor(grp$geno, levels = c("CON", "XDP", "XDP_non_edit"))
grp$Genotype = factor(grp$Genotype, levels = c("CON", "XDP"))
#grp$Genotype = factor(grp$Genotype, levels = c("CON", "XDP", "dSVA"))
table(grp$Genotype)

grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$ParentLine = factor(grp$ParentLine)
grp$MIN = factor(grp$MIN)

mat = raw[,rownames(grp)]
CPM = t(1000000 * t(mat) / grp$libsz)
target_col = "Genotype"
condition_list = lapply(names(table(grp[[target_col]])), FUN = function(x){which(grp[[target_col]] == x)})
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list)
mat = mat[which(expressedGenes), ]
message(nrow(mat)," genes expressed.")

ms = "Exp"
#mf = "geno"
mf = "Genotype"
#pf = paste0("XDP_non_edit", "_vs_CON")
pf = "baseline_vs_CON"

#method0 = "deseq2"
#test0 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method0, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")
#resultsNames(test0)
#res0 = results(test0, contrast = c(0, 1, rep(0, 7)))

method1 = "glm"
test1 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method1, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")

write.table(test1@vcov, paste0("model_params_vcov_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test1@colData, paste0("model_params_colData_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@beta, paste0("model_params_beta_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@stderr, paste0("model_params_stderr_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@pval, paste0("model_params_pval_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@fdr, paste0("model_params_fdr_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@deSV, paste0("model_params_deSV_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@resid, paste0("model_params_resid_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@status, paste0("model_params_warnings_",method1,"_",ms,"_",mf,"_",pf,"_new.txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#res1 = get_fdr(res = test1, term = "geno", class = pf)
res1 = get_fdr(res = test1, term = "Genotype", class = "XDP_vs_CON")
#length(res1[which(res1$FDR.genoXDP < 0.1), "LFC.genoXDP"])
ci1 = get_ci(res = test1, term = "GenotypeXDP", class = "XDP_vs_CON")
#ci1d = get_ci(res = test1, term = "GenotypedSVA", class = "dSVA_vs_CON")

#write.table(res1, file = "res_FDR_LFC_glm_Exp_geno_baseline_new.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(res1, file = paste0("res_FDR_LFC_",method1,"_",ms,"_",mf,"_","XDP_vs_CON","_new.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(ci1, file = paste0("res_LFC_CI_",method1,"_",ms,"_",mf,"_","XDP_vs_CON","_new.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
#write.table(ci1d, file = "res_LFC_CI_glm_Exp_Genotype_dSVA_vs_CON_new2.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#method2 = "glm"
#test2 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method2, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth"), remove_terms = c("BatchGrowth"), random_terms = NULL, distro_family = "nb")
#res1 = get_fdr(res = test1, term = "geno", class = pf)
#res2 = get_fdr(res = test2, term = "Genotype", class = pf)


#rnaPCA(test0@metadata$deSV, grp, "geno")
#rnaPCA(test1@deSV, grp, "geno")
rnaPCA(test1@deSV, grp, "Genotype")
#rnaPCA(test2@deSV, grp, "Genotype")
#ggsave("pdf/correctedPCA_glm_Exp_geno_baseline_new.pdf")
ggsave("pdf/correctedPCA_glm_Exp_Genotype_baseline_new.pdf")
#rnaPCA(test2@deSV, grp, "geno")
rnaPCA(test1@deSV, grp, "geno")
ggsave("pdf/correctedPCA_glm_Exp_Genotype_baseline_groupBygeno_new.pdf")


df = cbind(Exp = test1@deSV[key_name,], test1@colData)
ggplot(df, aes(x = Genotype, y = Exp)) + geom_boxplot(outlier.shape = NA) + geom_noBG() + xlab("") + ylab("Expression (log-scale)")
ggsave("pdf/boxplot_correctedExpression_TAF1_glm_Exp_Genotype_baseline_new.pdf")
ggplot(df, aes(x = geno, y = Exp)) + geom_boxplot(outlier.shape = NA) + geom_noBG() + xlab("") + ylab("Expression (log-scale)") + scale_x_discrete(labels=c("CON", "XDP_unexposed", "XDP_unedited"))
ggsave("pdf/boxplot_correctedExpression_TAF1_glm_Exp_Genotype_baseline_groupBygeno_new.pdf")


#===================================
#XDP as baseline
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))


ms = "Exp"
mf = "Genotype"
pf = "baseline_vs_XDP"

method1v = "glm"
test1v = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method1v, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")
res1v = get_fdr(res = test1v, term = "Genotype", class = pf)
#ci1v = get_ci(res = test1v, term = "GenotypeCON", class = "CON_vs_XDP")
ci1dv = get_ci(res = test1v, term = "GenotypedSVA", class = "dSVA_vs_XDP")

write.table(res1v, file = "res_FDR_LFC_glm_Exp_Genotype_baseline_vs_XDP_new2.txt", col.names = T, row.names = F, quote = F, sep = "\t")
#write.table(ci1v, file = "res_LFC_CI_glm_Exp_Genotype_CON_vs_XDP_new2.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(ci1dv, file = "res_LFC_CI_glm_Exp_Genotype_dSVA_vs_XDP_new2.txt", col.names = T, row.names = F, quote = F, sep = "\t")

rnaPCA(test1@deSV, grp, "Genotype")
ggsave("pdf/correctedPCA_glm_Exp_Genotype_baseline_new2.pdf")

save(test1, test1v, file = "RData/glm_Exp_Genotype_baseline_new2.RData")


#===================================

godb_list = list(
  GO_BP = "~/ref/gsea/human_gsea_go_bp.rds",
  GO_CC = "~/ref/gsea/human_gsea_go_cc.rds",
  GO_MF = "~/ref/gsea/human_gsea_go_mf.rds",
  PATHWAYS = "~/ref/gsea/human_gsea_pathways.rds",
  REG = "~/ref/gsea/human_gsea_reg.rds",
  HPO = "~/ref/gsea/human_gsea_hpo.rds"
)
ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")

deg = rownames(res0[res0$padj<0.1,])
bg = rownames(res0)
deg_gs = ref[deg, 1]
bg_gs = ref[bg, 1]
deg_go = go_test(deg_gs, bg_gs, godb_list)
head(deg_go[,1:6],10)




#===================================
setwd("~/projects/xdp_aso/nsc/output/")

key_name = "ENSG00000147133"

fc_cut = 1.2
sig_cut = 0.05

res_880 = read.table("res_FDR_LFC_glm_Exp_Trt_XDP880_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_131 = read.table("res_FDR_LFC_glm_Exp_Trt_XDP131_vs_CON.txt", header = T, sep = "\t", check.names = F)
a880_uc = as.vector(res_880[which(res_880[["pval.Treatment"]] >= sig_cut), "name"])
a131_uc = as.vector(res_131[which(res_131[["pval.Treatment"]] >= sig_cut), "name"])

con_mode = "baseline0"
if (con_mode == "baseline1"){
  res_0 = read.table("res_FDR_LFC_glm_Exp_Genotype_baseline_vs_CON_new2.txt", header = T, sep = "\t", check.names = F)
  xdp_deg = as.vector(res_0[which(res_0[["pval.GenotypeXDP"]] < sig_cut & abs(res_0[["LFC.GenotypeXDP"]]) >= log(fc_cut)), "name"])
} else if (con_mode == "baseline0") {
  res_0 = read.table("res_FDR_LFC_glm_Exp_Genotype_XDP_vs_CON_new.txt", header = T, sep = "\t", check.names = F)
  xdp_deg = as.vector(res_0[which(res_0[["pval.GenotypeXDP"]] < sig_cut & abs(res_0[["LFC.GenotypeXDP"]]) >= log(fc_cut)), "name"])
} else {
  res_0 = read.table("res_FDR_LFC_glm_Exp_Trt_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)
  xdp_deg = as.vector(res_0[which(res_0[["pval.Treatment"]] < sig_cut & abs(res_0[["LFC.Treatment"]]) >= log(fc_cut)), "name"])
}


dsva_mode = "all-in-one"
if (dsva_mode == "baseline"){
  res_dsva = read.table("res_FDR_LFC_glm_Exp_Genotype_baseline_vs_CON_new2.txt", header = T, sep = "\t", check.names = F)
  res_dsva = res_dsva[, c(2, 4, 6, 7, 8)]
  dsva_uc = as.vector(res_dsva[which(res_dsva[["pval.GenotypedSVA"]] >= sig_cut), "name"])
} else {
  res_dsva = read.table("res_FDR_LFC_glm_Exp_Trt_dSVA_vs_CON.txt", header = T, sep = "\t", check.names = F)
  dsva_uc = as.vector(res_dsva[which(res_dsva[["pval.Treatment"]] >= sig_cut), "name"])
}


length(intersect(xdp_deg, dsva_uc)) / length(xdp_deg)
length(intersect(xdp_deg, a880_uc)) / length(xdp_deg)
length(intersect(xdp_deg, a131_uc)) / length(xdp_deg)



df = c()
ci0 = read.table("res_LFC_CI_glm_Exp_Genotype_XDP_vs_CON_new.txt", header = T, sep = "\t", check.names = F)
taf1_base = ci0[as.vector(ci0$name) == key_name, ]
taf1_base$group = "XDP"
df = rbind(df, taf1_base)

for (x in c("dSVA", "XDP880", "XDP131", "XDP307", "XDP881", "XDP877", "XDP879", "XDP876")){
  ci1 = read.table(paste0("res_LFC_CI_glm_Exp_Trt_", x, "_vs_CON.txt"), header = T, sep = "\t", check.names = F)
  colnames(ci1) = c("LFC", "margin", "name", "class")
  taf1 = ci1[as.vector(ci1$name) == key_name, ]
  taf1$group = x
  df = rbind(df, taf1)
}

df$group = factor(df$group, levels = as.vector(df$group))
ggplot(data = df, aes(x = group , y = exp(LFC), ymin = exp(LFC-margin), ymax = exp(LFC+margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = exp(df[df$group=="dSVA","LFC"]), linetype="dashed", color = "purple") +
  geom_hline(yintercept = exp(df[df$group=="XDP","LFC"]), linetype="dashed", color = "red") +
  geom_hline(yintercept = 1, linetype="dashed", color = "black") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0.65, 1.15), breaks = seq(0.65, 1.15, 0.1))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_rnaseq_TAF1_FC_to_untreated_CON_corrected_Trt_new.pdf", width = 5.7, height = 4.2)

library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
library(ggbeeswarm)
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")

#################################################
#Use untreated XDP as baseline
#################################################
#################################################
#Building AS models for baseline: CON/dSVA vs XDP
#################################################
#i32AS
key_name = "X_70644089_70659646_1"
grp = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
#grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$LibType != "ILL"), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", "CON", "dSVA"))
grp$BatchGrowth = factor(grp$BatchGrowth)
table(paste0(grp$Genotype,grp$ParentLine))

mat0 = read.table("ASO_Capseq_R1.i32splicing.txt", header = T, row.names = 1, sep = "\t")
mat = mat0[,rownames(grp)]
norm = log1p(t(t(mat)/grp$sizeFactor))
#rnaPCA(norm,grp,c("Genotype"),ntop = 5)
new = cbind(Exp = unlist(norm[1, ]), grp)
write.table(new, file = "ASO_Capseq_R1.metadata_with_i32splicing.txt", row.names = T, col.names = T, quote = F, sep = "\t")
ggplot(new, aes(Genotype, log(exp(Exp)))) + geom_boxplot(outlier.alpha = 0) +  geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_raw_new.jpg")

test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glm", estimate_offset = F, offset_str = "sizeFactor", apply_sva = F, test_terms = c("Genotype", "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "nb")

norm1 = test@deSV
#rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
new1 = cbind(Exp = unlist(norm[1, ]), grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_corrected_new.jpg")
ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_corrected_v2_new.jpg")
new2 = new1[new1$Genotype=="XDP"|(new1$Genotype=="CON" & new1$Exp == 0)|(new1$Genotype=="dSVA" & new1$Exp == 0),]
ggplot(new2, aes(Genotype, log2(exp(Exp)))) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_corrected_outlierSamplesRemoved_new.jpg")

get_fdr(res = test, term = "Genotype", class = "test", keyword = key_name)

bl1 = get_ci(res = test, term = "GenotypedSVA", class = "dSVA")
bl1 = bl1[1,]
bl2 = get_ci(res = test, term = "GenotypeCON", class = "CON")
bl2 = bl2[1,]
bl = rbind(bl1, bl2)


#################################################
#Building AS models for treated and untreated XDP
#################################################
#Treatment
key_name = "X_70644089_70659646_1"
trt = c()
stats = c()
for (aso in c("126", "131", "134", "307", "849", "852", "876", "877", "879", "880", "881", "883", "MAL")){
  message("ASO",aso)
  grp1 = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  grp1 = grp1[,c(-5, -6)]
  #grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
  grp1 = grp1[which(grp1$LibType != "ILL"), ]
  grp1 = grp1[which(grp1$Treatment %in% c("NO", aso)), ]
  grp1 = grp1[which(grp1$Genotype %in% c("XDP")), ]
  #grp$Treatment = factor(grp$Treatment,levels=c("NO", "307", "876", "877", "879", "880", "881"))
  grp1$Treatment = factor(grp1$Treatment,levels=c("NO", aso))
  grp1$MIN = factor(grp1$MIN)
  #table(paste0(grp$Treatment,"_",grp$ParentLine))
  
  grp2 = read.table("ASO_Capseq_R2.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  grp2 = grp2[,c(-5, -6)]
  #grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
  grp2 = grp2[which(grp2$LibType != "ILL"), ]
  grp2 = grp2[which(grp2$Treatment %in% c("NO", aso)), ]
  grp2 = grp2[which(grp2$Genotype %in% c("XDP")), ]
  #grp$Treatment = factor(grp$Treatment,levels=c("NO", "307", "876", "877", "879", "880", "881"))
  grp2$Treatment = factor(grp2$Treatment,levels=c("NO", aso))
  grp2$MIN = factor(grp2$MIN)
  
  rownames(grp1) = paste0(rownames(grp1), "_R1")
  rownames(grp2) = paste0(rownames(grp2), "_R2")
  
  mat1 = read.table("ASO_Capseq_R1.i32splicing.txt", header = T, row.names = 1, sep = "\t")
  mat2 = read.table("ASO_Capseq_R2.i32splicing.txt", header = T, row.names = 1, sep = "\t")
  colnames(mat1) = paste0(colnames(mat1), "_R1")
  colnames(mat2) = paste0(colnames(mat2), "_R2")
  
  mat1 = mat1[,rownames(grp1)]
  mat2 = mat2[,rownames(grp2)]
  mat = cbind(mat1, mat2)
  
  grp = rbind(grp1, grp2)
  grp0 = read.table("ASO_Capseq_R1_2.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  grp$sizeFactor = grp0[rownames(grp),]$sizeFactor
  
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glmm", estimate_offset = F, offset_str = "sizeFactor", apply_sva = F, test_terms = c("Treatment", "MIN"), remove_terms = NULL, random_terms = "MIN", distro_family = "nb")
  
  norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
  new1 = cbind(Exp = exp(unlist(norm1[1, ])), grp)
  ggplot(new1, aes(Treatment, Exp)) + geom_violin(width = 1.0) + geom_boxplot(outlier.alpha = 0, width = 0.1) + geom_noBG() + labs(x = "", y = "Normalized Expression (normal-scale)", color = "Donor")
  ggsave(paste0("pdf/boxplot_capseq_i32AS_ASO_",aso,"_effect_corrected.jpg"))
  
  tmp = get_fdr(res = test, term = "Treatment", class = aso, keyword = key_name)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(res = test, term = paste0("Treatment", aso), class = aso)[1,])
}

stats[order(stats$LFC.Treatment), ]


##################################################
#Comparisons of Treatment effects to untreated XDP
##################################################
df = rbind(bl,trt)
df$color = "Finite error"
df[is.na(df$margin),"color"] = "Infinite error clipped"
df[is.na(df$margin),"margin"] = 0.5
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON", "dSVA", as.vector(stats[order(stats$LFC.Treatment),"class"])))

ggplot(data = df, aes(x = class , y = exp(LFC), ymin = exp(LFC - margin), ymax = exp(LFC + margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(aes(color = color), size=0.5, width = 0.3) + geom_color(c("black","red")) +
  geom_hline(yintercept = exp(df[df$class=="dSVA","LFC"]), linetype="dashed", color = "orange") +
  geom_hline(yintercept = exp(df[df$class=="CON","LFC"]), linetype="dashed", color = "blue") +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0.0001, 1.5))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_capseq_FCi32AS_to_untreated_XDP_corrected.pdf")















#################################################
#Use untreated CON as baseline
#################################################
#################################################
#Building AS models for baseline: XDP vs CON
#################################################
#i32AS
key_name = "X_70644089_70659646_1"
grp = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
#grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$LibType != "ILL"), ]
grp = grp[which(grp$geno %in% c("XDP_non_edit", "XDP", "CON")), ]
grp$Genotype = factor(grp$Genotype, levels = c("CON", "XDP"))
grp$Genotype1 = grp$geno
#grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_unedited"
grp[grp$geno=="XDP", "Genotype1"] = "XDP_unexposed"
grp$Genotype1 = factor(grp$Genotype1, levels = c("CON", "XDP_unexposed", "XDP_unedited"))
grp$BatchGrowth = factor(grp$BatchGrowth)
table(paste0(grp$Genotype,grp$ParentLine))

mat0 = read.table("ASO_Capseq_R1.i32splicing.txt", header = T, row.names = 1, sep = "\t")
mat = mat0[,rownames(grp)]
norm = log1p(t(t(mat)/grp$sizeFactor))
#rnaPCA(norm,grp,c("Genotype"),ntop = 5)
new = cbind(Exp = unlist(norm[1, ]), grp)
write.table(new, file = "ASO_Capseq_R1.metadata_with_i32splicing_baseline.txt", row.names = T, col.names = T, quote = F, sep = "\t")
ggplot(new, aes(Genotype, log(exp(Exp)))) + geom_boxplot(outlier.alpha = 0) +  geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_raw_new.jpg")

ms = "AS"
mf = "Genotype"
pf = "baseline_vs_CON"
method1 = "glm"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method1, estimate_offset = F, offset_str = "sizeFactor", apply_sva = F, test_terms = c("Genotype", "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "nb")

write.table(test@vcov, paste0("model_params_vcov_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_", method1, "_", ms, "_", mf, "_", pf, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")




norm1 = test@deSV
#rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
new1 = cbind(Exp = unlist(norm[1, ]), grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_corrected_new.jpg")
ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_corrected_v2_new.jpg")
new2 = new1[new1$Genotype=="XDP"|(new1$Genotype=="CON" & new1$Exp == 0)|(new1$Genotype=="dSVA" & new1$Exp == 0),]
ggplot(new2, aes(Genotype, log2(exp(Exp)))) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
ggsave("pdf/boxplot_capseq_i32AS_XDP_vs_CON_corrected_outlierSamplesRemoved_new.jpg")

aa = get_fdr(res = test, term = "Genotype", class = "XDP_vs_CON", keyword = key_name)
write.table(aa, file = paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
bb = get_ci(res = test, term = "GenotypeXDP", class = "XDP_vs_CON")
write.table(bb, file = paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

#bl1 = get_ci(res = test, term = "GenotypedSVA", class = "dSVA_vs_CON")
#bl1 = bl1[1,]
bl = bb[1, ]

df_as = c()

as_cpm = exp(norm1)[1, ]
as_grp0 = grp[colnames(norm1), c("Genotype", "MIN", "Clone")]
df_as1 = data.frame(CPM = unname(as_cpm), as_grp0)
df_as = rbind(df_as, df_as1)

#####################################################
#Building AS models for treated XDP and untreated CON
#####################################################
#Treatment
key_name = "X_70644089_70659646_1"

sel = c("CONNO", "XDPNO", "dSVANO", "XDP126", "XDP131", "XDP134", "XDP307", "XDP849", "XDP852", "XDP876", "XDP877", "XDP879", "XDP880", "XDP881", "XDP883", "XDPMAL")
grp1 = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp1 = grp1[,c(-5, -6)]
grp1 = grp1[which(grp1$LibType != "ILL"), ]
grp1$Trt = paste0(grp1$Genotype, grp1$Treatment)
grp1 = grp1[which(grp1$Trt %in% sel),]
grp1$Trt = factor(grp1$Trt, levels = sel)
grp1$MIN = factor(grp1$MIN)

grp2 = read.table("ASO_Capseq_R2.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp2 = grp2[,c(-5, -6)]
grp2 = grp2[which(grp2$LibType != "ILL"), ]
grp2$Trt = paste0(grp2$Genotype, grp2$Treatment)
grp2 = grp2[which(grp2$Trt %in% sel),]
grp2$Trt = factor(grp2$Trt, levels = sel)
grp2$MIN = factor(grp2$MIN)

if (nrow(grp1) > 0){
  rownames(grp1) = paste0(rownames(grp1), "_R1")
}
if (nrow(grp2) > 0){
  rownames(grp2) = paste0(rownames(grp2), "_R2")
}

mat1 = read.table("ASO_Capseq_R1.i32splicing.txt", header = T, row.names = 1, sep = "\t")
mat2 = read.table("ASO_Capseq_R2.i32splicing.txt", header = T, row.names = 1, sep = "\t")
colnames(mat1) = paste0(colnames(mat1), "_R1")
colnames(mat2) = paste0(colnames(mat2), "_R2")

if (nrow(grp1) > 0){
  mat1 = mat1[,rownames(grp1)]
} else {
  mat1 = c()
}
if (nrow(grp2) > 0){
  mat2 = mat2[,rownames(grp2)]
} else {
  mat2 = c()
}
if (is.null(mat1) || is.null(mat2)){
  mat = rbind(mat1, mat2)
} else {
  mat = cbind(mat1, mat2)
}
write.table(mat, file = "ASO_Capseq_R1_2.i32splicing.final.4modelingTrt.txt", row.names = T, col.names = T, quote = F, sep = "\t")

grp = rbind(grp1, grp2)
grp0 = read.table("ASO_Capseq_R1_2.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp$sizeFactor = grp0[rownames(grp),]$sizeFactor
grp$Trt = gsub("NO$", "", as.vector(grp$Trt))
sel0 = gsub("NO$", "", sel)
grp$Trt = factor(grp$Trt, levels = sel0)
write.table(grp, "ASO_Capseq_R1_2.metadata.final.4modelingTrt.txt", row.names = T, col.names = T, quote = F, sep = "\t")


grp = read.table("ASO_Capseq_R1_2.metadata.final.4modelingTrt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = read.table("ASO_Capseq_R1_2.i32splicing.final.4modelingTrt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat[, rownames(grp)]








ms = "AS"
method1 = "glm"
mf = "Trt"
pf = "ref_CON"

test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method1, estimate_offset = F, offset_str = "sizeFactor", apply_sva = F, test_terms = c(mf), remove_terms = NULL, random_terms = "MIN", distro_family = "nb")

write.table(test@vcov, paste0("model_params_vcov_",method1,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method1,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")


  
norm2 = test@deSV

as_cpm = exp(norm2)[1, ]
as_grp0 = grp[colnames(norm2), c("Treated", "MIN", "Clone")]
sel = which(grp$Treated == paste0("XDP", aso))
df_as2 = data.frame(CPM = unname(as_cpm)[sel], as_grp0[sel, ])
colnames(df_as2)[2] = "Genotype"
df_as = rbind(df_as, df_as2)


#rnaPCA(norm2,grp,c("Genotype"),ntop = 5)
#new1 = cbind(Exp = exp(unlist(norm1[1, ])), grp)
#ggplot(new1, aes(Treated, Exp)) + geom_violin(width = 1.0) + geom_boxplot(outlier.alpha = 0, width = 0.1) + geom_noBG() + labs(x = "", y = "Normalized Expression (normal-scale)", color = "Donor")
#ggsave(paste0("pdf/boxplot_capseq_i32AS_ASO_",aso,"_effect_corrected.jpg"))

#tmp = get_fdr(res = test, term = "Treated", class = paste0("XDP", aso, "_vs_CON"), keyword = key_name)
#colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
#stats = rbind(stats, tmp)
#trt = rbind(trt, get_ci(res = test, term = paste0("TreatedXDP", aso), class = paste0("XDP", aso, "_vs_CON"))[1,])

  

  
trt = c()
stats = c()
for (v in c("dSVA", "XDP126", "XDP131", "XDP134", "XDP307", "XDP849", "XDP852", "XDP876", "XDP877", "XDP879", "XDP880", "XDP881", "XDP883", "XDPMAL")){
  mc = paste0(v, "_vs_CON")
  aa = get_fdr(test, term = paste0("Trt", v), class = mc, keyword = key_name)
  bb = get_ci(test, term = paste0("Trt", v), class = mc)
  trt = rbind(trt, bb[1,])
  write.table(aa, file = paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",mc,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(bb, file = paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",mc,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  colnames(aa) = c("FDR.Trt", "LFC.Trt", "pval.Trt", "name", "class")
  stats = rbind(stats, aa)
}





stats[order(stats$LFC.Treatment), ]


##################################################
#Comparisons of Treatment effects to untreated CON
##################################################
df = rbind(bl,trt)
df = as.data.frame(df)
df$class = gsub("_vs_CON", "", df$class)
df$class = factor(df$class,levels = c("XDP", gsub("_vs_CON", "", as.vector(stats[order(stats$LFC.Trt),"class"]))))

ggplot(data = df, aes(x = class , y = exp(LFC), ymin = exp(LFC - margin), ymax = exp(LFC + margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) + 
  #geom_color(c("black","red")) +
  geom_hline(yintercept = exp(df[df$class=="dSVA","LFC"]), linetype="dashed", color = "purple") +
  geom_hline(yintercept = exp(df[df$class=="XDP","LFC"]), linetype="dashed", color = "red") +
  geom_hline(yintercept = 1, linetype="dashed", color = "blue") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(0.0001, 1.5))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_capseq_FCi32AS_to_untreated_CON_corrected.pdf")

df22 = rbind(df[1, ], df[3:nrow(df), ][order(df[3:nrow(df), "LFC"]), ], df[2, ])
rownames(df22) = NULL
df22$trt = c("dSVA", gsub("XDP", "ASO", df22$class[2:(nrow(df22)-1)]), "XDP_untreated")
df22$trt = factor(df22$trt, levels = as.vector(df22$trt))
df22$Family = "Family134"
df22[which(df22$trt %in% c("ASO876", "ASO131", "ASO877", "ASO849")), "Family"] = "Family131"
df22[which(df22$trt %in% c("ASO307", "ASO126")), "Family"] = "Family126"
df22[which(df22$trt %in% c("ASOMAL")), "Family"] = "Negative control"
df22[which(df22$trt %in% c("dSVA")), "Family"] = "dSVA"
df22[which(df22$trt %in% c("XDP_untreated")), "Family"] = "Untreated"
df22$Family = factor(df22$Family, levels = c("dSVA", "Family134", "Family131", "Family126", "Negative control", "Untreated"))
ggplot(data = df22, aes(x = trt , y = exp(LFC), fill = Family)) +
  geom_bar(stat = "identity") +
  #geom_errorbar(size=0.5, width = 0.3) + 
  geom_fill(c("purple","red", "blue", "orange", "grey", "black")) +
  #geom_hline(yintercept = exp(df[df$class=="dSVA","LFC"]), linetype="dashed", color = "purple") +
  #geom_hline(yintercept = exp(df[df$class=="XDP","LFC"]), linetype="dashed", color = "red") +
  geom_hline(yintercept = 1, linetype="dashed", color = "brown") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(0.0001, 1.5))+
  xlab("") + ylab("Fold Change to Controls")

ggsave("pdf/R01_barplot_capseq_FCi32AS_to_untreated_CON_corrected.pdf")



df_as$Genotype = factor(df_as$Genotype, levels = c("CON", "XDP", "dSVA", "XDP880", "XDP879", "XDP307", "XDP134", "XDP852", "XDP876", "XDP126", "XDP131", "XDP881", "XDP877", "XDP849", "XDP883", "XDPMAL"))
ggplot(data = df_as, aes(x = Genotype , y = CPM, color = MIN)) +
  geom_quasirandom() +
  #geom_errorbar(size=0.5, width = 0.3) + 
  #geom_color(c("black","red")) +
  geom_hline(yintercept = mean(df_as[df_as$Genotype == "dSVA", "CPM"]), linetype="dashed", color = "purple") +
  geom_hline(yintercept = mean(df_as[df_as$Genotype == "XDP", "CPM"]), linetype="dashed", color = "red") +
  geom_hline(yintercept = mean(df_as[df_as$Genotype == "CON", "CPM"]), linetype="dashed", color = "blue") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(0.0001, 1.5))+
  xlab("") + ylab("i32-AS Normalized Expression")
ggsave("pdf/beeswarm_capseq_Expi32AS_to_untreated_CON_corrected_normal_scale.pdf")

ggplot(data = df_as, aes(x = Genotype , y = log10(CPM), color = MIN)) +
  geom_quasirandom() +
  #geom_errorbar(size=0.5, width = 0.3) + 
  #geom_color(c("black","red")) +
  #geom_hline(yintercept = mean(log10(df_as[df_as$Genotype == "dSVA", "CPM"])), linetype="dashed", color = "purple") +
  #geom_hline(yintercept = mean(log10(df_as[df_as$Genotype == "XDP", "CPM"])), linetype="dashed", color = "red") +
  #geom_hline(yintercept = mean(log10(df_as[df_as$Genotype == "CON", "CPM"])), linetype="dashed", color = "blue") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(0.0001, 1.5))+
  xlab("") + ylab("i32-AS Normalized Expression\n(log10-transformed)")
ggsave("pdf/beeswarm_capseq_Expi32AS_to_untreated_CON_corrected_log10_scale.pdf")

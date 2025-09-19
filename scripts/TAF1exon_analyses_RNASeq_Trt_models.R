library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
library(ggplot2)
library(ggforce)
library(stringr)
source("~/lib/R/countsFilter.R")
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/lib/R/goplot.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/glmm.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/go.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/plot.R")


setwd("~/projects/xdp_aso/nsc/output/")


####################################################################################################################################################################################
#CON_NO as reference level
####################################################################################################################################################################################
sel = c("CON", "XDP", "dSVA", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307", "CON131", "CON880", "CON881", "CON877", "CON879", "CON876", "CON307")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO", "XDP_880", "CON_880")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO")
#sel = c("CON_NO", "dSVA_NO", "XDP_880")
grp = read.table("model_params_colData_glm_Exp_Trt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
#grp$Trt = paste0(grp$Genotype, grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)

mat0 = read.table("ASO_RNAseq_R1.TAF1exons.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat0[, rownames(grp)]

ms = "TAF1exon"
method = "glm"
mf = "Trt"
#pf = "vs_CON_NO"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = "buildDEmodel.sizeFactor", apply_sva = F, test_terms = c(mf, colnames(grp)[grep("SV", colnames(grp))]), remove_terms = colnames(grp)[grep("SV", colnames(grp))], random_terms = NULL, distro_family = "nb")

for (x in sel[2:length(sel)]){
  pf = paste0(x, "_vs_CON")
  stat = get_fdr(test, term = paste0("Trt", x), class = paste0(x, "_vs_CON"))
  colnames(stat) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  ci = get_ci(test, term = paste0("Trt", x), class = paste0(x, "_vs_CON"), ci = 0.9)
  write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}


####################################################################################################################################################################################
#CON_880 as reference level
####################################################################################################################################################################################
sel = c("CON880", "XDPNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307", "CONNO", "CON131")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO", "XDP_880", "CON_880")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO")
#sel = c("CON_NO", "dSVA_NO", "XDP_880")
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp$Trt = paste0(grp$Genotype, grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)

mat0 = read.table("ASO_RNAseq_R1.TAF1exons.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat0[, rownames(grp)]

grp$Trt = gsub("NO$", "", as.vector(grp$Trt))
sel0 = gsub("NO$", "", sel)
grp$Trt = factor(grp$Trt, levels = sel0)

grp0 = read.table("model_params_colData_glm_Exp_Trt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp$sf = grp0[rownames(grp), "buildDEmodel.sizeFactor"]

ms = "TAF1exon"
method = "glm"
mf = "Trt"
#pf = "vs_CON_NO"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = "sf", apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")

for (x in sel0[2:length(sel0)]){
  pf = paste0(x, "_vs_CON880")
  stat = get_fdr(test, term = paste0("Trt", x), class = paste0(x, "_vs_CON880"))
  colnames(stat) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  ci = get_ci(test, term = paste0("Trt", x), class = paste0(x, "_vs_CON880"), ci = 0.9)
  write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}


####################################################################################################################################################################################
#CON_880 as reference level
####################################################################################################################################################################################
sel = c("CON131", "XDPNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307", "CONNO", "CON880")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO", "XDP_880", "CON_880")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO")
#sel = c("CON_NO", "dSVA_NO", "XDP_880")
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp$Trt = paste0(grp$Genotype, grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)

mat0 = read.table("ASO_RNAseq_R1.TAF1exons.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat0[, rownames(grp)]

grp$Trt = gsub("NO$", "", as.vector(grp$Trt))
sel0 = gsub("NO$", "", sel)
grp$Trt = factor(grp$Trt, levels = sel0)

grp0 = read.table("model_params_colData_glm_Exp_Trt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp$sf = grp0[rownames(grp), "buildDEmodel.sizeFactor"]

ms = "TAF1exon"
method = "glm"
mf = "Trt"
#pf = "vs_CON_NO"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = "sf", apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")

for (x in sel0[2:length(sel0)]){
  pf = paste0(x, "_vs_CON131")
  stat = get_fdr(test, term = paste0("Trt", x), class = paste0(x, "_vs_CON131"))
  colnames(stat) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  ci = get_ci(test, term = paste0("Trt", x), class = paste0(x, "_vs_CON131"), ci = 0.9)
  write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}











df = c()
tmp1 = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_","XDP_vs_CON",".txt"), header = T, sep = "\t", check.names = F)
tmp2 = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_","dSVA_vs_CON",".txt"), header = T, sep = "\t", check.names = F)
tmp3 = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_","XDP880_vs_CON",".txt"), header = T, sep = "\t", check.names = F)
tmp4 = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_","XDP131_vs_CON",".txt"), header = T, sep = "\t", check.names = F)
df = rbind(df, tmp1, tmp2, tmp3, tmp4)

df$class = sapply(strsplit(as.vector(df$class), "_"), "[[", 1)
df$class = factor(df$class, levels = c("XDP", "dSVA", "XDP880", "XDP131"))
df$name = gsub("TAF1:0", "", as.vector(df$name))

ggplot(data = df, aes(x = name , y = exp(LFC), group = class, color = class)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype="dashed", color = "black") +
  geom_vline(xintercept = "32", linetype="dashed", color = "black") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0, 2))+
  xlab("") + ylab("Fold Change") + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scale_color_manual(values = c("darkgrey", "orange", "blue","red"))

ggsave("pdf/curve_TAF1exon_XDP_dSVA_XDP880_XDP131.pdf", width = 7.1, height = 4.24)




ggplot(data = df, aes(x = name , y = exp(LFC), ymin = exp(LFC-margin), ymax = exp(LFC+margin), color = class)) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3, position = "dodge") +
  geom_hline(yintercept = 1, linetype="dashed", color = "black") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0, 2))+
  xlab("") + ylab("Fold Change") + theme(axis.text.x = element_text(angle = 45, hjust=1)) 
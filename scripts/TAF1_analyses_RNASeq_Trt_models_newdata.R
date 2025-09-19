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
#sel = c("CONNO", "XDPNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307", "CON131", "CON880", "CON881", "CON877", "CON879", "CON876", "CON307")
sel = c("CONNO", "XDPNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")

#sel = c("CON_NO", "XDP_NO", "dSVA_NO", "XDP_880", "CON_880")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO")
#sel = c("CON_NO", "dSVA_NO", "XDP_880")
mat0 = read.table("ASO_RNAseq_new.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp = read.table("ASO_RNAseq_new.569_namecorrected.metadata.tab", header = T, row.names = 1, sep = "\t")
grp = grp[colnames(mat0), ]
table(grp$Genotype)
grp$Trt = paste0(grp$Genotype, grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)

libsz = read.table("ASO_RNAseq_new.libSize.txt", header = F, check.names = F, row.names = 1)
grp$libsz = as.vector(libsz[rownames(grp), "V2"])

mat = mat0[, rownames(grp)]
CPM = t(1000000 * t(mat) / grp$libsz)
condition_list = lapply(sel, FUN = function(x){return(which(grp$Trt == x))})
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list)
length(which(expressedGenes))
mat = mat[which(expressedGenes), ]

grp$Trt = gsub("NO$", "", as.vector(grp$Trt))
sel0 = gsub("NO$", "", sel)
grp$Trt = factor(grp$Trt, levels = sel0)

ms = "Exp"
method = "glm"
mf = "Trt"
pf = "ref_CON"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")

#saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,".rds"))
write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

key_name = "ENSG00000147133"
test@beta[key_name, ]

norm = test@deSV
rnaPCA(norm,grp,c("Genotype"), textLabel = NA, ntop=100)
new = cbind(Exp = norm[key_name, ], grp)
ggplot(new, aes(Trt, Exp)) + geom_boxplot() + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)")

for (x in sel0[2:10]){
  pf = paste0(x, "_vs_CON")
  stat = get_fdr(beta = "model_params_beta_glm_Exp_Trt.txt", pval = "model_params_pval_glm_Exp_Trt.txt", fdr = "model_params_fdr_glm_Exp_Trt.txt", term = paste0("Trt", x, "$"), class = paste0(x, "_vs_CON"))
  colnames(stat) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  print(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"))
  
  ci = get_ci(beta = "model_params_beta_glm_Exp_Trt.txt", stderr = "model_params_stderr_glm_Exp_Trt.txt", pval = "model_params_pval_glm_Exp_Trt.txt", fdr = "model_params_fdr_glm_Exp_Trt.txt", term = paste0("Trt", x), class = paste0(x, "_vs_CON"))
  write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  print(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"))
}

x = "XDP"
pf = paste0(x, "_vs_CON")
stat = get_fdr(beta = "model_params_beta_glm_Exp_Trt.txt",
               pval = "model_params_pval_glm_Exp_Trt.txt",
               fdr = "model_params_fdr_glm_Exp_Trt.txt",
               term = paste0("Trt", x,"$"), class = paste0(x, "_vs_CON"))
colnames(stat) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

ci = get_ci(beta = "model_params_beta_glm_Exp_Trt.txt",
            stderr = "model_params_stderr_glm_Exp_Trt.txt",
            pval = "model_params_pval_glm_Exp_Trt.txt",
            fdr = "model_params_fdr_glm_Exp_Trt.txt",
            term = paste0("Trt", x), class = paste0(x, "_vs_CON"))
write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")



bl = c()
for (x in c("XDP", "dSVA")){
  #tmp = read.table(paste0("res_LFC_CI_glm_Exp_Genotype_", x, "_vs_CON.txt"), header = T)
  tmp = read.table(paste0("res_LFC_CI_glm_Exp_Trt_", x, "_vs_CON.txt"), header = T)
  tmp = tmp[which(tmp$name == key_name), ]
  bl = rbind(bl, tmp)
}
rownames(bl) = NULL

trt = c()
for (x in sel[4:10]){
  tmp = read.table(paste0("res_LFC_CI_glm_Exp_Trt_", x, "_vs_CON.txt"), header = T)
  tmp = tmp[which(tmp$name == key_name), ]
  trt = rbind(trt, tmp)
}
rownames(trt) = NULL

df = rbind(bl, trt)
new_class = gsub("_vs_CON", "", as.vector(df$class))
new_class[3:length(new_class)] = gsub("XDP", "", new_class[3:length(new_class)])
df$class = new_class
df$class = factor(df$class, levels = new_class)

ggplot(data = df, aes(x = class , y = exp(LFC), ymin = exp(LFC-margin), ymax = exp(LFC+margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = exp(df[df$class=="dSVA","LFC"]), linetype="dashed", color = "purple") +
  geom_hline(yintercept = exp(df[df$class=="XDP","LFC"]), linetype="dashed", color = "red") +
  geom_hline(yintercept = 1, linetype="dashed", color = "black") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0.65, 1.15), breaks = seq(0.65, 1.15, 0.05))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_rnaseq_FC_to_untreated_XDP_corrected_new_Trt.pdf", width = 5.7, height = 4.2)














####################################################################################################################################################################################
#XDP as reference level
####################################################################################################################################################################################
#sel = c("XDPNO", "CONNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307", "CON131", "CON880", "CON881", "CON877", "CON879", "CON876", "CON307")
sel = c("XDPNO", "CONNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")

mat0 = read.table("ASO_RNAseq_new.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp = read.table("ASO_RNAseq_new.569_namecorrected.metadata.tab", header = T, row.names = 1, sep = "\t")
grp = grp[colnames(mat0), ]
table(grp$Genotype)
grp$Trt = paste0(grp$Genotype, grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)


libsz = read.table("ASO_RNAseq_new.libSize.txt", header = F, check.names = F, row.names = 1)
grp$libsz = as.vector(libsz[rownames(grp), "V2"])

mat = mat0[, rownames(grp)]
CPM = t(1000000 * t(mat) / grp$libsz)
condition_list = lapply(sel, FUN = function(x){return(which(grp$Trt == x))})
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list)
length(which(expressedGenes))
mat = mat[which(expressedGenes), ]

grp$Trt = gsub("NO$", "", as.vector(grp$Trt))
sel0 = gsub("NO$", "", sel)
grp$Trt = factor(grp$Trt, levels = sel0)

ms = "Exp"
method = "glm"
mf = "Trt"
pf = "ref_XDP"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")

#saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,".rds"))
write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

key_name = "ENSG00000147133"
test@beta[key_name, ]

norm = test@deSV
rnaPCA(norm,grp,c("Trt"), textLabel = NA, ntop=100)
new = cbind(Exp = norm[key_name, ], grp)
ggplot(new, aes(Trt, Exp)) + geom_boxplot() + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)")

for (x in sel0[2:length(sel0)]){
  pf = paste0(x, "_vs_XDP")
  stat = get_fdr(beta = "model_params_beta_glm_Exp_Trt_ref_XDP.txt",
                 pval = "model_params_pval_glm_Exp_Trt_ref_XDP.txt",
                 fdr = "model_params_fdr_glm_Exp_Trt_ref_XDP.txt",
                 term = paste0("Trt", x), class = paste0(x, "_vs_XDP"))
  colnames(stat) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  print(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"))
  
  ci = get_ci(beta = "model_params_beta_glm_Exp_Trt_ref_XDP.txt",
              stderr = "model_params_stderr_glm_Exp_Trt_ref_XDP.txt",
              pval = "model_params_pval_glm_Exp_Trt_ref_XDP.txt",
              fdr = "model_params_fdr_glm_Exp_Trt_ref_XDP.txt",
              term = paste0("Trt", x), class = paste0(x, "_vs_XDP"))
  write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  print(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"))
}


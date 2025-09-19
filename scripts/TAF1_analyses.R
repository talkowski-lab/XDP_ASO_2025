library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
source("~/lib/R/countsFilter.R")
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")


#############################################
#Use untreated XDP as baseline
#############################################
#############################################
#Build Expression models for CON/dSVA vs XDP
#############################################
#TAF1 expression
key_name = "ENSG00000147133"
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
#grp = grp[which(grp$Treatment %in% c("NO", "131", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$Genotype %in% c("XDP", "dSVA")), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "dSVA"))
#grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
#grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", "CON", "dSVA"))
grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", "dSVA"))
#grp$Treatment = factor(grp$Treatment, levels = c("NO", "307", "876", "877", "879", "880", "881"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$ParentLine = factor(grp$ParentLine)
grp$MIN = factor(grp$MIN)
grp$MIN1 = as.vector(grp$MIN)
grp[grp$MIN1=="35326", "MIN1"] = "35613"
grp$MIN1 = factor(grp$MIN1)
grp$MIN2 = as.vector(grp$MIN)
grp[grp$Genotype=="CON", "MIN2"] = "0"
grp$MIN2 = factor(grp$MIN2)

mat0 = read.table("ASO_RNAseq_R1.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat0[,rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
#condition_list = list(which(grp$Genotype1=="CON"), which(grp$Genotype1=="XDP"), which(grp$Genotype1=="XDP_noEdit"), which(grp$Genotype1=="dSVA"))
condition_list = list(which(grp$Genotype1=="dSVA"), which(grp$Genotype1=="XDP"), which(grp$Genotype1=="XDP_noEdit"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
mat = mat[which(expressedGenes), ]
nrow(mat) #17688

#norm = log(CPM)
#rnaPCA(norm,grp,c("Genotype"), textLabel = "ParentLine", ntop=100)
#ggsave("pdf/PCA_rnaseq_TAF1_XDP_vs_CON_raw_new.jpg")
#new = cbind(Exp = norm[key_name, ], grp)
#ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
#ggsave("pdf/boxplot_rnaseq_TAF1_XDP_vs_CON_raw_new.jpg")


ms = "Exp"
method = "glmm"
mf = "Genotype"
date()
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth", "MIN"), remove_terms = "BatchGrowth", random_terms = c("MIN"), distro_family = "nb")
date()

write.table(test@colData, paste0("RData/models_",method,"_",ms,"_",mf,".rds"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt = get_fdr(test,"Genotype","test")
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl1 = get_ci(test,"GenotypeCON", "CON")
bl2 = get_ci(test,"GenotypedSVA", "dSVA")
bl = as.data.frame(rbind(bl1,bl2))
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

#norm1 = test@deSV
#rnaPCA(norm1,grp,c("Genotype"), textLabel = "MIN", ntop=500)
#ggsave("pdf/PCA_rnaseq_TAF1_XDP_vs_CON_corrected_new.jpg")
#new1 = cbind(Exp = norm1["ENSG00000147133", ], grp)
#ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
#ggsave("pdf/boxplot_rnaseq_TAF1_XDP_vs_CON_corrected_new.jpg")
#ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
#ggsave("pdf/boxplot_rnaseq_TAF1_XDP_vs_CON_corrected_v2_new.jpg")


ms = "Exp"
method = "glmm"
mf = "Genotype1"
test1 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_model_str = "(1|ParentLine) + (1|BatchGrowth)", offset_str = "sizeFactor")

write.table(test1@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt00 = get_fdr(test1,"Genotype1","test")
write.table(tt00, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl01 = get_ci(test1,"Genotype1XDP_noEdit", "XDP_noEdit")
bl02 = get_ci(test1,"Genotype1CON", "CON")
bl03 = get_ci(test1,"Genotype1dSVA", "dSVA")
bl00 = as.data.frame(rbind(bl01,bl02,bl03))
write.table(bl00, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test, test1)
gc(verbose = F)


######################################################
#Build Expression models for treated and untreated XDP
######################################################
#Treatment
mat0 = read.table("ASO_RNAseq_R1.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)

trt = c()
stats = c() 
#aso = "307"
for (aso in c("131", "307", "876", "877", "879", "880", "881")){
  message(paste0("ASO",aso))
  grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  #grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
  grp = grp[which(grp$Treatment %in% c("NO", aso)), ]
  grp = grp[which(grp$Genotype %in% c("XDP")), ]
  #grp$Treatment = factor(grp$Treatment,levels=c("NO", "307", "876", "877", "879", "880", "881"))
  grp$Treatment = factor(grp$Treatment,levels=c("NO", aso))
  grp$BatchGrowth = factor(grp$BatchGrowth)
  grp$BatchTreat = factor(grp$BatchTreat)
  grp$MIN = factor(grp$MIN)
  grp$ParentLine = factor(grp$ParentLine)
  #table(paste0(grp$Treatment,"_",grp$ParentLine))
  
  #mat = round(mat /10)
  mat = mat0[,rownames(grp)]
  CPM = t(t(mat) / grp$sizeFactor)
  condition_list = list(which(grp$Treatment=="NO"), which(grp$Treatment==aso))
  expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
  mat = mat[which(expressedGenes), ]
  message(nrow(mat)," genes expressed.")
  
  #norm = log1p(t(t(mat)/as.vector(grp$sizeFactor)))
  #rnaPCA(norm,grp,c("Treatment"), ntop = 100)
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_raw.jpg"))
  #new = cbind(Exp = norm["ENSG00000147133", ], grp)
  #ggplot(new, aes(Treatment, log2(exp(Exp)))) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_ASO_",aso,"_effect_raw.jpg"))
  
  ms = "Exp"
  method = "glmm"
  mf = "Treatment"
  #plan("multisession", workers = 4)
  #date()
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, apply_sva = T, test_terms = c(mf, "MIN"), remove_terms = NULL, random_terms = "MIN", offset_str = "sizeFactor")
  #date()
  write.table(test@colData, paste0("RData/models_",method,"_",ms,"_",mf,"_",aso,".rds"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  tmp = get_fdr(test, "Treatment", aso, NULL)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(test,paste0("Treatment",aso), aso))
  
  #norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Treatment"))
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_corrected.jpg"))
  #new1 = cbind(Exp = norm1["ENSG00000147133", ], grp)
  #ggplot(new1, aes(Treatment, log2(exp(Exp)))) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_ASO_",aso,"_effect_corrected.jpg"))
  
  rm(test, tmp)
  gc(verbose = F)
}

write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


######################################################
#Comparisons of Treatment effects to untreated XDP
######################################################
key_name = "ENSG00000147133"
ms = "Exp"
method = "glmm"
mf = "Treatment"

stats = read.table(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), header = T, sep = "\t")
stats = stats[stats$name == key_name, ]
stats[order(stats[,2], decreasing = T), ]

trt = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), header = T, sep = "\t")
trt = trt[trt$name == key_name, ]

method = "glm"
mf = "Genotype"
bl = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"),header = T, sep = "\t")
bl = bl[bl$name == key_name, ]

df = rbind(bl,trt)
df[df==Inf] = 0.05
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON","dSVA", "131", "880", "881", "877", "879", "876", "307"))

ggplot(data = df, aes(x = class , y = exp(LFC), ymin = exp(LFC-margin), ymax = exp(LFC+margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = exp(df[df$class=="dSVA","LFC"]), linetype="dashed", color = "orange") +
  geom_hline(yintercept = exp(df[df$class=="CON","LFC"]), linetype="dashed", color = "blue") +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(-0.013, 0.005))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_rnaseq_FC_to_untreated_XDP_corrected_new.jpg")








######################################################
#Comparisons of Treatment effects to untreated CON
######################################################
key_name = "ENSG00000147133"





sel = c("CON_NO", "XDP_NO", "dSVA_NO", "XDP_131", "XDP_880", "XDP_881", "XDP_877", "XDP_879", "XDP_876", "XDP_307")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO", "XDP_880", "CON_880")
#sel = c("CON_NO", "XDP_NO", "dSVA_NO")
#sel = c("CON_NO", "dSVA_NO", "XDP_880")
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp$Trt = paste0(grp$Genotype, "_", grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)

mat = mat0[, rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
condition_list = lapply(sel, FUN = function(x){return(which(grp$Trt == x))})
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
length(which(expressedGenes))
mat = mat[which(expressedGenes), ]

test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glm", estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c("Trt"), remove_terms = NULL, random_terms = NULL, distro_family = "nb")

bl = c()
for (x in c("XDP", "dSVA")){
  tmp = read.table(paste0("res_LFC_CI_glm_Exp_Genotype_", x, "_vs_CON.txt"), header = T)
  tmp = tmp[which(tmp$name == key_name), ]
  bl = rbind(bl, tmp)
}
rownames(bl) = NULL

trt = c()
for (x in sel[4:length(sel)]){
  tmp = get_ci(test, term = paste0("Trt", x), class = x)
  tmp = tmp[which(tmp$name == key_name), ]
  trt = rbind(trt, tmp)
}
rownames(trt) = NULL

df = rbind(bl, trt)
new_class = c("XDP", "dSVA", gsub("XDP_", "", df$class[3:length(df$class)]))
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
  scale_y_continuous(limits = c(0.7, 1.15), breaks = seq(0.7, 1.15, 0.05))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_rnaseq_FC_to_untreated_XDP_corrected_new.pdf", width = 5.7, height = 4.2)

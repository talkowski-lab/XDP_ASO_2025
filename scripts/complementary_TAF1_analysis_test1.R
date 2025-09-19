library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/countsFilter.R")
source("~/lib/R/geom_noBG.R")
source("~/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")


#############################################
#RNASeq
#############################################
#TAF1 expression
key_name = "ENSG00000147133"
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
#grp = grp[which(grp$Treatment %in% c("NO", "131", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
#grp = grp[which(grp$Genotype %in% c("XDP", "dSVA", "CON")), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", "CON", "dSVA"))
#grp$Treatment = factor(grp$Treatment, levels = c("NO", "307", "876", "877", "879", "880", "881"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)

mat0 = read.table("ASO_RNAseq_R1.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat0[,rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
condition_list = list(which(grp$Genotype1=="CON"), which(grp$Genotype1=="XDP"), which(grp$Genotype1=="XDP_noEdit"), which(grp$Genotype1=="dSVA"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
mat = mat[which(expressedGenes), ]
#nrow(mat) #17688

#norm = log(CPM)
#rnaPCA(norm,grp,c("Genotype"), textLabel = "ParentLine", ntop=100)
#ggsave("pdf/PCA_rnaseq_TAF1_XDP_vs_CON_raw.jpg")
#new = cbind(Exp = norm[key_name, ], grp)
#ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
#ggsave("pdf/boxplot_rnaseq_TAF1_XDP_vs_CON_raw.jpg")


ms = "Exp"
method = "glmm"
mf = "Genotype"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_model_str = "(1|ParentLine)", offset_str = "sizeFactor")

#write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model1_params_warnings_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#tt = get_fdr(test,"Genotype","test")
#write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

#bl1 = get_ci(test,"GenotypeCON", "CON")
#bl2 = get_ci(test,"GenotypedSVA", "dSVA")
#bl = as.data.frame(rbind(bl1,bl2))
#write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

#norm1 = test@deSV
#rnaPCA(norm1,grp,c("Genotype"), textLabel = "ParentLine", ntop=500)
#ggsave("pdf/PCA_rnaseq_TAF1_XDP_vs_CON_corrected.jpg")
#new1 = cbind(Exp = norm1["ENSG00000147133", ], grp)
#ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
#ggsave("pdf/boxplot_rnaseq_TAF1_XDP_vs_CON_corrected.jpg")
#ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Normalized Expression (log-scale)", color = "Donor")
#ggsave("pdf/boxplot_rnaseq_TAF1_XDP_vs_CON_corrected_v2.jpg")


ms = "Exp"
method = "glmm"
mf = "Genotype1"
test1 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_model_str = "(1|ParentLine)", offset_str = "sizeFactor")

#write.table(test1@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test1@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test1@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test1@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(test1@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@status, paste0("model1_params_warnings_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#tt00 = get_fdr(test1,"Genotype1","test")
#write.table(tt00, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

#bl01 = get_ci(test1,"Genotype1XDP_noEdit", "XDP_noEdit")
#bl02 = get_ci(test1,"Genotype1CON", "CON")
#bl03 = get_ci(test1,"Genotype1dSVA", "dSVA")
#bl00 = as.data.frame(rbind(bl01,bl02,bl03))
#write.table(bl00, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

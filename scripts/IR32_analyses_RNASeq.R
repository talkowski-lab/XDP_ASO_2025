library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")


mat0 = read.table("ASO_RNAseq_R1.IR.txt", header = T, row.names = 1, sep = "\t", check.names = F)

#################################################
#Building IR models for baseline: CON/dSVA vs XDP
#################################################
#IR
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
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


mat = mat0[,rownames(grp)]
mat = 1000 * mat / get_intron_length(rownames(mat))

norm = mat
rnaPCA(norm,grp,c("Genotype"),ntop = 500)
ggsave("pdf/PCA_rnaseq_IR_XDP_vs_CON_raw_new.jpg")
new = cbind(Exp = unlist(norm[key_name, ]), grp)
ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_raw_new.jpg")


ms = "IR"
method = "glm"
mf = "Genotype"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = NULL, apply_sva = F, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "gaussian")

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

tt = get_fdr(test,"Genotype","test")
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl1 = get_ci(test,"GenotypeCON", "CON")
bl2 = get_ci(test,"GenotypedSVA", "dSVA")
bl = as.data.frame(rbind(bl1,bl2))
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

norm1 = test@deSV
rnaPCA(norm1,grp,c("Genotype"),ntop = 500)
ggsave("pdf/PCA_rnaseq_IR_XDP_vs_CON_corrected_new.jpg")
new1 = cbind(Exp = norm1[key_name, ], grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_corrected_new.jpg")
ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG()+ labs(x = "", y = "IR Ratio", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_corrected_v2_new.jpg")


ms = "IR"
method = "glm"
mf = "Genotype1"
test1 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = NULL, apply_sva = F, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "gaussian")

#saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,".rds"))
write.table(test1@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test1@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
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


#################################################
#Building IR models for baseline: XDP/dSVA vs CON
#################################################
#IR
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
#grp = grp[which(grp$Treatment %in% c("NO", "131", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
#grp = grp[which(grp$Genotype %in% c("XDP", "dSVA", "CON")), ]
grp$Genotype = factor(grp$Genotype, levels = c("CON", "XDP", "dSVA"))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp[grp$geno=="XDP", "Genotype1"] = "XDP_vanilla"
grp$Genotype1 = factor(grp$Genotype1, levels = c("CON", "XDP_vanilla", "XDP_noEdit", "dSVA"))
#grp$Treatment = factor(grp$Treatment, levels = c("NO", "307", "876", "877", "879", "880", "881"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)

mat = mat0[,rownames(grp)]
mat = 1000 * mat / get_intron_length(rownames(mat))

norm = mat
rnaPCA(norm,grp,c("Genotype"),ntop = 500)
ggsave("pdf/PCA_rnaseq_IR_XDP_vs_CON_raw_new.jpg")
new = cbind(Exp = unlist(norm[key_name, ]), grp)
ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_raw_new.jpg")


ms = "IR"
method = "glm"
mf = "Genotype"
pf = "baseline_vs_CON"
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = NULL, apply_sva = F, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "gaussian")

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

tt = get_fdr(test, term = "Genotype",class = "baseline_vs_CON")
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

tt1 = tt[, c(1,3,5,7,8)]
tt1$class = "XDP_vs_CON"
write.table(tt1, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_","XDP_vs_CON",".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

tt2 = tt[, c(2,4,6,7,8)]
tt2$class = "dSVA_vs_CON"
write.table(tt2, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_","dSVA_vs_CON",".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl1 = get_ci(test, term = "GenotypeXDP", class = "XDP_vs_CON")
bl2 = get_ci(test, term = "GenotypedSVA", class = "dSVA_vs_CON")
bl = as.data.frame(rbind(bl1,bl2))
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

norm1 = test@deSV
rnaPCA(norm1,grp,c("Genotype"),ntop = 500)
ggsave("pdf/PCA_rnaseq_IR_XDP_vs_CON_corrected_new.jpg")
new1 = cbind(Exp = norm1[key_name, ], grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_corrected_new.jpg")
ggplot(new1, aes(Genotype1, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG()+ labs(x = "", y = "IR Ratio", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_corrected_v2_new.jpg")

rm(test)
gc(verbose = F)


################################################
#Building IR models for treated vs untreated XDP
################################################
#mat0 = read.table("ASO_RNAseq_R1.IR.txt", header = T, row.names = 1, sep = "\t", check.names = F)
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"

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
  #table(paste0(grp$Treatment,"_",grp$ParentLine))
  
  mat = mat0[,rownames(grp)]
  mat = 1000 * mat / get_intron_length(rownames(mat))
  
  #norm = mat
  #rnaPCA(norm,grp,c("Treatment"), ntop = 100)
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_raw_new.jpg"))
  #new = cbind(Exp = unlist(norm[key_name, ]), grp)
  #ggplot(new, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_IR_ASO_",aso,"_effect_raw_new.jpg"))
  
  ms = "IR"
  method = "glm"
  #method = "glmm"
  mf = "Treatment"
  pf = paste0("ASO", aso, "_vs_XDP")
  #plan("multisession", workers = 4)
  #date()
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = NULL, apply_sva = F, test_terms = c(mf, "BatchTreat"), remove_terms = "BatchTreat", random_terms = NULL, distro_family = "gaussian")
  #test = buildDEmodel(mat, grp, method, F, c(mf), NULL, "(1|ParentLine) + (1|BatchTreat)", offset_str = NULL, distro_family = "gaussian")
  #date()
  #plan("sequential")
  #saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",aso,".rds"))
  write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  tmp = get_fdr(test, term = "Treatment", class = pf)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(test, term = paste0("Treatment", aso), class = pf))
  
  #norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Treatment"),ntop = 100)
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_corrected_new.jpg"))
  #new1 = cbind(Exp = unlist(norm1[key_name, ]), grp)
  #ggplot(new1, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_IR_ASO_",aso,"_effect_corrected_new.jpg"))
  
  rm(test, tmp)
  gc(verbose = F)
}

write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


####################################################
#Building IR models for treated XDP vs untreated CON
####################################################
#mat0 = read.table("ASO_RNAseq_R1.IR.txt", header = T, row.names = 1, sep = "\t", check.names = F)
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"

#aso = "307"
for (aso in c("131", "307", "876", "877", "879", "880", "881")){
  trt = c()
  stats = c()
  message(paste0("ASO",aso))
  grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  #grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
  grp = grp[which(grp$Treatment %in% c("NO", aso)), ]
  grp = grp[which(grp$Genotype %in% c("XDP", "CON")), ]
  #grp$Treatment = factor(grp$Treatment,levels=c("NO", "307", "876", "877", "879", "880", "881"))
  grp$Treatment = factor(grp$Treatment,levels=c("NO", aso))
  grp$BatchGrowth = factor(grp$BatchGrowth)
  grp$BatchTreat = factor(grp$BatchTreat)
  grp$MIN = factor(grp$MIN)
  #table(paste0(grp$Treatment,"_",grp$ParentLine))
  grp$ParentLine = factor(grp$ParentLine)
  grp$Treated = paste0(grp$Genotype, "_", grp$Treatment)
  
  grp = grp[which((grp$Treated == "CON_NO") | (grp$Treated == paste0("XDP_", aso))), ]
  grp[grp$Treated == "CON_NO", "Treated"] = "CON"
  grp[grp$Treated == paste0("XDP_", aso), "Treated"] = paste0("XDP", aso)
  grp$Treated = factor(grp$Treated, levels = c("CON", paste0("XDP", aso)))
  
  mat = mat0[,rownames(grp)]
  mat = 1000 * mat / get_intron_length(rownames(mat))
  
  #norm = mat
  #rnaPCA(norm,grp,c("Treatment"), ntop = 100)
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_raw_new.jpg"))
  #new = cbind(Exp = unlist(norm[key_name, ]), grp)
  #ggplot(new, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_IR_ASO_",aso,"_effect_raw_new.jpg"))
  
  ms = "IR"
  method = "glm"
  #method = "glmm"
  mf = "Treated"
  pf = paste0("XDP", aso, "_vs_CON")
  #plan("multisession", workers = 4)
  #date()
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = F, offset_str = NULL, apply_sva = F, test_terms = c(mf, "BatchTreat"), remove_terms = "BatchTreat", random_terms = NULL, distro_family = "gaussian")
  #test = buildDEmodel(mat, grp, method, F, c(mf), NULL, "(1|ParentLine) + (1|BatchTreat)", offset_str = NULL, distro_family = "gaussian")
  #date()
  #plan("sequential")
  #saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",pf,".rds"))
  write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  tmp = get_fdr(test, term = "Treated", class = pf)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(test, term = paste0("TreatedXDP", aso), class = pf))
  
  #norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Treatment"),ntop = 100)
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_corrected_new.jpg"))
  #new1 = cbind(Exp = unlist(norm1[key_name, ]), grp)
  #ggplot(new1, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_IR_ASO_",aso,"_effect_corrected_new.jpg"))
  
  rm(test, tmp)
  gc(verbose = F)
  
  write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}


##################################################
#Comparisons of Treatment effects to untreated XDP
##################################################
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
ms = "IR"
method = "glm"
mf = "Treatment"

stats = read.table(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), header = T, sep = "\t")
stats = stats[stats$name == key_name, ]
stats[order(stats[,2]),]

trt = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), header = T, sep = "\t")
trt = trt[trt$name == key_name, ]

mf = "Genotype"
bl = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"),header = T, sep = "\t")
bl = bl[bl$name == key_name, ]

df = rbind(bl,trt)
df$color = "Finite error"
df[is.infinite(df$margin),"color"] = "Infinite error clipped"
df[is.infinite(df$margin),"margin"] = 1.1 * max(df[!is.infinite(df$margin),"margin"])
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON","dSVA", "880", "131", "307", "877", "876", "881", "879"))

ggplot(data = df, aes(x = class , y = LFC, ymin = LFC-margin, ymax = LFC+margin)) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) + 
  #geom_color(c("black","red")) +
  geom_hline(yintercept = df[df$class=="dSVA","LFC"], linetype="dashed", color = "orange") +
  geom_hline(yintercept = df[df$class=="CON","LFC"], linetype="dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(-0.013, 0.005))+
  xlab("") + ylab("IR Change")

ggsave("pdf/errorbar_rnaseq_IRchange_to_untreated_XDP_corrected_new.pdf")




##################################################
#Comparisons of Treatment effects to untreated CON
##################################################
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
ms = "IR"
method = "glm"
mf = "Treated"

stats = c()
for (aso in c("131", "307", "876", "877", "879", "880", "881")){
  pf = paste0("XDP", aso, "_vs_CON")
  tmp = read.table(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, sep = "\t")
  tmp = tmp[tmp$name == key_name, ]
  stats = rbind(stats, tmp)
}
stats[order(stats[,2]),]

trt = c()
for (aso in c("131", "307", "876", "877", "879", "880", "881")){
  pf = paste0("XDP", aso, "_vs_CON")
  tmp = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, sep = "\t")
  tmp = tmp[tmp$name == key_name, ]
  trt = rbind(trt, tmp)
}

mf = "Genotype"
pf = "baseline_vs_CON"
bl = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"),header = T, sep = "\t")
bl = bl[bl$name == key_name, ]

df = rbind(bl,trt)
df = as.data.frame(df)
df$class = gsub("_vs_CON", "", df$class)
df$class = factor(df$class,levels = c("XDP","dSVA", "XDP880", "XDP131", "XDP307", "XDP877", "XDP876", "XDP881", "XDP879"))

ggplot(data = df, aes(x = class , y = LFC, ymin = LFC-margin, ymax = LFC+margin)) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) + 
  #geom_color(c("black","red")) +
  geom_hline(yintercept = df[df$class=="dSVA","LFC"], linetype="dashed", color = "purple") +
  geom_hline(yintercept = df[df$class=="XDP","LFC"], linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype="dashed", color = "blue") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(-0.013, 0.005))+
  xlab("") + ylab("IR Change per Kb")

ggsave("pdf/errorbar_rnaseq_IRchange_to_untreated_CON_corrected_new.pdf")


df22 = rbind(df[2, ], df[3:nrow(df), ][order(df[3:nrow(df), "LFC"]), ], df[1, ])
rownames(df22) = NULL
df22$trt = c("dSVA", gsub("XDP", "ASO", df22$class[2:(nrow(df22)-1)]), "XDP_untreated")
df22$trt = factor(df22$trt, levels = as.vector(df22$trt))
df22$Family = "Family134"
df22[which(df22$trt %in% c("ASO876", "ASO131", "ASO877", "ASO849")), "Family"] = "Family131"
df22[which(df22$trt %in% c("ASO307", "ASO126")), "Family"] = "Family126"
#df22[which(df22$trt %in% c("ASOMAL")), "Family"] = "Negative control"
df22[which(df22$trt %in% c("dSVA")), "Family"] = "dSVA"
df22[which(df22$trt %in% c("XDP_untreated")), "Family"] = "Untreated"
df22$Family = factor(df22$Family, levels = c("dSVA", "Family134", "Family131", "Family126", "Untreated"))
ggplot(data = df22, aes(x = trt , y = LFC, fill = Family)) +
  geom_bar(stat = "identity") +
  #geom_errorbar(size=0.5, width = 0.3) + 
  geom_fill(c("purple","red", "blue", "orange", "black")) +
  #geom_hline(yintercept = exp(df[df$class=="dSVA","LFC"]), linetype="dashed", color = "purple") +
  #geom_hline(yintercept = exp(df[df$class=="XDP","LFC"]), linetype="dashed", color = "red") +
  geom_hline(yintercept = 0, linetype="dashed", color = "brown") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(0.0001, 1.5))+
  xlab("") + ylab("IR i32 Change to Control")

ggsave(paste0("pdf/R01_barplot_rnaseq_IRchange_to_untreated_XDP_corrected_new_Trt_R1.pdf"))

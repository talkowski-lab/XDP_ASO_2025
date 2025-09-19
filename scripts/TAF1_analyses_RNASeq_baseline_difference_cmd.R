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
#RNASeq
#############################################
key_name = "ENSG00000147133"
mat0 = read.table("ASO_RNAseq_R1.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)


#Untreated dSVA vs XDP
lvl = "dSVA"
message(paste0("Untreated ", lvl, " vs XDP"))
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$Genotype %in% c("XDP", lvl)), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", lvl))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", lvl))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$ParentLine = factor(grp$ParentLine)
grp$MIN = factor(grp$MIN)

mat = mat0[,rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
condition_list = list(which(grp$Genotype1 == lvl), which(grp$Genotype1 == "XDP"), which(grp$Genotype1 == "XDP_noEdit"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
mat = mat[which(expressedGenes), ]
message(nrow(mat)," genes expressed.")

ms = "Exp"
method = "glmm"
mf = "Genotype"
pf = paste0(lvl, "_vs_XDP")
date()
plan(multicore, workers = 18)
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth", "MIN"), remove_terms = "BatchGrowth", random_terms = c("MIN"), distro_family = "nb")
plan(sequential)
date()

saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",pf,".rds"))
write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt = get_fdr(res = test, term = "Genotype", class = "test")
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl = get_ci(res = test, term = paste0("Genotype", lvl), class = lvl)
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test)
gc(verbose = F)


#Untreated CON vs XDP
lvl = "CON"
message(paste0("Untreated ", lvl, " vs XDP"))
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$Genotype %in% c("XDP", lvl)), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", lvl))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", lvl))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$ParentLine = factor(grp$ParentLine)
grp$MIN = factor(grp$MIN)

mat = mat0[,rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
condition_list = list(which(grp$Genotype1 == lvl), which(grp$Genotype1 == "XDP"), which(grp$Genotype1 == "XDP_noEdit"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
mat = mat[which(expressedGenes), ]
message(nrow(mat)," genes expressed.")

ms = "Exp"
method = "glm"
mf = "Genotype"
pf = paste0(lvl, "_vs_XDP")
date()
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "nb")
date()

saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",pf,".rds"))
write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#tt = get_fdr(beta = paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), pval = paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), fdr = paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), term = "Genotype", class = "test")
tt = get_fdr(res = test, term = "Genotype", class = "test")
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl = get_ci(res = test, term = paste0("Genotype", lvl), class = lvl)
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test)
gc(verbose = F)

pf_new = paste0("XDP_vs_", lvl)
tt = read.table(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, sep = "\t", check.names = F)
colnames(tt) = c(paste0("FDR.", mf, "XDP"), paste0("LFC.", mf, "XDP"), paste0("pval.", mf, "XDP"), "name", "class")
tt$LFC.GenotypeXDP = (-1) * tt$LFC.GenotypeXDP
tt$class = "XDP_vs_CON"
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf_new,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl = read.table(paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, sep = "\t", check.names = F)
bl$LFC = (-1) * bl$LFC
bl$class = "XDP_vs_CON"
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf_new,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


#(Untreated CON vs XDP_vanilla) and (Untreated XDP_noEdit vs XDP_vanilla)
ms = "Exp"
method = "glm"
mf = "Genotype1"
pf = paste0(lvl, "_vs_XDP")
date()
test1 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "nb")
date()

rnaPCA(test1@deSV, grp, c("Genotype1"), ntop = 500)

saveRDS(test1@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",pf,".rds"))
write.table(test1@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test1@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt00 = get_fdr(res = test1, term = "Genotype1", class = "test")
write.table(tt00, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl01 = get_ci(res = test1, term = "Genotype1XDP_noEdit", class = "XDP_noEdit")
bl02 = get_ci(res = test1, term = "Genotype1CON", class = "CON")
bl00 = as.data.frame(rbind(bl01,bl02))
write.table(bl00, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test1)
gc(verbose = F)


#(Untreated XDP_vanilla vs CON) and (Untreated XDP_noEdit vs CON)
lvl = "CON"
message(paste0("Untreated ", "XDP vs ", lvl))
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$Genotype %in% c("XDP", lvl)), ]
grp$Genotype = factor(grp$Genotype, levels = c(lvl, "XDP"))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp[grp$geno=="XDP", "Genotype1"] = "XDP_vanilla"
grp$Genotype1 = factor(grp$Genotype1, levels = c(lvl, "XDP_vanilla", "XDP_noEdit"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$ParentLine = factor(grp$ParentLine)
grp$MIN = factor(grp$MIN)

mat = mat0[,rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
condition_list = list(which(grp$Genotype1 == lvl), which(grp$Genotype1 == "XDP_vanilla"), which(grp$Genotype1 == "XDP_noEdit"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
mat = mat[which(expressedGenes), ]
message(nrow(mat)," genes expressed.")

ms = "Exp"
method = "glm"
mf = "Genotype1"
pf = paste0("XDP_vs_", lvl)
date()
test2 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "nb")
date()

rnaPCA(test2@deSV, grp, c("Genotype1"), ntop = 500)

saveRDS(test2@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",pf,".rds"))
write.table(test2@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test2@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test2@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt11 = get_fdr(res = test2, term = "Genotype1", class = "test")
write.table(tt11, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

tt11 = read.table(paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, sep = "\t", check.names = F)
tt11v = tt11[, c(1,3,5,7,8)]
tt11v$class = "XDPvanilla_vs_CON"
pf_new1 = "XDPvanilla_vs_CON"
write.table(tt11v, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf_new1,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
tt11n = tt11[, c(2,4,6,7,8)]
tt11n$class = "XDPnoEdit_vs_CON"
pf_new2 = "XDPnoEdit_vs_CON"
write.table(tt11n, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf_new2,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl001 = get_ci(res = test2, term = "Genotype1XDP_noEdit", class = "XDP_noEdit_vs_CON")
bl002 = get_ci(res = test2, term = "Genotype1XDP_vanilla", class = "XDP_vanilla_vs_CON")
bl11 = as.data.frame(rbind(bl001,bl002))
write.table(bl11, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test2)
gc(verbose = F)





#Untreated dSVA vs CON
lvl = "dSVA"
message(paste0("Untreated ", lvl, " vs CON"))
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$Genotype %in% c("CON", lvl)), ]
grp$Genotype = factor(grp$Genotype, levels = c("CON", lvl))
grp$Genotype1 = grp$geno
#grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
#grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
#grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", lvl))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$ParentLine = factor(grp$ParentLine)
grp$MIN = factor(grp$MIN)

mat = mat0[,rownames(grp)]
CPM = t(t(mat) / grp$sizeFactor)
condition_list = list(which(grp$Genotype == lvl), which(grp$Genotype == "CON"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
mat = mat[which(expressedGenes), ]
message(nrow(mat)," genes expressed.")

ms = "Exp"
method = "glm"
mf = "Genotype"
pf = paste0(lvl, "_vs_CON")
date()
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchGrowth"), remove_terms = "BatchGrowth", random_terms = NULL, distro_family = "nb")
date()

rnaPCA(test@deSV, grp, "Genotype", ntop = 500)

saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",pf,".rds"))
write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt = get_fdr(res = test, term = "Genotype", class = pf)
write.table(tt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl = get_ci(res = test, term = paste0("Genotype", lvl), class = pf)
write.table(bl, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test)
gc(verbose = F)

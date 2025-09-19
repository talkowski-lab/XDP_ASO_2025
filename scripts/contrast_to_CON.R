source("~/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")


#############################################
#RNASeq
#############################################
#Untreated dSVA vs CON
#############################################
ms = "Exp"
method = "crossModel"
mf = "Genotype"
lvl1 = "dSVA"
lvl2 = "CON"
pf = paste0(lvl1, "_vs_", lvl2)
trt = cross_model_z_test(term1 = "Genotype",
                         term2 = NULL,
                         lvl1 = lvl1,
                         lvl2 = lvl2,
                         cond1 = "model_params_colData_glmm_Exp_Genotype_dSVA_vs_XDP.txt",
                         cond2 = "model_params_colData_glm_Exp_Genotype_CON_vs_XDP.txt",
                         beta1 = "model_params_beta_glmm_Exp_Genotype_dSVA_vs_XDP.txt",
                         beta2 = "model_params_beta_glm_Exp_Genotype_CON_vs_XDP.txt",
                         stderr1 = "model_params_stderr_glmm_Exp_Genotype_dSVA_vs_XDP.txt",
                         stderr2 = "model_params_stderr_glm_Exp_Genotype_CON_vs_XDP.txt")

write.table(trt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


ms = "IR"
method = "glmGaussianZsum"
mf = "Genotype"
lvl1 = "dSVA"
lvl2 = "CON"
pf = paste0(lvl1, "_vs_", lvl2)
trt = cross_model_z_test(term1 = "Genotype",
                         term2 = NULL,
                         lvl1 = lvl1,
                         lvl2 = lvl2,
                         n1 = 30,
                         n2 = 8,
                         beta1 = "model_params_beta_glm_IR_Genotype.txt",
                         beta2 = "model_params_beta_glm_IR_Genotype.txt",
                         stderr1 = "model_params_stderr_glm_IR_Genotype.txt",
                         stderr2 = "model_params_stderr_glm_IR_Genotype.txt")

write.table(trt, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


#############################################
#ASO880-treated vs CON
#############################################
ms = "Exp"
method = "crossModel"
mf1 = "Treatment"
mf2 = "Genotype"
lvl1 = "880"
lvl2 = "CON"
pf = paste0(mf1, lvl1, "_vs_", mf2, lvl2)
trt = cross_model_z_test(term1 = mf1,
                         term2 = mf2,
                         lvl1 = lvl1,
                         lvl2 = lvl2,
                         cond1 = "model_params_colData_glmm_Exp_Treatment_880.txt",
                         cond2 = "model_params_colData_glm_Exp_Genotype_CON_vs_XDP.txt",
                         beta1 = "model_params_beta_glmm_Exp_Treatment_880.txt",
                         beta2 = "model_params_beta_glm_Exp_Genotype_CON_vs_XDP.txt",
                         stderr1 = "model_params_stderr_glmm_Exp_Treatment_880.txt",
                         stderr2 = "model_params_stderr_glm_Exp_Genotype_CON_vs_XDP.txt")

write.table(trt, paste0("res_FDR_LFC_",method,"_",ms,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


ms = "IR"
method = "crossModel"
mf1 = "Treatment"
mf2 = "Genotype"
lvl1 = "880"
lvl2 = "CON"
pf = paste0(mf1, lvl1, "_vs_", mf2, lvl2)
trt = cross_model_z_test(term1 = mf1,
                         term2 = mf2,
                         lvl1 = lvl1,
                         lvl2 = lvl2,
                         n1 = 32,
                         n2 = 8,
                         beta1 = "model_params_beta_glm_IR_Treatment_880.txt",
                         beta2 = "model_params_beta_glm_IR_Genotype.txt",
                         stderr1 = "model_params_stderr_glm_IR_Treatment_880.txt",
                         stderr2 = "model_params_stderr_glm_IR_Genotype.txt")

write.table(trt, paste0("res_FDR_LFC_",method,"_",ms,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


#############################################
#ASO131-treated vs CON
#############################################
ms = "Exp"
method = "crossModel"
mf1 = "Treatment"
mf2 = "Genotype"
lvl1 = "131"
lvl2 = "CON"
pf = paste0(mf1, lvl1, "_vs_", mf2, lvl2)
trt = cross_model_z_test(term1 = mf1,
                         term2 = mf2,
                         lvl1 = lvl1,
                         lvl2 = lvl2,
                         cond1 = "model_params_colData_glmm_Exp_Treatment_131.txt",
                         cond2 = "model_params_colData_glm_Exp_Genotype_CON_vs_XDP.txt",
                         beta1 = "model_params_beta_glmm_Exp_Treatment_131.txt",
                         beta2 = "model_params_beta_glm_Exp_Genotype_CON_vs_XDP.txt",
                         stderr1 = "model_params_stderr_glmm_Exp_Treatment_131.txt",
                         stderr2 = "model_params_stderr_glm_Exp_Genotype_CON_vs_XDP.txt")

write.table(trt, paste0("res_FDR_LFC_",method,"_",ms,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


ms = "IR"
method = "crossModel"
mf1 = "Treatment"
mf2 = "Genotype"
lvl1 = "131"
lvl2 = "CON"
pf = paste0(mf1, lvl1, "_vs_", mf2, lvl2)
trt = cross_model_z_test(term1 = mf1,
                         term2 = mf2,
                         lvl1 = lvl1,
                         lvl2 = lvl2,
                         n1 = 33,
                         n2 = 8,
                         beta1 = "model_params_beta_glm_IR_Treatment_131.txt",
                         beta2 = "model_params_beta_glm_IR_Genotype.txt",
                         stderr1 = "model_params_stderr_glm_IR_Treatment_131.txt",
                         stderr2 = "model_params_stderr_glm_IR_Genotype.txt")

write.table(trt, paste0("res_FDR_LFC_",method,"_",ms,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")



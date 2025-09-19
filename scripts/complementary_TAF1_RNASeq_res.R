source("/home/dadi/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")


targets = c("dSVA", "131", "307", "876", "877", "879", "880", "881")
for (x in targets){
  if (x == "dSVA"){
    mf = "Genotype"
    trt = "dSVA_vs_XDP"
  } else {
    mf = "Treatment"
    trt = x
  }
  message(mf, " ", x, " vs XDP")
  
  beta_file = paste0("model_params_beta_glmm_Exp_", mf, "_", trt, ".txt")
  pval_file = paste0("model_params_pval_glmm_Exp_", mf, "_", trt, ".txt")
  fdr_file  = paste0("model_params_fdr_glmm_Exp_", mf, "_", trt, ".txt")
  stderr_file  = paste0("model_params_stderr_glmm_Exp_", mf, "_", trt, ".txt")
  
  stats = get_fdr(beta = beta_file, pval = pval_file, fdr = fdr_file, term = mf, class = paste0(mf, x, "_vs_XDP"))
  if (x != "dSVA"){colnames(stats) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")}
  ci = get_ci(beta = beta_file, stderr = stderr_file, pval = pval_file, fdr = fdr_file, term = paste0(mf, x), class = paste0(mf, x, "_vs_XDP"))
  
  write.table(stats, paste0("res_FDR_LFC_glmm_Exp_", mf ,"_", trt, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(ci, paste0("res_LFC_CI_glmm_Exp_", mf, "_", trt, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  rm(mf, trt, beta_file, pval_file, fdr_file, stderr_file, stats, ci)
  gc(verbose = F)
}

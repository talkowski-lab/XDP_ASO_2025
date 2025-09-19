library(ggplot2)


setwd("~/projects/xdp_aso/nsc/output/")


merge_vcov_from_model_list = function(model_list){
  tt = lapply(model_list, vcov)
  tt1 = lapply(names(tt), FUN = function(x){tmp = as.data.frame(as.matrix(tt[[x]])); rownames(tmp) = NULL; tmp[["name"]] = x; return(tmp)}) 
  tt2 = do.call(rbind, tt1)
  return(tt2)
}


lfs = list.files("RData")
model_files = lfs[grep("models_", lfs)]
for (x in model_files){
  message(x)
  models0 = readRDS(paste0("RData/", x))
  txt = gsub(".rds", ".txt", x)
  vcov_file = gsub("models_", "model_params_vcov_", txt)
  pval_file = gsub("models_", "model_params_pval_", txt)
  resid_file = gsub("models_", "model_params_resid_", txt)
  desv_file = gsub("models_", "model_params_deSV_", txt)
  
  deSV = read.table(desv_file, header = T, sep = "\t", check.names = F)
  
  model_failed = sapply(models0, is.null)
  models = models0[which(!model_failed)]
  vcov_all = merge_vcov_from_model_list(models)
  p_all = t(sapply(models, function(x){coef(summary(x))[,4]}))
  resid_all = t(sapply(models, resid, type = "response"))
  rownames(p_all) = names(models)
  rownames(resid_all) = names(models)
  colnames(resid_all) = colnames(deSV)
  
  write.table(vcov_all, vcov_file, row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(as.matrix(p_all), pval_file, row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(as.matrix(resid_all), resid_file, row.names = T, col.names = T, quote = F, sep = "\t")
  
  rm(models0, models, vcov_all, p_all, resid_all, deSV)
  gc(verbose = F)
}
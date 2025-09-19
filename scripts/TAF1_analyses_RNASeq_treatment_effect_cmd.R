library(optparse)
library(lme4)
library(pbkrtest)
library(Hmisc)
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
option_list = list(make_option(c("-a", "--aso"), type = "character", default = NULL, help = "ASO #", metavar = "character"),
                   make_option(c("-r", "--reference"), type = "character", default = NULL, help = "One of 'XDP' for untreated XDP, 'CON' for untreated control, or 'CON_ASO' for treated control as the reference level", metavar = "character"))
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

aso = opt$aso
ref = opt$reference
if (is.null(aso)){stop("An ASO number (i.e. three digits) must be provided via '-a' option.")}
if ((aso %in% c("880", "131", "307", "876", "877", "879", "881")) == FALSE){stop("The ASO number must be one of 880, 131, 307, 876, 877, 879 and 881.")}
if (is.null(ref)){stop("Reference level must be provided via '-a' option. Either 'XDP' for untreated XDP or 'CON' for untreated control as the reference level.")}
if ((ref %in% c("XDP", "CON", "CON_ASO")) == FALSE){stop("The reference level must be one of 'XDP', 'CON' or 'CON_ASO'.")}

key_name = "ENSG00000147133"
mat0 = read.table("ASO_RNAseq_R1.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)


#ASO treated vs untreated XDP
if (ref == "XDP"){
  trt = c()
  stats = c()
  plan(multicore, workers = 18)
  message(paste0("ASO", aso, "-treated vs untreated XDP"))
  
  grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  grp = grp[which(grp$Treatment %in% c("NO", aso)), ]
  grp = grp[which(grp$Genotype %in% c("XDP")), ]
  grp$Treatment = factor(grp$Treatment,levels=c("NO", aso))
  grp$BatchGrowth = factor(grp$BatchGrowth)
  grp$BatchTreat = factor(grp$BatchTreat)
  grp$MIN = factor(grp$MIN)
  grp$ParentLine = factor(grp$ParentLine)
  
  mat = mat0[,rownames(grp)]
  CPM = t(t(mat) / grp$sizeFactor)
  condition_list = list(which(grp$Treatment == "NO"), which(grp$Treatment == aso))
  expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
  mat = mat[which(expressedGenes), ]
  message(nrow(mat)," genes expressed.")
  
  ms = "Exp"
  method = "glmm"
  mf = "Treatment"
  date()
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchTreat", "MIN"), remove_terms = "BatchTreat", random_terms = "MIN", distro_family = "nb")
  date()
  
  saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,"_",aso,".rds"))
  write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  tmp = get_fdr(res = test, term = "Treatment", class = aso, keyword = NULL)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(res = test, term = paste0("Treatment",aso), class = aso))
  
  rm(test, tmp)
  gc(verbose = F)
  
  plan(sequential)
  
  write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
} else {
  trt = c()
  stats = c()
  if (ref == "CON"){
    message(paste0("ASO", aso, "-treated XDP vs untreated CON"))
  } else {
    message(paste0("ASO", aso, "-treated XDP vs ASO", aso, "-treated CON"))
  }
  
  grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  grp = grp[which(grp$Treatment %in% c("NO", aso)), ]
  grp = grp[which(grp$Genotype %in% c("XDP", "CON")), ]
  grp$Treatment = factor(grp$Treatment,levels=c("NO", aso))
  grp$BatchGrowth = factor(grp$BatchGrowth)
  grp$BatchTreat = factor(grp$BatchTreat)
  grp$MIN = factor(grp$MIN)
  grp$ParentLine = factor(grp$ParentLine)
  grp$Treated = paste0(grp$Genotype, "_", grp$Treatment)
  if (ref == "CON"){
    grp = grp[which((grp$Treated == "CON_NO") | (grp$Treated == paste0("XDP_", aso))), ]
    grp[grp$Treated == "CON_NO", "Treated"] = "CON"
    grp[grp$Treated == paste0("XDP_", aso), "Treated"] = paste0("XDP", aso)
    grp$Treated = factor(grp$Treated, levels = c("CON", paste0("XDP", aso)))
  } else {
    grp = grp[which((grp$Treated == paste0("CON_", aso)) | (grp$Treated == paste0("XDP_", aso))), ]
    grp[grp$Treated == paste0("CON_", aso), "Treated"] = paste0("CON", aso)
    grp[grp$Treated == paste0("XDP_", aso), "Treated"] = paste0("XDP", aso)
    grp$Treated = factor(grp$Treated, levels = c(paste0("CON", aso), paste0("XDP", aso)))
  }
  
  mat = mat0[,rownames(grp)]
  CPM = t(t(mat) / grp$sizeFactor)
  if (ref == "CON"){
    condition_list = list(which(grp$Treated == "CON"), which(grp$Treated == paste0("XDP", aso)))
  } else {
    condition_list = list(which(grp$Treated == paste0("CON", aso)), which(grp$Treated == paste0("XDP", aso)))
  }
  expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
  mat = mat[which(expressedGenes), ]
  message(nrow(mat)," genes expressed.")
  
  ms = "Exp"
  method = "glm"
  mf = "Treated"
  if (ref == "CON"){
    pf = paste0("XDP", aso, "_vs_CON")
  } else {
    pf = paste0("XDP", aso, "_vs_CON", aso)
  }
  date()
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf, "BatchTreat"), remove_terms = "BatchTreat", random_terms = NULL, distro_family = "nb")
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
  
  tmp = get_fdr(res = test, term = "Treated", class = pf, keyword = NULL)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "pval.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(res = test, term = paste0("TreatedXDP", aso), class = pf))
  
  rm(test, tmp)
  gc(verbose = F)
  
  write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}

message("All Done.")

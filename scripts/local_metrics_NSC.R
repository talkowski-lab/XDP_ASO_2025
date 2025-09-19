library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
library(ggbeeswarm)
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/projects/xdp_aso/nsc/scripts/glmm.R")


setwd("~/projects/xdp_aso/nsc/output/")


#####################################################
#Building AS models for treated XDP and untreated CON
#####################################################
#RNACapSeq
key_name = "X_70644089_70659646_1"

grp = read.table("ASO_Capseq_R1_2.metadata.final.4modelingTrt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = read.table("ASO_Capseq_R1_2.i32splicing.final.4modelingTrt.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat[, rownames(grp)]

ms = "AS"
method1 = "glm"
mf = "Trt"

for (v in c("dSVA", "XDP126", "XDP131", "XDP134", "XDP307", "XDP849", "XDP852", "XDP876", "XDP877", "XDP879", "XDP880", "XDP881", "XDP883", "XDPMAL")){
  target = v
  mf1 = paste0(mf, "_OnebyOne_", target)
  message(mf1)
  
  ref1 = "CON"
  ref2 = "XDP"
  grp1 = grp[which(grp[[mf]] == ref1 | grp[[mf]] == ref2 |grp[[mf]] == target), ]
  mat1 = mat[, rownames(grp1)]
  
  grp1[[mf]] = factor(grp1[[mf]], levels = c(ref1, ref2, target))
  
  test = buildDEmodel(feature_matrix = mat1, sampleTable = grp1, method = method1, estimate_offset = F, offset_str = "sizeFactor", apply_sva = F, test_terms = c(mf), remove_terms = NULL, random_terms = "MIN", distro_family = "nb")
  
  asi32 = t(exp(test@deSV))
  asi32 = cbind(asi32, grp1)
  write.table(asi32, paste0("Supplementary_Table_SX_NSC_RNACapSeq_AS-i32_corrected_counts_",method1,"_",ms,"_",mf1,"_vs_CON_vs_XDP.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  

  if (target  == "dSVA"){
    asi32$MIN = factor(asi32$MIN)
    ggplot(asi32, aes(x = Trt, y = X_70644089_70659646_1, fill = Trt)) + geom_boxplot() + geom_jitter() + geom_fill(c("#377EB8", "#E41A1C", "#7570B3")) + ylab("i32-AS Corrected Counts") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/boxplot_NSC_i32AS_corrected_counts_CON_vs_XDP_vs_dSVA_OnebyOne_by_clone.pdf"), width = 3, height = 9)
    
    asi32_ind = aggregate(X_70644089_70659646_1~MIN+Trt, asi32, FUN = median)
    ggplot(asi32_ind, aes(x = Trt, y = X_70644089_70659646_1, fill = Trt)) + geom_boxplot() + geom_jitter() + geom_fill(c("#377EB8", "#E41A1C", "#7570B3")) + ylim(c(0, 600)) + ylab("i32-AS Corrected Counts") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/boxplot_NSC_i32AS_corrected_counts_CON_vs_XDP_vs_dSVA_OnebyOne_by_MIN.pdf"), width = 3, height = 9)
    ggplot(asi32_ind, aes(x = Trt, y = X_70644089_70659646_1, fill = Trt)) + geom_violin() + geom_jitter() + geom_fill(c("#377EB8", "#E41A1C", "#7570B3")) + ylim(c(0, 600)) + ylab("i32-AS Corrected Counts") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/violinplot_NSC_i32AS_corrected_counts_CON_vs_XDP_vs_dSVA_OnebyOne_by_MIN.pdf"), width = 3, height = 9)
    ggplot(asi32_ind, aes(x = Trt, y = log10(X_70644089_70659646_1), fill = Trt)) + geom_violin() + geom_jitter() + geom_fill(c("#377EB8", "#7570B3", "#E41A1C"))+ ylim(c(0, 3)) + ylab("i32-AS Corrected Counts\n(log10-transformed)") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/violinplot_log10_NSC_i32AS_corrected_counts_CON_vs_XDP_vs_dSVA_OnebyOne_by_MIN.pdf"), width = 3, height = 9)
    
    ggplot(asi32, aes(x = Trt, y = log10(X_70644089_70659646_1), fill = Trt)) + geom_boxplot() + geom_jitter() + geom_fill(c("#377EB8", "#E41A1C", "#7570B3")) + ylab("i32-AS Corrected Counts\n(log10-transformed)") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/boxplot_log10_NSC_i32AS_corrected_counts_CON_vs_XDP_vs_dSVA_OnebyOne_by_clone.pdf"))
    
    ggplot(asi32, aes(x = Trt, y = log10(X_70644089_70659646_1), fill = Trt)) + geom_violin(width=1.2) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + geom_jitter() + geom_fill(c("#377EB8", "#E41A1C", "#7570B3")) + ylab("i32-AS Corrected Counts\n(log10-transformed)") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/violinboxplot_log10_NSC_i32AS_corrected_counts_CON_vs_XDP_vs_dSVA_OnebyOne_by_clone.pdf"))
    
    
    1 - mean(asi32[asi32$Trt == "dSVA", key_name])/mean(asi32[asi32$Trt == "XDP", key_name])
    
    asi32_sub = asi32[which(asi32$geno != "XDP_dSVA"), ]
    asi32_sub$geno = factor(as.vector(asi32_sub$geno))
    ais32_sub$MIN = factor(as.vector(asi32_sub$MIN))
    
    ggplot(asi32_sub, aes(x = geno, y = X_70644089_70659646_1)) + geom_boxplot() + geom_jitter(aes(color = MIN)) + ylab("i32-AS Corrected Counts") + xlab("") + geom_noBG()
    ggsave(paste0("pdf/boxplot_NSC_i32AS_corrected_counts_CON_vs_XDPNaive_vsXDPUnedit_by_clone.pdf"))
    t.test(asi32_sub[asi32_sub$geno == "CON", key_name], asi32_sub[asi32_sub$geno == "XDP", key_name])$p.value            # 4.634723e-07
    t.test(asi32_sub[asi32_sub$geno == "CON", key_name], asi32_sub[asi32_sub$geno == "XDP_non_edit", key_name])$p.value   # 2.0726e-07
    t.test(asi32_sub[asi32_sub$geno == "XDP", key_name], asi32_sub[asi32_sub$geno == "XDP_non_edit", key_name])$p.value   # 0.4474019
    
    
    ggplot(asi32_sub, aes(x = X_70644089_70659646_1, fill = Trt)) + geom_density() + xlab("i32-AS Corrected Counts (log10-transformed)") + ylab("Density") + geom_noBG()
    
  }
  
  write.table(test@vcov, paste0("model_params_vcov_",method1,"_",ms,"_",mf1,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@pval, paste0("model_params_pval_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@resid, paste0("model_params_resid_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method1,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  mc = paste0("OnebyOne_", v, "_vs_CON")
  aa = get_fdr(test, term = paste0("Trt", v), class = mc, keyword = key_name)
  bb = get_ci(test, term = paste0("Trt", v), class = mc)
  write.table(aa, file = paste0("res_FDR_LFC_",method1,"_",ms,"_",mf,"_",mc,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(bb, file = paste0("res_LFC_CI_",method1,"_",ms,"_",mf,"_",mc,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}


#RNASeq Exp
key_name = "ENSG00000147133"
sel = c("CONNO", "XDPNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")

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

for (v in c("dSVA", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")){
  target = v
  mf1 = paste0(mf, "_OnebyOne_", target)
  message(mf1)
  
  ref1 = "CON"
  ref2 = "XDP"
  grp1 = grp[which(grp[[mf]] == ref1 | grp[[mf]] == ref2 |grp[[mf]] == target), ]
  mat1 = mat[, rownames(grp1)]
  
  grp1[[mf]] = factor(grp1[[mf]], levels = c(ref1, ref2, target))
  
  test = buildDEmodel(feature_matrix = mat1, sampleTable = grp1, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")
  
  #saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,".rds"))
  write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf1,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf1,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  mc = paste0("OnebyOne_", v, "_vs_CON")
  aa = get_fdr(test, term = paste0("Trt", v), class = mc)
  bb = get_ci(test, term = paste0("Trt", v), class = mc)
  write.table(aa, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",mc,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(bb, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",mc,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}


#RNASeq IR
version_name = "avg"
if (version_name == "avg"){version_num = "_avg"} else {version_num = ""}

key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
sel = c("CONNO", "XDPNO", "dSVANO", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")

mat0 = read.table(paste0("ASO_RNAseq_", version_name, ".IR.txt"), header = T, row.names = 1, sep = "\t", check.names = F)

#remove outlier untreated XDP samples
bad_idx = c(grep("NO_ASO_XDP_33363C_n_15", colnames(mat0)), grep("NO_ASO_XDP_non_edit_35833A_ISO_SVA_1B5", colnames(mat0)))
mat0 = mat0[, -bad_idx]

meta0 = read.table("ASO_RNAseq_new.569_namecorrected.metadata.tab", header = T, check.names = F, sep = "\t", row.names = 1)

grp = meta0[colnames(mat0), ]
table(grp$Genotype)
grp$Trt = paste0(grp$Genotype, grp$Treatment)
grp = grp[which(grp$Trt %in% sel),]
grp$Trt = factor(grp$Trt, levels = sel)
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)
grp$SequencingBatch = factor(grp$SequencingBatch)

mat = mat0[,rownames(grp)]
mat = 1000 * mat / get_intron_length(rownames(mat))

grp$Trt = gsub("NO$", "", as.vector(grp$Trt))
sel0 = gsub("NO$", "", sel)
grp$Trt = factor(grp$Trt, levels = sel0)


ms = "IR"
method = "glm"
mf = "Trt"
mem_size_in_mb = 1024
options(future.globals.maxSize= mem_size_in_mb * (1024^2))

for (v in c("dSVA", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")){
  target = v
  mf1 = paste0(mf, "_OnebyOne_", target)
  message(mf1)
  
  ref1 = "CON"
  ref2 = "XDP"
  grp1 = grp[which(grp[[mf]] == ref1 | grp[[mf]] == ref2 |grp[[mf]] == target), ]
  mat1 = mat[, rownames(grp1)]
  
  grp1[[mf]] = factor(grp1[[mf]], levels = c(ref1, ref2, target))
  grp1[["BatchTreat"]] = factor(grp1[["BatchTreat"]])
  
  test = buildDEmodel(feature_matrix = mat1, sampleTable = grp1, method = method, estimate_offset = F, offset_str = NULL, apply_sva = F, test_terms = c(mf, "BatchTreat"), remove_terms = c("BatchTreat"), random_terms = NULL, distro_family = "gaussian")
  
  #saveRDS(test@models, paste0("RData/models_",method,"_",ms,"_",mf,".rds"))
  write.table(test@vcov, paste0("model_params_vcov_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table(test@colData, paste0("model_params_colData_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@pval, paste0("model_params_pval_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@resid, paste0("model_params_resid_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@status, paste0("model_params_warnings_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

  mc = paste0("OnebyOne_", v, "_vs_CON")
  stat = get_fdr(test, term = paste0("Trt", v), class = mc)
  write.table(stat, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",mc,"_new", version_num, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  ci = get_ci(test, term = paste0("Trt", v), class = mc)
  write.table(ci, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",mc,"_new", version_num, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
  
  rm(test, mat1, grp1, stat, ci)
  gc()
}


for (v in c("dSVA", "XDP131", "XDP880", "XDP881", "XDP877", "XDP879", "XDP876", "XDP307")){
  target = v
  mf1 = paste0(mf, "_OnebyOne_", target)
  message(mf1)
  
  ref1 = "CON"
  ref2 = "XDP"
  grp1 = grp[which(grp[[mf]] == ref1 | grp[[mf]] == ref2 |grp[[mf]] == target), ]
  mat1 = mat[, rownames(grp1)]
  
  grp1[[mf]] = factor(grp1[[mf]], levels = c(ref1, ref2, target))
  grp1[["BatchTreat"]] = factor(grp1[["BatchTreat"]])
  
  
  tmp = read.table(paste0("model_params_deSV_",method,"_",ms,"_",mf1,"_new", version_num, ".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
  tmp = tmp[key_name, ]
  tmp = as.matrix(tmp)
  iri32 = t(tmp)
  iri32 = cbind(iri32, grp1)
  write.table(iri32, paste0("Supplementary_Table_SX_NSC_RNASeq_IR-i32_corrected_values_",method1,"_",ms,"_",mf1,"_vs_CON_vs_XDP.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
}



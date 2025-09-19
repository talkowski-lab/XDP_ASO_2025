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
#Cell paper study
makenames=function(x){
  if (x[2]=="Fibroblasts"){
    x[2]="Fibro"
  }
  else if (x[2]=="Neurons"){
    x[2]="Neuron"
  }
  paste(x[2],x[3],"rnaseq",x[4],x[1],"deseq2.count",sep="_")
}

key_name = "ENSG00000147133"
mat0_p1 = read.table("cell_paper_all_counts_raw.txt", header = T, row.names = 1, sep = "\t", check.names = F)
tissue = "NSC"

#Cell Paper XDP vs CON
mytab=read.table("cellpaper_metadata.txt", header=T, sep="\t")
sampleNames=apply(mytab, 1, makenames)
rownames(mytab)=sampleNames
mysz= read.table("cellpaper_sizeFactors.txt")
#tmp=strsplit(as.vector(mysz$V1),"_deseq2.count")
#libNames=sapply(tmp,function(x) x[1])
rownames(mysz)=as.vector(mysz$V1)
mytab$fileName=sampleNames

#make the full design matrix
group=mytab[,c("Sample", "fileName", "CellType", "Genotype", "Batch", "MINnumber", "Sex")]
row.names(group) = NULL
group$Sample = sampleNames
group$MINnumber = paste0("ID_", group$MINnumber)
group$CellType = as.vector(group$CellType)
group$CellType[group$CellType == "Fibroblasts"] = "Fibro"
group$CellType[group$CellType == "Neurons"] = "Neuron"
group$lib.size = mysz[sampleNames, 2]

group$Sample = factor(group$Sample)
group$CellType = factor(group$CellType)
group$Genotype = factor(group$Genotype, levels = c("Control", "Carrier", "XDP"))
group$Batch = factor(group$Batch)
group$MINnumber = factor(group$MINnumber)
group$Sex = factor(group$Sex, levels=c("M", "F"))

grp = group[group$CellType == tissue,]
grp = grp[which(grp$Genotype %in% c("Control", "XDP")), ]
grp$Sample = factor(grp$Sample)
grp$CellType = factor(grp$CellType)
grp$Genotype = factor(grp$Genotype, levels=c("Control", "XDP"))

if (tissue == "Fibro"){
  grp[grp$Batch == "mergedFibro" | grp$Batch == "mergedFibro2", "Batch"] = "Oct15"
} else if (tissue == "Neuron"){
  grp[grp$Batch == "NeuronRapid", "Batch"] = "Aug16"
}
grp$Batch=factor(grp$Batch)
grp$MINnumber = factor(grp$MINnumber)
grp$Sex=factor(grp$Sex, levels=c("M", "F"))

grp$MIN = gsub("ID_", "", grp$MINnumber)

grp1 = grp
rownames(grp1) = as.vector(grp1$fileName)
mat_p1 = mat0_p1[, rownames(grp1)]
CPM_p1 = t(t(mat_p1) / grp1$lib.size)
#condition_list = list(which(grp$Genotype == "Control"), which(grp$Genotype == "XDP"))
#expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 10, pct = 0.5)
#message(length(which(expressedGenes)), " genes expressed.")
#mat = mat[which(expressedGenes), ]












#######################################################################
#Current Study: Untreated dSVA vs XDP
mat0_p2 = read.table("ASO_RNAseq_R1.exp.txt", header = T, row.names = 1, sep = "\t", check.names = F)

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

grp2 = grp
mat_p2 = mat0_p2[, rownames(grp2)]
CPM_p2 = t(t(mat_p2) / grp2$sizeFactor)

mat_p1 = mat_p1[rownames(mat_p2), ]
CPM_p1 = CPM_p1[rownames(mat_p2), ]
mat = cbind(mat_p1, mat_p2)
CPM = cbind(CPM_p1, CPM_p2)

grp1$sizeFactor = grp1$lib.size
grp1$BatchStudy = "CellPaper"
grp1$Genotype = as.vector(grp1$Genotype)
grp1[which(grp1$Genotype == "Control"), "Genotype"] = "CON"
grp2$BatchStudy = "ASO_dSVA"
grp2$Genotype = as.vector(grp2$Genotype)

grp = rbind(grp1[, c("Genotype", "BatchStudy", "MIN", "sizeFactor")], grp2[, c("Genotype", "BatchStudy", "MIN", "sizeFactor")])
grp$Genotype = factor(grp$Genotype, levels = c("CON", "XDP"))
grp$BatchStudy = factor(grp$BatchStudy)
grp$MIN = factor(grp$MIN)
pedigree = c("22", "12", "12", "12", "17", "17", "20", "26", "32", "35", "22")
names(pedigree) = names(table(grp$MIN))
grp$Pedigree = unname(pedigree[as.vector(grp$MIN)])
grp$Pedigree = factor(grp$Pedigree)
  
condition_list = list(which(grp$Genotype == "CON"), which(grp$Genotype == "XDP"))
expressedGenes = apply(CPM, 1, countsFilter, conditionList = condition_list, cut = 5, pct = 0.5)
message(length(which(expressedGenes))," genes expressed.")
mat = mat[which(expressedGenes), ]


#######################################################################
trt = c()
stats = c()

ms = "Exp"
method = "glmm"
mf = "Genotype"
pf = "XDP_vs_CON_mergedData"
message(paste0("Start building models..."))
date()
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = NULL, distro_family = "nb")
#plan(multicore, workers = 18)
#test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, estimate_offset = T, offset_str = NULL, apply_sva = T, test_terms = c(mf), remove_terms = NULL, random_terms = "MIN", distro_family = "nb")
#plan(sequential)
date()

#test@beta[key_name, ]
#norm = test@deSV
#rnaPCA(norm, grp, c("Genotype"), ntop = 100)
#rnaPCA(norm, grp, c("MIN"), ntop = 500)
#new = cbind(Exp = unlist(norm[key_name, ]), grp)
#ggplot(new, aes(Genotype, Exp)) + geom_boxplot() + geom_noBG() + labs(x = "", y = "Expression (log-scale)", color = "Donor")


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

tmp = get_fdr(res = test, term = "Genotype", class = "XDP_vs_CON", keyword = NULL)
colnames(tmp) = c("FDR.GenotypeXDP", "LFC.GenotypeXDP", "pval.GenotypeXDP", "name", "class")
stats = rbind(stats, tmp)
trt = rbind(trt, get_ci(res = test, term = "GenotypeXDP", class = "XDP_vs_CON"))

rm(test, tmp)
gc(verbose = F)

write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,"_",pf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


message("All Done.")


ms = "Exp"
method = "glmm"
mf = "Genotype"
pf = "XDP_vs_CON_mergedData"
grp = read.table(paste0("model_params_colData_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, row.names = 1, sep = "\t", check.names =  F)
pedigree = c("22", "12", "12", "12", "17", "17", "20", "26", "32", "35", "22")
names(pedigree) = names(table(grp$MIN))
grp$Pedigree = unname(pedigree[as.vector(factor(grp$MIN))])
grp$Pedigree = factor(grp$Pedigree)

norm = read.table(paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",pf,".txt"), header = T, row.names = 1, sep = "\t", check.names =  F)
rnaPCA(norm, grp, "BatchStudy", ntop = 500, pointSize = 1.5, returnData = F)
rnaPCA(norm, grp, "Genotype", ntop = 500, pointSize = 1.5, returnData = F)
rnaPCA(norm, grp, "Pedigree", ntop = 500, pointSize = 1.5, returnData = F)
rnaPCA(norm, grp, c("Genotype","BatchStudy"), ntop = 500, pointSize = 1.5, returnData = F)

norm0 = log(t(t(mat) / grp$buildDEmodel.sizeFactor))
rnaPCA(norm0, grp, "BatchStudy", ntop = 500, pointSize = 1.5)
rnaPCA(norm0, grp, "Genotype", ntop = 500, pointSize = 1.5)
rnaPCA(norm0, grp, "Pedigree", ntop = 500, pointSize = 1.5)

#new = cbind(Exp = unlist(norm[key_name, ]), grp)
#ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Expression (log-scale)", color = "Donor")

norm0 = norm
norm = norm[,grep(".count", colnames(norm0))]
rs = rowSums(norm)
norm_clean = norm[which(!is.infinite(rs) & !is.na(rs)), ]
sampleDists <- dist(t(norm_clean))
hc = hclust(sampleDists, method = "ward.D2")

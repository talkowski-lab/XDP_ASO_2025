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


#############################################
#RNASeq
#############################################
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

mat0 = read.table("ASO_RNAseq_R1.IR.txt", header = T, row.names = 1, sep = "\t", check.names = F)
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
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, apply_sva = F, test_terms = c(mf,"BatchGrowth"), remove_terms = c("BatchGrowth"), random_model_str = "(1|ParentLine)", offset_str = NULL, distro_family = "gaussian")
write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

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
test1 = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = method, apply_sva = F, test_terms = c(mf,"BatchGrowth"), remove_terms = c("BatchGrowth"), random_model_str = "(1|ParentLine)", offset_str = NULL, distro_family = "gaussian")
write.table(test1@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(test1@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

tt00 = get_fdr(test1,"Genotype1","test")
write.table(tt00, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

bl01 = get_ci(test1,"Genotype1XDP_noEdit", "XDP_noEdit")
bl02 = get_ci(test1,"Genotype1CON", "CON")
bl03 = get_ci(test1,"Genotype1dSVA", "dSVA")
bl00 = as.data.frame(rbind(bl01,bl02,bl03))
write.table(bl00, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

rm(test, test1)
gc(verbose = F)


#Treatment
mat0 = read.table("ASO_RNAseq_R1.IR.txt", header = T, row.names = 1, sep = "\t", check.names = F)

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
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_raw.jpg"))
  #new = cbind(Exp = unlist(norm[key_name, ]), grp)
  #ggplot(new, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_IR_ASO_",aso,"_effect_raw.jpg"))
  
  ms = "IR"
  method = "glm"
  #method = "glmm"
  mf = "Treatment"
  #plan("multisession", workers = 4)
  #date()
  test = buildDEmodel(mat, grp, method, F, c(mf,"BatchTreat"), c("BatchTreat"), "(1|ParentLine)", offset_str = NULL, distro_family = "gaussian")
  #test = buildDEmodel(mat, grp, method, F, c(mf), NULL, "(1|ParentLine) + (1|BatchTreat)", offset_str = NULL, distro_family = "gaussian")
  #date()
  #plan("sequential")
  write.table(test@beta, paste0("model_params_beta_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@stderr, paste0("model_params_stderr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@fdr, paste0("model_params_fdr_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  write.table(test@deSV, paste0("model_params_deSV_",method,"_",ms,"_",mf,"_",aso,".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
  
  tmp = get_fdr(test, "Treatment", aso, NULL)
  colnames(tmp) = c("FDR.Treatment", "LFC.Treatment", "name", "class")
  stats = rbind(stats, tmp)
  trt = rbind(trt, get_ci(test,paste0("Treatment",aso), aso))
  
  #norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Treatment"),ntop = 100)
  #ggsave(paste0("pdf/PCA_rnaseq_ASO_",aso,"_effect_corrected.jpg"))
  #new1 = cbind(Exp = unlist(norm1[key_name, ]), grp)
  #ggplot(new1, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio", color = "Donor")
  #ggsave(paste0("pdf/boxplot_rnaseq_IR_ASO_",aso,"_effect_corrected.jpg"))
  
  rm(test, tmp)
  gc(verbose = F)
}

write.table(stats, paste0("res_FDR_LFC_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
write.table(trt, paste0("res_LFC_CI_",method,"_",ms,"_",mf,".txt"), row.names = F, col.names = T, quote = F, sep = "\t")


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
df[df==Inf] = 0.05
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON","dSVA", "880", "131", "307", "877", "876", "881", "879"))

ggplot(data = df, aes(x = class , y = LFC, ymin = LFC-margin, ymax = LFC+margin)) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = df[df$class=="dSVA","LFC"], linetype="dashed", color = "orange") +
  geom_hline(yintercept = df[df$class=="CON","LFC"], linetype="dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  #scale_y_continuous(limits = c(-0.013, 0.005))+
  xlab("") + ylab("IR Change")

ggsave("pdf/errorbar_rnaseq_IRchange_to_untreated_XDP_corrected_new.jpg")


#############################################
#CapSeq
#############################################
#IR
grp = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
#grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))
grp$BatchGrowth = factor(grp$BatchGrowth)
table(paste0(grp$Genotype,grp$ParentLine))

mat = read.table("ASO_Capseq_R1.IR.txt", header = T, row.names = 1, sep = "\t")
mat = mat[,rownames(grp)]
mat = 1000 * mat / get_intron_length(rownames(mat))
norm = mat
#rnaPCA(norm,grp,c("Genotype"),ntop = 5)
new = cbind(Exp = unlist(norm[1, ]), grp)
ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG()
ggsave("pdf/boxplot_capseq_IR_XDP_vs_CON_raw.jpg")

mat = rbind(mat, mat)
test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glm", apply_sva = F, test_terms = c("Genotype","BatchGrowth"), remove_terms = c("BatchGrowth"), random_model_str = "(1|ParentLine)", offset_str = NULL, distro_family = "gaussian")

norm1 = test@deSV
#rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
new1 = cbind(Exp = norm1[1, ], grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG()
ggsave("pdf/boxplot_capseq_IR_XDP_vs_CON_corrected.jpg")

bl1 = get_ci(test,"GenotypedSVA", "dSVA")
bl1 = bl1[1,]
bl2 = get_ci(test,"GenotypeCON", "CON")
bl2 = bl2[1,]
bl = rbind(bl1,bl2)


#Treatment
trt = c()
for (aso in c("307", "876", "877", "879", "880", "881")){
  grp = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  #grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
  grp = grp[which(grp$Treatment %in% c("NO", aso)), ]
  grp = grp[which(grp$Genotype %in% c("XDP")), ]
  #grp$Treatment = factor(grp$Treatment,levels=c("NO", "307", "876", "877", "879", "880", "881"))
  grp$Treatment = factor(grp$Treatment,levels=c("NO", aso))
  grp$BatchGrowth = factor(grp$BatchGrowth)
  grp$MIN = factor(grp$MIN)
  table(paste0(grp$Treatment,"_",grp$ParentLine))
  
  mat = read.table("ASO_Capseq_R1.IR.txt", header = T, row.names = 1, sep = "\t")
  mat = mat[,rownames(grp)]
  mat = 1000 * mat / get_intron_length(rownames(mat))
  norm = mat
  #rnaPCA(norm,grp,c("ParentLine"),ntop = 5)
  new = cbind(Exp = unlist(norm[1, ]), grp)
  ggplot(new, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG()
  ggsave(paste0("pdf/boxplot_capseq_IR_ASO_",aso,"_effect_raw.jpg"))
  
  mat = rbind(mat, mat)
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glm", apply_sva = F, test_terms = c("Treatment","BatchTreat"), remove_terms = c("BatchTreat"), random_model_str = "(1|ParentLine)", offset_str = NULL, distro_family = "gaussian")
  
  norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
  new1 = cbind(Exp = unlist(norm1[1, ]), grp)
  ggplot(new1, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG()
  ggsave(paste0("pdf/boxplot_capseq_IR_ASO_",aso,"_effect_corrected.jpg"))
  
  trt = rbind(trt, get_ci(test,paste0("Treatment",aso), aso)[1,])
}

df = rbind(bl,trt)
df[df==Inf] = 0.01
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON","dSVA", "307", "876", "877", "879", "880", "881"))

ggplot(data = df, aes(x = class , y = logFC, ymin = logFC-1.96*margin, ymax = logFC+1.96*margin)) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = df[df$class=="dSVA","logFC"], linetype="dashed", color = "orange") +
  geom_hline(yintercept = df[df$class=="CON","logFC"], linetype="dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(-0.1, 0.1))+
  xlab("") + ylab("IR Change")

ggsave("pdf/errorbar_capseq_IRchange_to_untreated_XDP_corrected.jpg")
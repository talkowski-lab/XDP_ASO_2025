library(ggplot2)
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")


setwd("~/projects/xdp_aso/nsc/output/")


#############################################
#IR
#############################################
#############################################
#Untreated
#############################################
#Baseline
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
grp = read.table("model_params_colData_glm_IR_Genotype.txt", header = T, row.names = 1, sep = "\t")
grp$MIN = factor(grp$MIN)
grp[which(grp$Genotype1 == "XDP_noEdit"), "Genotype1"] = "XDP_unedited"
grp[which(grp$Genotype1 == "XDP"), "Genotype1"] = "XDP_unexposed" 
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))

mat = read.table("model_params_deSV_glm_IR_Genotype.txt", header = T, row.names = 1, sep = "\t")
rnaPCA(mat,grp,c("Genotype"),ntop = 500)
ggsave("pdf/PCA_rnaseq_IR_XDP_vs_CON_corrected_new.jpg")

new = cbind(Exp = unlist(mat[key_name, ]), grp)
ggplot(new, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "IR Ratio per Kb", color = "Donor")
ggsave("pdf/boxplot_rnaseq_IR_XDP_vs_CON_corrected_new.jpg")



# Rescue
bl1 = read.table("res_FDR_LFC_glm_IR_Genotype_CON_vs_XDP.txt", header = T, sep = "\t")
rownames(bl1) = as.vector(bl1$name)
bl2 = read.table("res_FDR_LFC_glm_IR_Genotype_dSVA_vs_XDP.txt", header = T, sep = "\t")
rownames(bl2) = as.vector(bl2$name)
bl3 = read.table("res_FDR_LFC_glmGaussianZsum_IR_Genotype_dSVA_vs_CON.txt", header = T, sep = "\t")
rownames(bl3) = as.vector(bl3$name)
colnames(bl3) = paste0(colnames(bl3),"1")
shared_genes = intersect(intersect(rownames(bl1), rownames(bl2)), rownames(bl3))
shared_genes = shared_genes[grep("clean",shared_genes)]
new = cbind(bl1[shared_genes, 1:2], bl2[shared_genes, 1:2], bl3[shared_genes, 1:2])
new=as.matrix(new)
new[which(is.na(new))] = 1
new=as.data.frame(new)
new$status = "Other"
new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.GenotypedSVA < 0.1) & (sign(new$LFC.GenotypeCON) == sign(new$LFC.GenotypedSVA)), "status"] = "Rescued"
new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.GenotypedSVA < 0.1) & (sign(new$LFC.GenotypeCON) != sign(new$LFC.GenotypedSVA)), "status"] = "Un-rescued"
new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.GenotypedSVA >= 0.1), "status"] = "Un-rescued"
new[(new$FDR.GenotypedSVA1 < 0.1) & (new$status != "Rescued"), "status"] = "Off-target"
table(new$status)
t0 = rownames(new)[new$status == "Rescued"]
t1 = rownames(new)[new$status == "Un-rescued"]
t2 = rownames(new)[new$status == "Off-target"]
t3 = rownames(new)
write.table(as.data.frame(unique(t0)), file = paste0("dSVA_Rescued_IR_IRFinderID.txt"), row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(t1)), file = paste0("dSVA_Unrescued_IR_IRFinderID.txt"), row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(t2)), file = paste0("dSVA_Offtarget_IR_IRFinderID.txt"), row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(sapply(strsplit(t0,"/"),"[[",4))), file = paste0("dSVA_Rescued_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(sapply(strsplit(t1,"/"),"[[",4))), file = paste0("dSVA_Unrescued_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(sapply(strsplit(t2,"/"),"[[",4))), file = paste0("dSVA_Offtarget_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
write.table(as.data.frame(unique(sapply(strsplit(t3,"/"),"[[",4))), file = paste0("dSVA_background_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)


bl = read.table("res_FDR_LFC_glm_IR_Genotype.txt", header = T, sep = "\t")
bl$val0 = -log10(bl$FDR.GenotypeCON) * sign(bl$LFC.GenotypeCON)
bl$val1 = -log10(bl$FDR.GenotypedSVA) * sign(bl$LFC.GenotypedSVA)
bl$color = "Other IR events"
bl[bl$name == key_name, "color"] = "TAF1 intron 32"
bl = bl[order(bl$val0, decreasing = F), ]
bl = bl[grep("clean", bl$name), ]

ggplot(bl, aes(x = val0, y = val1, color = color)) + geom_point() + 
  xlim(c(-20,20)) + ylim(c(-20,20)) + 
  geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = c("red","black","red")) + 
  geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed", color = c("red","black","red")) + 
  xlab("Significance of IR Changes in Control (directional)") + 
  ylab("Significance of IR Changes in dSVA (directional)") +
  geom_color(c("blue","red")) + geom_noBG()
ggsave("pdf/scatter_IR_change_untreated_genotypes_new.jpg")


#############################################
#Treated
#############################################
key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
trt = read.table(paste0("res_LFC_CI_glm_IR_Treatment.txt"), header = T, sep = "\t")
trt0 = trt[trt$name == key_name, ]
geno = read.table(paste0("res_LFC_CI_glm_IR_Genotype.txt"), header = T, sep = "\t")
geno0 = geno[geno$name == key_name, ]

df = rbind(geno0, trt0)
df[df==Inf] = 0.05
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON","dSVA", "880", "131", "307", "877", "876", "881", "879"))

ggplot(data = df, aes(x = class , y = LFC, ymin = LFC-margin, ymax = LFC+margin)) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = df[df$class=="dSVA","LFC"], linetype="dashed", color = "orange") +
  geom_hline(yintercept = df[df$class=="CON","LFC"], linetype="dashed", color = "blue") +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
  xlab("") + ylab("IR Change") + 
  geom_noBG()

ggsave("pdf/errorbar_rnaseq_IRchange_to_untreated_XDP_corrected_new.jpg")


stat = read.table("res_FDR_LFC_glm_IR_Treatment.txt", header = T, sep = "\t")
for (aso in  c("880", "131", "307", "877", "876", "881", "879")){
  grp = read.table(paste0("model_params_colData_glm_IR_Treatment_",aso,".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
  grp$MIN = factor(grp$MIN)
  mat = read.table(paste0("model_params_deSV_glm_IR_Treatment_",aso,".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
  rnaPCA(as.matrix(mat),grp,c("Treatment"),ntop = 500)
  ggsave(paste0("pdf/PCA_rnaseq_IR_untreated_XDP_vs_ASO", aso, "_new.jpg"))
  
  bl1 = read.table(paste0("res_FDR_LFC_glm_IR_Genotype_CON_vs_XDP.txt"), header = T, sep = "\t")
  rownames(bl1) = as.vector(bl1$name)
  bl2 = stat[which(stat$class == aso),]
  rownames(bl2) = as.vector(bl2$name)
  bl3 = read.table(paste0("res_FDR_LFC_crossModel_IR_Treatment",aso,"_vs_GenotypeCON.txt"), header = T, sep = "\t")
  rownames(bl3) = as.vector(bl3$name)
  shared_genes = intersect(intersect(rownames(bl1), rownames(bl2)), rownames(bl3))
  shared_genes = shared_genes[grep("clean",shared_genes)]
  new = cbind(bl1[shared_genes, 1:2], bl2[shared_genes, 1:2], bl3[shared_genes, 1:2])
  new=as.matrix(new)
  new[which(is.na(new))] = 1
  new=as.data.frame(new)
  new$status = "Other"
  new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.Treatment < 0.1) & (sign(new$LFC.GenotypeCON) == sign(new$LFC.Treatment)), "status"] = "Rescued"
  new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.Treatment < 0.1) & (sign(new$LFC.GenotypeCON) != sign(new$LFC.Treatment)), "status"] = "Un-rescued"
  new[(new$FDR.GenotypeCON < 0.1) & (new$FDR.Treatment >= 0.1), "status"] = "Un-rescued"
  new[(new[[paste0("FDR.Treatment",aso)]] < 0.1) & (new$status != "Rescued"), "status"] = "Off-target"
  table(new$status)
  t0 = rownames(new)[new$status == "Rescued"]
  t1 = rownames(new)[new$status == "Un-rescued"]
  t2 = rownames(new)[new$status == "Off-target"]
  t3 = rownames(new)
  write.table(as.data.frame(unique(t0)), file = paste0("ASO",aso,"_Rescued_IR_IRFinderID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(unique(t1)), file = paste0("ASO",aso,"_Unrescued_IR_IRFinderID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(unique(t2)), file = paste0("ASO",aso,"_Offtarget_IR_IRFinderID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(unique(sapply(strsplit(t0,"/"),"[[",4))), file = paste0("ASO",aso,"_Rescued_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(unique(sapply(strsplit(t1,"/"),"[[",4))), file = paste0("ASO",aso,"_Unrescued_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(unique(sapply(strsplit(t2,"/"),"[[",4))), file = paste0("ASO",aso,"_Offtarget_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  write.table(as.data.frame(unique(sapply(strsplit(t3,"/"),"[[",4))), file = paste0("ASO",aso,"_background_IRhost_EnsemblID.txt"), row.names = F, col.names = F, quote = F)
  
  ov = intersect(bl$name, stat$name)
  tt = stat[stat$class==aso, ]
  rownames(tt) = as.vector(tt$name)
  tt = tt[ov, ]
  tt$val0 = -log10(tt$FDR.Treatment) * sign(tt$LFC.Treatment)
  bb = bl[ov, ]
  tt$val1 = -log10(bb$FDR.GenotypedSVA) * sign(bb$LFC.GenotypedSVA)
  
  tt$color = "Other IR events"
  tt[tt$name == key_name, "color"] = "TAF1 intron 32"
  tt = tt[order(tt$val0, decreasing = F), ]
  tt = tt[grep("clean", tt$name), ]
  
  #ggplot(tt, aes(x = LFC.Treatment, y = -log10(FDR.Treatment), color = color)) + geom_point() +
  #  geom_hline(yintercept = 1, linetype = "dashed", color = c("black")) + 
  #  geom_vline(xintercept = 0, linetype = "dashed", color = c("black")) + 
  #  xlab("IR Change per KB") +
  #  geom_color(c("white","red")) + 
  #  geom_noBG()
  #ggsave(paste0("pdf/vocalno_rnaseq_IR_untreated_XDP_vs_ASO", aso, "_new.jpg"))
  
  ggplot(tt, aes(x = val0, y = val1, color = color)) + geom_point() + 
    xlim(c(-20,20)) + ylim(c(-20,20)) + 
    geom_hline(yintercept = c(-1, 0, 1), linetype = "dashed", color = c("red","black","red")) + 
    geom_vline(xintercept = c(-1, 0, 1), linetype = "dashed", color = c("red","black","red")) + 
    xlab(paste0("Significance of IR Changes in ASO", aso, " (directional)")) + 
    ylab("Significance of IR Changes in dSVA (directional)") +
    geom_color(c("blue","red")) + geom_noBG()
  ggsave(paste0("pdf/scatter_IR_change_treated_ASO", aso, "_new.jpg"))
}


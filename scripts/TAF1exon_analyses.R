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
#CapSeq
#############################################
#TAF1-exons
grp = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
#grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))
grp$BatchGrowth = factor(grp$BatchGrowth)
table(paste0(grp$Genotype,grp$ParentLine))

mat = read.table("ASO_Capseq_R1.TAF1exons.txt", header = T, row.names = 1, sep = "\t")
mat = mat[,rownames(grp)]
norm = log1p(t(t(mat)/grp$sizeFactor))
#rnaPCA(norm,grp,c("Genotype"),ntop = 5)
new = c()
for (i in 1:nrow(mat)){
  new = rbind(new,cbind(Exp = unlist(norm[1, ]),Exon = paste0("E",i), grp))
}
new$Exon = factor(new$Exon,levels = paste0("E",1:38))
#ggplot(new, aes(Exon, Exp,color=Genotype)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(col = factor(MIN))) + geom_noBG()
ggplot(new, aes(Exon, Exp,color=Genotype)) + geom_boxplot(outlier.alpha = 0)  + geom_noBG()
ggsave("pdf/boxplot_capseq_TAF1exons_XDP_vs_CON_raw.pdf")

test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glmm", apply_sva = F, test_terms = c("Genotype","BatchGrowth"), remove_terms = c("BatchGrowth"), random_model_str = "(1|ParentLine)")

norm1 = test@deSV
#rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
new1 = cbind(Exp = unlist(norm[1, ]), grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(col = factor(MIN))) + geom_noBG()
ggsave("pdf/boxplot_capseq_TAF1exons_XDP_vs_CON_corrected.pdf")

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
  
  mat = read.table("ASO_Capseq_R1.i32splicing.txt", header = T, row.names = 1, sep = "\t")
  mat = mat[,rownames(grp)]
  norm = log1p(t(t(mat)/grp$sizeFactor))
  #rnaPCA(norm,grp,c("ParentLine"),ntop = 5)
  new = cbind(Exp = unlist(norm[1, ]), grp)
  ggplot(new, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(col = factor(MIN))) + geom_noBG()
  ggsave(paste0("pdf/boxplot_capseq_TAF1exons_ASO_",aso,"_effect_raw.pdf"))
  
  test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glmm", apply_sva = F, test_terms = c("Treatment","BatchTreat"), remove_terms = c("BatchTreat"), random_model_str = "(1|ParentLine)")
  
  norm1 = test@deSV
  #rnaPCA(norm1,grp,c("Genotype"),ntop = 5)
  new1 = cbind(Exp = unlist(norm1[1, ]), grp)
  ggplot(new1, aes(Treatment, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(col = factor(MIN))) + geom_noBG()
  ggsave(paste0("pdf/boxplot_capseq_TAF1exons_ASO_",aso,"_effect_corrected.pdf"))
  
  trt = rbind(trt, get_ci(test,paste0("Treatment",aso), aso)[1,])
}

df = rbind(bl,trt)
df[df==Inf] = 0.01
df = as.data.frame(df)
df$class = factor(df$class,levels = c("CON","dSVA", "307", "876", "877", "879", "880", "881"))

ggplot(data = df, aes(x = class , y = exp(logFC), ymin = exp(logFC-1.96*margin), ymax = exp(logFC+1.96*margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = exp(df[df$class=="dSVA","logFC"]), linetype="dashed", color = "orange") +
  geom_hline(yintercept = exp(df[df$class=="CON","logFC"]), linetype="dashed", color = "blue") +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0.001, 1.2))+
  xlab("") + ylab("Fold Change")

ggsave("pdf/errorbar_capseq_TAF1exons_to_untreated_XDP_corrected.pdf")


#############################################
#RNASeq
#############################################
#TAF1-exons
grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t")
#grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
#grp = grp[which(grp$Genotype %in% c("XDP", "dSVA", "CON")), ]
grp$Genotype = factor(grp$Genotype, levels = c("XDP", "CON", "dSVA"))
#grp$Treatment = factor(grp$Treatment, levels = c("NO", "307", "876", "877", "879", "880", "881"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)

mat = read.table("ASO_RNAseq_R1.TAF1exons.txt", header = T, row.names = 1, sep = "\t", check.names = F)
mat = mat[,rownames(grp)]
norm = log1p(t(t(mat)/grp$sizeFactor))
new = c()
for (i in 1:nrow(mat)){
  new = rbind(new,cbind(Exp = unlist(norm[i, ]),Exon = paste0("E",i), grp))
}
new$Exon = factor(new$Exon,levels = paste0("E",1:38))
ggplot(new, aes(Exon, log2(exp(Exp)),color=Genotype)) + geom_boxplot(outlier.alpha = 0)  + geom_noBG() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
ggsave("pdf/boxplot_rnaseq_TAF1exons_XDP_vs_CON_raw.pdf", width = 21)

test = buildDEmodel(feature_matrix = mat, sampleTable = grp, method = "glm", apply_sva = F, test_terms = c("Genotype","BatchGrowth"), remove_terms = c("BatchGrowth"), random_model_str = "(1|ParentLine)", offset_str = "sizeFactor", distro_family = "poisson")

norm1 = log1p(exp(test@deSV))
new1 = c()
for (i in 1:nrow(norm1)){
  new1 = rbind(new1,cbind(Exp = unlist(norm1[i, ]),Exon = paste0("E",i), grp))
}
new1$Exon = factor(new1$Exon,levels = paste0("E",1:38))
ggplot(new1, aes(Exon, log2(exp(Exp)),color=Genotype)) + geom_boxplot(outlier.alpha = 0)  + geom_noBG() + theme(axis.text.x = element_text(angle = 45, hjust=1)) + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
ggsave("pdf/boxplot_rnaseq_TAF1exons_XDP_vs_CON_corrected.pdf", width = 21)

tt = get_fdr(test,"Genotype","test")


bl1 = get_ci(test,"GenotypedSVA", "dSVA")
bl1$exon = rownames(bl1)
bl2 = get_ci(test,"GenotypeCON", "CON")
bl2$exon = rownames(bl2)

bl = as.data.frame(rbind(bl1,bl2))


#Treatment
trt = c()
stats = c()
for (aso in c("131", "307", "876", "877", "879", "880", "881")){
  grp = read.table("ASO_RNAseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  #grp = grp[which(grp$Treatment %in% c("NO", "307", "876", "877", "879", "880", "881")), ]
  grp = grp[which(grp$Treatment %in% c("NO", aso)), ]
  grp = grp[which(grp$Genotype %in% c("XDP")), ]
  #grp$Treatment = factor(grp$Treatment,levels=c("NO", "307", "876", "877", "879", "880", "881"))
  grp$Treatment = factor(grp$Treatment,levels=c("NO", aso))
  grp$BatchGrowth = factor(grp$BatchGrowth)
  grp$MIN = factor(grp$MIN)
  table(paste0(grp$Treatment,"_",grp$ParentLine))
  
  mat = read.table("ASO_RNAseq_R1.TAF1exons.txt", header = T, row.names = 1, sep = "\t", check.names = F)
  #mat = round(mat /10)
  mat = mat[,rownames(grp)]
  norm = log1p(t(t(mat)/grp$sizeFactor))
  new = c()
  for (i in 1:nrow(mat)){
    new = rbind(new,cbind(Exp = unlist(norm[i, ]),Exon = paste0("E",i), grp))
  }
  new$Exon = factor(new$Exon,levels = paste0("E",1:38))
  #ggplot(new, aes(Exon, Exp,color=Genotype)) + geom_boxplot(outlier.alpha = 0) + geom_jitter(aes(col = factor(MIN))) + geom_noBG()
  ggplot(new, aes(Exon, log2(exp(Exp)),color=Treatment)) + geom_boxplot(outlier.alpha = 0)  + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
  ggsave(paste0("pdf/boxplot_rnaseq_TAF1exons_ASO_",aso,"_effect_raw.pdf"))
  
  test = buildDEmodel(mat, grp, "glmm", F, c("Treatment","BatchTreat"), c("BatchTreat"), "(1|ParentLine)")
  
  norm1 = test@deSV
  new1 = c()
  for (i in 1:nrow(norm1)){
    new1 = rbind(new1,cbind(Exp = unlist(norm1[i, ]),Exon = paste0("E",i), grp))
  }
  new1$Exon = factor(new1$Exon,levels = paste0("E",1:38))
  ggplot(new1, aes(Exon, log2(exp(Exp)),color=Treatment)) + geom_boxplot(outlier.alpha = 0)  + geom_noBG() + labs(x = "", y = "Normalized Expression (log2-scale)", color = "Donor")
  ggsave(paste0("pdf/boxplot_rnaseq_TAF1exons_ASO_",aso,"_effect_corrected.pdf"))
  
  stats = rbind(stats, get_fdr(test, "Treatment", aso))
  tmp = get_ci(test,paste0("Treatment",aso), aso)
  tmp$exon = rownames(tmp)
  trt = rbind(trt, tmp)
}

df = rbind(bl,trt)
df[df$margin>20,"margin"] = 0.5
df = as.data.frame(df)

df0 = df
df0 = df0[(df0$class %in% c("880","131")), ]
df = df[df$exon %in% paste0("TAF1:0",29:38),]
df$class = factor(df$class,levels = c("CON","dSVA", "880","131", "881", "877", "876", "879", "307"))

ggplot(data = df, aes(x = exon , y = exp(logFC), ymin = exp(logFC-margin), ymax = exp(logFC+margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0.6, 1.8))+
  xlab("") + ylab("Fold Change") + facet_wrap(~class,nrow=3) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_x_discrete(labels=c("TAF1:029" = "E29",
                            "TAF1:030" = "E30", 
                            "TAF1:031" = "E31", 
                            "TAF1:032" = "E32", 
                            "TAF1:033" = "E33", 
                            "TAF1:034" = "E34", 
                            "TAF1:035" = "E35", 
                            "TAF1:036" = "E36", 
                            "TAF1:037" = "E37", 
                            "TAF1:038" = "E38"))

ggsave("pdf/errorbar_rnaseq_TAF1exons29to38_to_untreated_XDP_corrected.pdf")


ggplot(data = df0, aes(x = exon , y = exp(logFC), ymin = exp(logFC-margin), ymax = exp(logFC+margin))) +
  geom_point(size=1,col = "red") +
  geom_errorbar(size=0.5, width = 0.3) +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") +
  theme_bw()+
  theme(axis.line.x = element_line(color = 'black'),axis.line.y = element_line(color = 'black'),axis.text.x = element_text(angle = 0, hjust = 0.5),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())+
  scale_y_continuous(limits = c(0.6, 1.8))+
  xlab("") + ylab("Fold Change") + facet_wrap(~class,nrow=2) + theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave("pdf/errorbar_rnaseq_TAF1exons_to_untreated_XDP_corrected_ASO_131_880.pdf")

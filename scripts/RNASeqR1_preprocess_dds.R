library(DESeq2)
library(matrixStats)
library(ggplot2)
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")


setwd("~/projects/xdp_aso/nsc/output/")

load("RData/RNAseq_ASO_dSVA_all.Rdata")
df0 = as.data.frame(colData(desTissue))
df0 = df0[,c("Genotype","Treatment","Parent_Line_CV","Clone_CV","Batch_Growing","Batch_Treatment")]
colnames(df0) = c("geno","Treatment","ParentLine","Clone","BatchGrowth","BatchTreat")
df0$Genotype = df0$geno
df0[which(df0$geno == "XDP_non_edit"),"Genotype"] = "XDP"
df0[which(df0$geno == "XDP_dSVA"),"Genotype"] = "dSVA"
df0$Treatment = gsub("ASO_","",df0$Treatment)
df0$Treatment = gsub("_ASO","",df0$Treatment)

table(df0$ParentLine)
df0$MIN = "32517"
df0[which(df0$ParentLine == "33109_2B" | df0$ParentLine == "33109_2G"),"MIN"] = "33109"
df0[which(df0$ParentLine == "33113_2D" | df0$ParentLine == "331132I"),"MIN"] = "33113"
df0[which(df0$ParentLine == "33114B" | df0$ParentLine == "33114C"),"MIN"] = "33114"
df0[which(df0$ParentLine == "33362C" | df0$ParentLine == "33362D"),"MIN"] = "33362"
df0[which(df0$ParentLine == "33363C" | df0$ParentLine == "33363D"),"MIN"] = "33363"
df0[which(df0$ParentLine == "33808D"),"MIN"] = "33808"
df0[which(df0$ParentLine == "34363A" | df0$ParentLine == "34363B"),"MIN"] = "34363"
df0[which(df0$ParentLine == "35326I"),"MIN"] = "35326"
df0[which(df0$ParentLine == "35613B"),"MIN"] = "35613"
df0[which(df0$ParentLine == "35833A"),"MIN"] = "35833"
dds = desTissue
colData(dds) = colData(dds)[,1,drop=F]
colData(dds) = cbind(colData(dds), df0)
saveRDS(dds,"RData/ASO_RNAseq_R1.dds_all_exp.rds")


dds = readRDS("RData/ASO_RNAseq_R1.dds_all_exp.rds")
cnt = counts(dds,normalized = F)
write.table(cnt, file="ASO_RNAseq_R1.exp.txt", row.names=T, col.names=T, quote=F, sep="\t")
write.table(df0, file="ASO_RNAseq_R1.metadata.simple.txt", row.names=T, col.names=T, quote=F, sep="\t")

exp = which(rowSums(counts(dds, normalized = F)) >= 5*ncol(dds))
exp_ids = rownames(dds)[exp]
dds = dds[exp_ids, ]

cnt = counts(dds,normalized = F)
write.table(cnt, file="ASO_RNAseq_R1.exp.filtered.txt", row.names=T, col.names=T, quote=F, sep="\t")

dds = estimateSizeFactors(dds)
df0 = as.data.frame(colData(dds))
df0 = df0[,-1]
write.table(df0, file="ASO_RNAseq_R1.metadata.final.txt", row.names=T, col.names=T, quote=F, sep="\t")


load("RData/RNAseq_ASO_dSVA_Exon.Rdata")
desTissue_Exon = desTissue_Exon[1:38, ]
dds1 = desTissue_Exon
colData(dds1) = colData(dds1)[,1,drop=F]
colData(dds1) = cbind(colData(dds1), df0)
saveRDS(dds1,"RData/ASO_RNAseq_R1.dds_all_TAF1exon.rds")

cnt1 = counts(dds1,normalized = F)
write.table(cnt1, file="ASO_RNAseq_R1.TAF1exons.txt", row.names=T, col.names=T, quote=F, sep="\t")


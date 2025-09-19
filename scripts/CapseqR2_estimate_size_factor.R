library(DESeq2)
library(matrixStats)
library(ggplot2)
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")


setwd("~/projects/xdp_aso/nsc/output/")


ref = read.table("GRCh37.75_ensemblID_geneSymbol_pairs.txt", row.names = 1)
colnames(ref) = "gs"

raw1 = read.table("ASO_Capseq_R1.exp.txt", header = T, sep = "\t", row.names = 1)
raw2 = read.table("ASO_Capseq_R2.exp.txt", header = T, sep = "\t", row.names = 1)
colnames(raw1) = paste0(colnames(raw1), "_R1")
colnames(raw2) = paste0(colnames(raw2), "_R2")

meta1 = read.table("ASO_Capseq_R1.metadata.simple.txt", header = T, sep = "\t")
meta2 = read.table("ASO_Capseq_R2.metadata.simple.txt", header = T, sep = "\t")
rownames(meta1) = paste0(rownames(meta1), "_R1")
rownames(meta2) = paste0(rownames(meta2), "_R2")

raw1 = raw1[,rownames(meta1)]
raw2 = raw2[,rownames(meta2)]
raw = cbind(raw1, raw2)

meta = rbind(meta1[, c(-5, -6)], meta2[, c(-5, -6)])

dds0 = DESeqDataSetFromMatrix(countData = raw, colData = meta, design = ~1)
exp = which(rowSums(counts(dds0, normalized = F)) >= 1.5*10^7)
exp_ids = rownames(dds0)[exp]
exp_genes = ref[exp_ids,"gs"]
sort(exp_genes)
dds = dds0[exp_ids,]
dds = estimateSizeFactors(dds)
df = as.data.frame(colData(dds))
df$Treated = "Yes"
df[which(df$Treatment=="NO"),"Treated"] = "No"

df$Treated = "Yes"
df[which(df$Treatment=="NO"),"Treated"] = "No"
write.table(df, file="ASO_Capseq_R1_2.metadata.final.txt", row.names=T, col.names=T, quote=F, sep="\t")

cnt = counts(dds, normalized = T)
write.table(df, file="ASO_Capseq_R1_2.exp.txt", row.names=T, col.names=T, quote=F, sep="\t")

library(DESeq2)
library(matrixStats)
library(ggplot2)
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")


setwd("~/projects/xdp_aso/nsc/output/")

ref = read.table("GRCh37.75_ensemblID_geneSymbol_pairs.txt", row.names = 1)
colnames(ref) = "gs"
raw = read.table("ASO_Capseq_R1.exp.txt", header = T, sep = "\t", row.names = 1)
meta = read.table("ASO_Capseq_R1.metadata.simple.txt", header = T, sep = "\t")
raw = raw[,rownames(meta)]
dds0 = DESeqDataSetFromMatrix(countData = raw, colData = meta, design = ~1)
exp = which(rowSums(counts(dds0, normalized = F)) >= 1*10^7)
exp_ids = rownames(dds0)[exp]
exp_genes = ref[exp_ids,"gs"]
sort(exp_genes)
dds = dds0[exp_ids,]
dds = estimateSizeFactors(dds)
df = as.data.frame(colData(dds))
df$Treated = "Yes"
df[which(df$Treatment=="NO"),"Treated"] = "No"
write.table(df, file="ASO_Capseq_R1.metadata.final.txt", row.names=T, col.names=T, quote=F, sep="\t")
meta1 = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, sep = "\t")
saveRDS(dds, "RData/ASO_Capseq_R1.dds_all.rds")

mat = raw[exp_ids,]
write.table(mat, file="ASO_Capseq_R1.exp.filtered.txt", row.names=T, col.names=T, quote=F, sep="\t")

rnaPCA(counts(dds,normalized = T), meta1, intgroup = c("Genotype","Treated"), ntop = 4)

cnt = counts(dds,normalized = T)
rownames(cnt) = exp_genes
mat = data.frame(rowSums(cnt), exp_genes)
colnames(mat) = c("counts", "gene")
mat$gene = factor(mat$gene, levels = c("ZMYM3", "NONO", "TAF1", "OGT", "PABPN1"))
ggplot(mat, aes(gene,counts)) + geom_col() + geom_noBG() + ylab("Normalized Expression") + xlab("")
ggsave("pdf/Capture_region_expression.pdf")


a = dds[,dds$LibType == "rRNAdep"]
cnt = counts(a,normalized = T)
rownames(cnt) = exp_genes
mat = data.frame(rowSums(cnt), exp_genes)
colnames(mat) = c("counts", "gene")
mat$gene = factor(mat$gene, levels = c("ZMYM3", "NONO", "TAF1", "OGT", "PABPN1"))
ggplot(mat, aes(gene,counts)) + geom_col() + geom_noBG() + ylab("Normalized Expression") + xlab("")
ggsave("pdf/Capture_region_expression_rRNAdep.pdf")


a = dds[,dds$LibType == "ILL"]
cnt = counts(a,normalized = T)
rownames(cnt) = exp_genes
mat = data.frame(rowSums(cnt), exp_genes)
colnames(mat) = c("counts", "gene")
mat$gene = factor(mat$gene, levels = c("ZMYM3", "NONO", "TAF1", "OGT", "PABPN1"))
ggplot(mat, aes(gene,counts)) + geom_col() + geom_noBG() + ylab("Normalized Expression") + xlab("")
ggsave("pdf/Capture_region_expression_ILL.pdf")

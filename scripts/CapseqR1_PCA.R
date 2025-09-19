library(DESeq2)
library(matrixStats)
library(ggplot2)
source("~/lib/R/geom_noBG.R")


setwd("~/projects/xdp_aso/nsc/output/")

ref = read.table("GRCh37.75_ensemblID_geneSymbol_pairs.txt", row.names = 1)
colnames(ref) = "gs"
raw = read.table("ASO_Capseq_R1.exp.txt", header = T, sep = "\t", row.names = 1)
meta = read.table("ASO_Capseq_R1.metadata.simple.txt", header = T, sep = "\t")
raw = raw[,rownames(meta)]
dds0 = DESeqDataSetFromMatrix(countData = raw, colData = meta, design = ~1)
exp = which(rowSums(counts(dds0, normalized = F)) >= 1*10^8)
exp_ids = rownames(dds0)[exp]
exp_genes = ref[exp_ids,"gs"]
sort(exp_genes)
dds = dds0[exp_ids,]
dds = estimateSizeFactors(dds)

grp = read.table("ASO_Capseq_R1.metadata.final.txt", header = T, row.names = 1, sep = "\t", check.names = F)
grp$sizeFactor = dds$sizeFactor
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_dSVA", "Genotype1"] = "dSVA"
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_noEdit"
grp$Genotype1 = factor(grp$Genotype1, levels = c("XDP", "XDP_noEdit", "CON", "dSVA"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)


mat = read.table("ASO_Capseq_R1.exp.filtered.txt", header = T, row.names = 1, sep = "\t")
mat = mat[rownames(dds),]
mat = mat[,rownames(grp)]
mat = mat[, grp$Treatment!="MAL"]
grp = grp[grp$Treatment!="MAL",]
norm = log1p(t(t(mat)/as.vector(grp$sizeFactor)))
tt = rnaPCA(norm,grp,c("Genotype"),ntop = 5,returnData = T)

mat = mat[,rownames(grp)]
rv <- rowVars(mat)
select <- order(rv, decreasing = TRUE)[seq_len(min(5, length(rv)))]
pca <- prcomp(t(mat[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

ggplot(data = tt, aes_string(x = paste0("PC1"), y = paste0("PC2"), color = "BatchGrowth")) + geom_point(size = 2) + xlab(paste0("PC1",": ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2",": ", round(percentVar[2] * 100), "% variance")) + labs(color = "BatchGrowth")
ggplot(data = tt, aes_string(x = paste0("PC1"), y = paste0("PC2"), color = "BatchTreat")) + geom_point(size = 2) + xlab(paste0("PC1",": ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2",": ", round(percentVar[2] * 100), "% variance")) + labs(color = "BatchTreat")
ggplot(data = tt, aes_string(x = paste0("PC1"), y = paste0("PC2"), color = "Genotype1")) + geom_point(size = 2) + xlab(paste0("PC1",": ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2",": ", round(percentVar[2] * 100), "% variance")) + labs(color = "Genotype")
ggplot(data = tt, aes_string(x = paste0("PC1"), y = paste0("PC2"), color = "Treatment")) + geom_point(size = 2) + xlab(paste0("PC1",": ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2",": ", round(percentVar[2] * 100), "% variance")) + labs(color = "Treatment")

library(ggplot2)


setwd("~/projects/xdp_aso/nsc/output/")


###############################################Merge Gene Expression Begins###############################################
matrix_original = read.table("ASO_RNAseq_R1.exp.txt", header = T, check.names = F) 
matrix_new = read.table("ASO_RNAseq_mergedSamples.exp.txt", header = T, check.names = F)

matrix_new = matrix_new[rownames(matrix_original), ]

for (sp in colnames(matrix_new)){
    matrix_original[[sp]] = matrix_new[[sp]]
}

write.table(matrix_original, file = "ASO_RNAseq_new.exp.txt", sep = "\t", row.names = T, col.names = T, quote = F)
##############################################Merge Gene Expression Ends################################################


########################################################################################################################
#TWO WAYS TO MERGE IR
#ONE TO USE THE IR RESULTS FROM THE MERGED ALIGNMENTS
#ONE TO TAKE THE AVERAGE IR FOR OVERLAPPED SAMPLES
########################################################################################################################
###################################################WAY ONE BEGINS#######################################################
ir_original = read.table("ASO_RNAseq_R1.IR.txt", header = T, check.names = F) 
ir_new = read.table("ASO_RNAseq_mergedSamples.IR.txt", header = T, check.names = F)

all(rownames(ir_original) == rownames(ir_new))

key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
idx = which(rownames(ir_original) == key_name)
il = get_intron_length(key_name)

for (sp in colnames(ir_new)){
  ir_original[[sp]] = ir_new[[sp]]
}

write.table(ir_original, file = "ASO_RNAseq_new.IR.txt", sep = "\t", row.names = T, col.names = T, quote = F)
####################################################WAY ONE ENDS########################################################


###################################################WAY TWO BEGINS#######################################################
ir_r1 = read.table("ASO_RNAseq_R1.IR.txt", header = T, check.names = F) 
ir_r2 = read.table("ASO_RNAseq_redo.IR.txt", header = T, check.names = F)

all(rownames(ir_r1) == rownames(ir_r2))

key_name = "X/70644088/70676646/TAF1/ENSG00000147133/clean/0/+"
idx = which(rownames(ir_r1) == key_name)
#il = get_intron_length(key_name)

ir_new = ir_r1
for (sp in colnames(ir_r2)){
  if (sp %in% colnames(ir_r1)){
    ir_new[[sp]] = (ir_r1[[sp]] + ir_r2[[sp]])/2
  } else {
    ir_new[[sp]] = ir_r2[[sp]]
  }
}
write.table(ir_new, file = "ASO_RNAseq_avg.IR.txt", sep = "\t", row.names = T, col.names = T, quote = F)

#Fake Evaluation
meta0 = read.table("ASO_RNAseq_new.569_namecorrected.metadata.tab", header = T, check.names = F, sep = "\t", row.names = 1)
grp = meta0[colnames(ir_new), ]
grp = grp[which(grp$Treatment %in% c("NO")), ]
grp = grp[which(grp$geno %in% c("XDP_non_edit", "XDP", "CON")), ]
grp$Genotype = factor(grp$Genotype, levels = c("CON", "XDP"))
grp$Genotype1 = grp$geno
grp[grp$geno=="XDP_non_edit", "Genotype1"] = "XDP_unedited"
grp[grp$geno=="XDP", "Genotype1"] = "XDP_unexposed"
grp$Genotype1 = factor(grp$Genotype1, levels = c("CON", "XDP_unexposed", "XDP_unedited"))
grp$BatchGrowth = factor(grp$BatchGrowth)
grp$BatchTreat = factor(grp$BatchTreat)

i32ir = unname(unlist(ir_new[idx, rownames(grp)]))
il = 70676646 - 70644088 + 1
grp$IR = 1000*i32ir/il

ggplot(grp, aes(x = Genotype, y =IR)) + geom_boxplot()
####################################################WAY TWO ENDS########################################################
library(dplyr)


setwd("~/projects/xdp_aso/nsc/output/")

df = read.table("ASO_Capseq_R1.metadata.txt", header = T, sep = "\t", comment.char = "")
sp = scan("ASO_Capseq_R1.samples.txt",what = "character")
ban = scan("ASO_Capseq_R1.banned.txt",what = "character")
sp_final = sp[!(sp %in% ban)]
df0 = df[which(df$newFileNameFinal %in% sp_final),]
rownames(df0) = as.vector(df0$newFileNameFinal)
df0 = df0[,c("Genotype","ASO_treatment","Parent.Line","Clone.2","Batch_Growing","Batch_Treatment","LibType_abbrev","Sample_num_within_treatment")]
colnames(df0) = c("geno","Treatment","ParentLine","Clone","BatchGrowth","BatchTreat","LibType","Sample_num_within_treatment")
df0$geno = gsub(" ","_",df0$geno)
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

write.table(df0, file="ASO_Capseq_R1.metadata.simple.txt", row.names=T, col.names=T, quote=F, sep="\t")





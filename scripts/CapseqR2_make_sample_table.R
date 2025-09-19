library(dplyr)


setwd("~/projects/xdp_aso/nsc/output/")

df = read.table("metaData_ASO_Capseq_R2.tab", header = T, sep = "\t", comment.char = "")
sp = df$Sample
ban = c("ASO_126_XDP_dSVA_33109_2G_1A9", "ASO_126_XDP_dSVA_33109_2G_1D3", "ASO_126_XDP_non_edit_33363C2_1B8", "ASO_126_XDP_non_edit_33808D2_1A8",
        "ASO_131_XDP_dSVA_33109_2G_1A9", "ASO_131_XDP_dSVA_33109_2G_1D3", "ASO_131_XDP_non_edit_33363C2_1B8", "ASO_131_XDP_non_edit_33808D2_1A8",
        "ASO_881_XDP_dSVA_34363B_1A10_2")
sp_final = sp[!(sp %in% ban)]
df0 = df[which(df$Sample %in% sp_final),]
rownames(df0) = as.vector(df0$Sample)
df0 = df0[,c("Genotype","ASO_treatment","Parent_Line","CloneFullText","batch","pool","lib_type","Sample_num_within_treatment")]
colnames(df0) = c("geno","Treatment","ParentLine","Clone","BatchSeq","BatchPool","LibType","Sample_num_within_treatment")
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

write.table(df0, file="ASO_Capseq_R2.metadata.simple.txt", row.names=T, col.names=T, quote=F, sep="\t")





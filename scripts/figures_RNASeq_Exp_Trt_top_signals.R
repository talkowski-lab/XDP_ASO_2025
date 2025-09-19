library(lme4)
library(pbkrtest)
library(Hmisc)
library(DESeq2)
library(sva)
library(ggplot2)
library(ggforce)
library(stringr)
source("~/lib/R/countsFilter.R")
source("~/lib/R/auto_deSV.R")
source("~/lib/R/rnaPCA.R")
source("~/lib/R/geom_noBG.R")
source("~/lib/R/goplot.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/glmm.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/go.R")
source("/home/dadi/projects/xdp_aso/nsc/scripts/plot.R")


setwd("~/projects/xdp_aso/nsc/output/")


####################################################################################################################################################################################
#This study
#Determine XDP DEG lists
#############################################
res_xv = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPvanilla_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_xn = read.table("res_FDR_LFC_glm_Exp_Genotype1_XDPnoEdit_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_xc = read.table("res_FDR_LFC_glm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)

rownames(res_xv) = as.vector(res_xv$name)
rownames(res_xn) = as.vector(res_xn$name)
rownames(res_xc) = as.vector(res_xc$name)

deg_xv = rownames(res_xv)[which(res_xv$FDR.Genotype1XDP_vanilla < 0.1)]
deg_xv_up = rownames(res_xv)[which((res_xv$FDR.Genotype1XDP_vanilla < 0.1) & (res_xv$LFC.Genotype1XDP_vanilla > 0))]
deg_xv_dn = rownames(res_xv)[which((res_xv$FDR.Genotype1XDP_vanilla < 0.1) & (res_xv$LFC.Genotype1XDP_vanilla < 0))]

deg_xv_p = rownames(res_xv)[which(res_xv$pval.Genotype1XDP_vanilla < 0.05)]
deg_xv_up_p = rownames(res_xv)[which((res_xv$pval.Genotype1XDP_vanilla < 0.05) & (res_xv$LFC.Genotype1XDP_vanilla > 0))]
deg_xv_dn_p = rownames(res_xv)[which((res_xv$pval.Genotype1XDP_vanilla < 0.05) & (res_xv$LFC.Genotype1XDP_vanilla < 0))]

deg_xn = rownames(res_xn)[which(res_xn$FDR.Genotype1XDP_noEdit < 0.1)]
deg_xn_up = rownames(res_xn)[which((res_xn$FDR.Genotype1XDP_noEdit < 0.1) & (res_xn$LFC.Genotype1XDP_noEdit > 0))]
deg_xn_dn = rownames(res_xn)[which((res_xn$FDR.Genotype1XDP_noEdit < 0.1) & (res_xn$LFC.Genotype1XDP_noEdit < 0))]

deg_xn_p = rownames(res_xn)[which(res_xn$pval.Genotype1XDP_noEdit < 0.05)]
deg_xn_up_p = rownames(res_xn)[which((res_xn$pval.Genotype1XDP_noEdit < 0.05) & (res_xn$LFC.Genotype1XDP_noEdit > 0))]
deg_xn_dn_p = rownames(res_xn)[which((res_xn$pval.Genotype1XDP_noEdit < 0.05) & (res_xn$LFC.Genotype1XDP_noEdit < 0))]

deg_xc = rownames(res_xc)[which(res_xc$FDR.GenotypeXDP < 0.1)]
deg_xc_up = rownames(res_xc)[which((res_xc$FDR.GenotypeXDP < 0.1) & (res_xc$LFC.GenotypeXDP > 0))]
deg_xc_dn = rownames(res_xc)[which((res_xc$FDR.GenotypeXDP < 0.1) & (res_xc$LFC.GenotypeXDP < 0))]

deg_xc_p = rownames(res_xc)[which(res_xc$pval.GenotypeXDP < 0.05)]
deg_xc_up_p = rownames(res_xc)[which((res_xc$pval.GenotypeXDP < 0.05) & (res_xc$LFC.GenotypeXDP > 0))]
deg_xc_dn_p = rownames(res_xc)[which((res_xc$pval.GenotypeXDP < 0.05) & (res_xc$LFC.GenotypeXDP < 0))]

deg_ov2 = Reduce(intersect, list(deg_xv, deg_xn))
deg_ov2_up = Reduce(intersect, list(deg_xv_up, deg_xn_up))
deg_ov2_dn = Reduce(intersect, list(deg_xv_dn, deg_xn_dn))
deg_ov2_dir = unique(c(deg_ov2_up, deg_ov2_dn))

deg_ov2_p = Reduce(intersect, list(deg_xv_p, deg_xn_p))
deg_ov2_up_p = Reduce(intersect, list(deg_xv_up_p, deg_xn_up_p))
deg_ov2_dn_p = Reduce(intersect, list(deg_xv_dn_p, deg_xn_dn_p))
deg_ov2_dir_p = unique(c(deg_ov2_up_p, deg_ov2_dn_p))
bg_ov2 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn)))

deg_ov3 = Reduce(intersect, list(deg_xv, deg_xn, deg_xc))
deg_ov3_up = Reduce(intersect, list(deg_xv_up, deg_xn_up, deg_xc_up))
deg_ov3_dn = Reduce(intersect, list(deg_xv_dn, deg_xn_dn, deg_xc_dn))
deg_ov3_dir = unique(c(deg_ov3_up, deg_ov3_dn))
bg_ov3 = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc)))
shared_this = Reduce(intersect, list(rownames(res_xv), rownames(res_xn), rownames(res_xc)))

#Permissive Set: XDP_vanilla vs CON at p-value significance
perm = deg_xc_p
perm_up = deg_xc_up_p
perm_dn = deg_xc_dn_p
perm_bg = rownames(res_xc)

#Replication Set: (DEGs of XDP_vanilla vs CON at p-value significance) overlap with (DEGs of XDP_noEdit vs CON at p-value significance)
repl = deg_ov2_dir_p
repl_up = deg_ov2_up_p
repl_dn = deg_ov2_dn_p
repl_bg = bg_ov2

#Stringent Set: (DEGs of XDP_vanilla vs CON at FDR significance) overlap with (DEGs of XDP_noEdit vs CON at FDR significance)
stri = deg_ov2_dir
stri_up = deg_ov2_up
stri_dn = deg_ov2_dn
stri_bg = bg_ov2

perm_top = rownames(res_xc[which((abs(res_xc[["LFC.GenotypeXDP"]]) >= log(1.2)) & (res_xc[["pval.GenotypeXDP"]] < 0.05)), ])


####################################################################################################################################################################################
#dSVA/ASO Signature and Rescue Effects
####################################################################################################################################################################################
####################################################################################################################################################################################
#dSVA/ASO vs CON DEG list and rescued genes
#############################################
#Baseline: ASO-treated XDP vs ASO-treatad CON
aso = "dSVA"
mf = "Treatment"
mf1 = "Trt"
if (aso != "dSVA"){
  aso0 = paste0("ASO", aso)
  aso1 = paste0("XDP", aso)
  sel1 = c("CON", "XDP", "dSVA", aso1, paste0("CON", aso))
  ggcol0 = c("brown", "blue", "orange", "white", "white")
  ggcol1 = c("brown", "blue", "orange", "red", "pink")
  ref_pf = ""
} else {
  aso0 = "dSVA"
  aso1 = "dSVA"
  sel1 = c("CON", "XDP", aso1)
  ggcol0 = c("brown", "blue", "white")
  ggcol1 = c("brown", "blue", "orange")
  ref_pf = ""
}

key_name = "ENSG00000147133"

grp1 = read.table(paste0("model_params_colData_glm_Exp_", mf1, ".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
grp1 = grp1[which(grp1$Trt %in% sel1), ]
grp1$Trt = factor(grp1$Trt, levels = sel1)

mat1 = read.table(paste0("model_params_deSV_glm_Exp_", mf1, ".txt"), header = T, row.names = 1, sep = "\t", check.names = F)
mat1 = mat1[, rownames(grp1)]

rnaPCA(mat1, grp1, c(mf1), ntop = 500) + scale_color_manual(values = ggcol0)
ggsave(paste0("pdf/PCA_rnaseq_Exp_", mf1, "_", aso1, "_vs_CON_corrected_new_white.jpg"), width = 5.7, height = 4.2)
rnaPCA(mat1, grp1, c(mf1), ntop = 500) + scale_color_manual(values = ggcol1)
ggsave(paste0("pdf/PCA_rnaseq_Exp_", mf1, "_", aso1, "_vs_CON_corrected_new.jpg"), width = 5.7, height = 4.2)


new1 = cbind(Exp = unlist(mat1[key_name, ]), grp1)
ggplot(new1, aes(Trt, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Expression (log-scale)")
ggsave(paste0("pdf/boxplot_rnaseq_Exp_", mf1, "_", aso1, "_vs_CON_corrected_new.jpg"), width = 5.7, height = 4.2)

#DEGs of dSVA/ASO vs CON
res_aso = read.table(paste0("res_FDR_LFC_glm_Exp_", mf1, "_", aso1, "_vs_CON", ref_pf, ".txt"), header = T, sep = "\t")
rownames(res_aso) = as.vector(res_aso$name)
res_aso[key_name, ]
bg_aso = rownames(res_aso)
deg_aso_up = rownames(res_aso[which(res_aso$FDR.Treatment < 0.1 & res_aso$LFC.Treatment >= log(1.2)), ])
deg_aso_dn = rownames(res_aso[which(res_aso$FDR.Treatment < 0.1 & res_aso$LFC.Treatment <= -log(1.2)), ])
deg_aso = unique(c(deg_aso_up, deg_aso_dn))
deg_aso_up_p = rownames(res_aso[which(res_aso$pval.Treatment < 0.05 & res_aso$LFC.Treatment >= log(1.2)), ])
deg_aso_dn_p = rownames(res_aso[which(res_aso$pval.Treatment < 0.05 & res_aso$LFC.Treatment <= -log(1.2)), ])
deg_aso_p = unique(c(deg_aso_up_p, deg_aso_dn_p))
deg_aso_nc = rownames(res_aso[which(res_aso$FDR.Treatment >= 0.1), ])
deg_aso_nc_p = rownames(res_aso[which(res_aso$pval.Treatment >= 0.05), ])
length(deg_aso_p)
length(deg_aso_nc_p)
length(deg_aso)
length(deg_aso_nc)

#Number of genes rescued by ASO
perm_top_rescue_aso = intersect(perm_top, deg_aso_nc_p)
length(perm_top_rescue_aso)
repl_rescue_aso = intersect(repl, deg_aso_nc_p)
length(repl_rescue_aso)
stri_rescue_aso = intersect(stri, deg_aso_nc)
length(stri_rescue_aso)

saveRDS(perm_top_rescue_aso, paste0("RData/perm_top_rescue_", aso0, ".rds"))
saveRDS(repl_rescue_aso, paste0("RData/repl_rescue_", aso0, ".rds"))
saveRDS(stri_rescue_aso, paste0("RData/stri_rescue_", aso0, ".rds"))

#Significance of genes rescued by ASO, compared to randomness
overlap_hgtest(perm_top, deg_aso_nc_p, intersect(bg_ov2, bg_aso))
overlap_hgtest(repl, deg_aso_nc_p, intersect(bg_ov2, bg_aso))
overlap_hgtest(stri, deg_aso_nc_p, intersect(bg_ov2, bg_aso))

#make plot to summarize rescue
df = data.frame(count = c(length(perm_top), length(repl), length(stri), length(perm_top_rescue_aso), length(repl_rescue_aso), length(stri_rescue_aso)), 
                category = c("Top Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", "Top Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures"),
                rescue = c("No", "No", "No", "Yes", "Yes", "Yes"),
                label = c(length(perm_top), length(repl), length(stri), round(100*length(perm_top_rescue_aso)/length(perm_top), 2), round(100*length(repl_rescue_aso)/length(repl), 2), round(100*length(stri_rescue_aso)/length(stri), 2)))

df[df$rescue == "No", "label"] = as.character(as.integer(df[df$rescue == "No", "label"]))
df[df$rescue == "Yes", "label"] = paste0(df[df$rescue == "Yes", "label"], "%")
df$class = df$category
df[df$rescue == "Yes", "class"] = paste0(aso0, "-rescued XDP Gene Signatures")
df$category = factor(df$category, levels = c("Top Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures"))
df$class = factor(df$class, levels = c("Top Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", paste0(aso0, "-rescued XDP Gene Signatures")))

ggplot(df, aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 1.5, position = position_dodge(0)) + geom_text(aes(y = count + 50)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave(paste0("pdf/barplot_XDPDEG_3Sets_and_", aso0, "_Rescued.jpg"))

ggplot(df[1:3, ], aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 0.75) + 
  geom_text(aes(y = count + 50)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("pdf/barplot_XDPDEG_3Sets.jpg")

ggplot(df[c(1, 4), ], aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 1, position = position_dodge(0)) + geom_text(aes(y = count + 100)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave(paste0("pdf/barplot_XDPDEG_1Set_and_", aso0, "_Rescued.jpg"), width = 3.7, height = 4.2)

ggplot(df[1, ], aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 0.5) + 
  geom_text(aes(y = count + 100)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("pdf/barplot_XDPDEG_1Set.jpg", width = 3.7, height = 4.2)


####################################################################################################################################################################################
#XDP Gene Signatures in dSVA/ASO-treated vs untreated XDP
####################################################
if (aso == "dSVA"){
  statfile_aso_vs_xdp = "res_FDR_LFC_glmm_Exp_Genotype_dSVA_vs_XDP.txt"
  term_word = "GenotypedSVA"
} else {
  statfile_aso_vs_xdp = paste0("res_FDR_LFC_glm_Exp_", mf1, "_", aso1, "_vs_XDP", ".txt")
  term_word = "Treatment"
}

res_trt = read.table(statfile_aso_vs_xdp, header = T, sep = "\t")
rownames(res_trt) = as.vector(res_trt$name)
res_trt[key_name, ]
bg_trt = rownames(res_trt)
deg_trt_up = rownames(res_trt[which(res_trt[[paste0("FDR.", term_word)]] < 0.1 & res_trt[[paste0("LFC.", term_word)]] > 0), ])
deg_trt_dn = rownames(res_trt[which(res_trt[[paste0("FDR.", term_word)]] < 0.1 & res_trt[[paste0("LFC.", term_word)]] < 0), ])
deg_trt = unique(c(deg_trt_up, deg_trt_dn))
deg_trt_up_p = rownames(res_trt[which(res_trt[[paste0("pval.", term_word)]] < 0.05 & res_trt[[paste0("LFC.", term_word)]] > 0), ])
deg_trt_dn_p = rownames(res_trt[which(res_trt[[paste0("pval.", term_word)]] < 0.05 & res_trt[[paste0("LFC.", term_word)]] < 0), ])
deg_trt_p = unique(c(deg_trt_up_p, deg_trt_dn_p))


#Average change of XDP DEGs vs their Average change against dSVA/ASO
res_xv_perm_top_up = res_xv[perm_top_up, ]
res_trt_perm_top_up = res_trt[intersect(perm_top_up, rownames(res_trt)), ]
res_xc_repl_up = res_xc[repl_up, ]
res_trt_repl_up = res_trt[intersect(repl_up, rownames(res_trt)), ]
res_xc_stri_up = res_xc[stri_up, ]
res_trt_stri_up = res_trt[intersect(stri_up, rownames(res_trt)), ]

res_xv_perm_top_dn = res_xv[perm_top_dn, ]
res_trt_perm_top_dn = res_trt[intersect(perm_top_dn, rownames(res_trt)), ]
res_xc_repl_dn = res_xc[repl_dn, ]
res_trt_repl_dn = res_trt[intersect(repl_dn, rownames(res_trt)), ]
res_xc_stri_dn = res_xc[stri_dn, ]
res_trt_stri_dn = res_trt[intersect(stri_dn, rownames(res_trt)), ]


#Average change of dSVA/ASO-rescued XDP DEGs vs their Average change against dSVA/ASO
perm_top_rescue_aso_up = intersect(perm_top_rescue_aso, perm_top_up)
perm_top_rescue_aso_dn = intersect(perm_top_rescue_aso, perm_top_dn)
repl_rescue_aso_up = intersect(repl_rescue_aso, repl_up)
repl_rescue_aso_dn = intersect(repl_rescue_aso, repl_dn)
stri_rescue_aso_up = intersect(stri_rescue_aso, stri_up)
stri_rescue_aso_dn = intersect(stri_rescue_aso, stri_dn)

res_xv_perm_top_rescue_aso_up = res_xv[perm_top_rescue_aso_up, ]
res_trt_perm_top_rescue_aso_up = res_trt[intersect(perm_top_rescue_aso_up, rownames(res_trt)), ]
res_xc_repl_rescue_aso_up = res_xc[repl_rescue_aso_up, ]
res_trt_repl_rescue_aso_up = res_trt[intersect(repl_rescue_aso_up, rownames(res_trt)), ]
res_xc_stri_rescue_aso_up = res_xc[stri_rescue_aso_up, ]
res_trt_stri_rescue_aso_up = res_trt[intersect(stri_rescue_aso_up, rownames(res_trt)), ]

res_xv_perm_top_rescue_aso_dn = res_xv[perm_top_rescue_aso_dn, ]
res_trt_perm_top_rescue_aso_dn = res_trt[intersect(perm_top_rescue_aso_dn, rownames(res_trt)), ]
res_xc_repl_rescue_aso_dn = res_xc[repl_rescue_aso_dn, ]
res_trt_repl_rescue_aso_dn = res_trt[intersect(repl_rescue_aso_dn, rownames(res_trt)), ]
res_xc_stri_rescue_aso_dn = res_xc[stri_rescue_aso_dn, ]
res_trt_stri_rescue_aso_dn = res_trt[intersect(stri_rescue_aso_dn, rownames(res_trt)), ]


df = data.frame(FC = 100 * exp(res_xc[perm_top_rescue_aso, "LFC.GenotypeXDP"]))
df$class = "Rescued Top Permissive XDP Signatures"
p = ggplot(df, aes(x = class, y = FC)) + 
  geom_sina(aes(color = class), size = 0.1, position = position_dodge(1), maxwidth = 0.5, scale = "area") + 
  geom_hline(yintercept = c(80, 120), color = "red") +
  scale_color_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c( "darkgrey")) +
  xlab("") + 
  ylab("Fold Change (%)") +
  geom_noBG() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  scale_y_continuous(limits = c(0, 200), breaks=seq(0, 200, 20)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1), strip.text = element_text(size = 8), legend.position="none") 
p
ggsave(paste0("pdf/sinaplot_FC_XDPDEG_rescue_1Set_against_", aso0, ".jpg"), width = 2.5, height = 4.2)

fres = rescue_fisher(res_xc[, paste0("LFC.GenotypeXDP"), drop = F], res_trt[,  paste0("LFC.", term_word), drop = F], perm_top_rescue_aso)
fres


####################################################################################################################################################################################
#This study
#ASO vs CON Pathway analyses
#############################################
godb_list = list(
  GO_BP = "~/ref/gsea/gsea_go_bp.rds",
  GO_CC = "~/ref/gsea/gsea_go_cc.rds",
  GO_MF = "~/ref/gsea/gsea_go_mf.rds",
  PATHWAYS = "~/ref/gsea/gsea_pathways.rds",
  REG = "~/ref/gsea/gsea_reg.rds",
  HPO = "~/ref/gsea/gsea_hpo.rds"
)
ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")

perm_top_gs = ref[perm_top, 1]
bg_gs = ref[rownames(res_xc), 1]
perm_top_go = go_test(perm_top_gs, bg_gs, godb_list)
head(perm_top_go[, 1:6], 3)


#deg_aso_gs = ref[deg_aso, 1]
#deg_aso_up_gs = ref[deg_aso_up, 1]
#deg_aso_dn_gs = ref[deg_aso_dn, 1]
bg_aso_gs = ref[bg_aso, 1]

deg_aso_p_gs = ref[deg_aso_p, 1]
#deg_aso_up_p_gs = ref[deg_aso_up_p, 1]
#deg_aso_dn_p_gs = ref[deg_aso_dn_p, 1]

#GO for DEGs at FDR level
#deg_aso_go = go_test(deg_aso_gs, bg_aso_gs, godb_list)
#deg_aso_up_go = go_test(deg_aso_up_gs, bg_aso_gs, godb_list)
#deg_aso_dn_go = go_test(deg_aso_dn_gs, bg_aso_gs, godb_list)
#write.table(deg_aso_go, file = paste0("GSEA_GO_", aso1, "vCON_","FDR_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(deg_aso_up_go, file = paste0("GSEA_GO_", aso1, "vCON_","FDR_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(deg_aso_dn_go, file = paste0("GSEA_GO_", aso1, "vCON_","FDR_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#GO for DEGs at Nominal level
deg_aso_p_go = go_test(deg_aso_p_gs, bg_aso_gs, godb_list)
#deg_aso_up_p_go = go_test(deg_aso_up_p_gs, bg_aso_gs, godb_list)
#deg_aso_dn_p_go = go_test(deg_aso_dn_p_gs, bg_aso_gs, godb_list)
#write.table(deg_aso_p_go, file = paste0("GSEA_GO_", aso1, "vCON_","pval_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(deg_aso_up_p_go, file = paste0("GSEA_GO_", aso1, "vCON_","pval_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
#write.table(deg_aso_dn_p_go, file = paste0("GSEA_GO_", aso1, "vCON_","pval_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#Read in GO results
#deg_aso_p_go = read.table(paste0("GSEA_GO_", aso1, "vCON_","pval_all",".txt"), header = T, row.names = 1, sep = "\t")
#deg_aso_up_p_go = read.table(paste0("GSEA_GO_", aso1, "vCON_","pval_up",".txt"), header = T, row.names = 1, sep = "\t")
#deg_aso_dn_p_go = read.table(paste0("GSEA_GO_", aso1, "vCON_","pval_dn",".txt"), header = T, row.names = 1, sep = "\t")
#deg_aso_go = read.table(paste0("GSEA_GO_", aso1, "vCON_","FDR_all",".txt"), header = T, row.names = 1, sep = "\t")
#deg_aso_up_go = read.table(paste0("GSEA_GO_", aso1, "vCON_","FDR_up",".txt"), header = T, row.names = 1, sep = "\t")
#deg_aso_dn_go = read.table(paste0("GSEA_GO_", aso1, "vCON_","FDR_dn",".txt"), header = T, row.names = 1, sep = "\t")

#force the same order
#deg_aso_p_go = deg_aso_p_go[rownames(deg_aso_p_go), ]
#deg_aso_up_p_go = deg_aso_up_p_go[rownames(deg_aso_p_go), ]
#deg_aso_dn_p_go = deg_aso_dn_p_go[rownames(deg_aso_p_go), ]
#deg_aso_go = deg_aso_go[rownames(deg_aso_p_go), ]
#deg_aso_up_go = deg_aso_up_go[rownames(deg_aso_p_go), ]
#deg_aso_dn_go = deg_aso_dn_go[rownames(deg_aso_p_go), ]


#extract terms of interest
aloy_vec = scan("go_terms_roi_v2.txt", what = "character")
keyword_vec = c(aloy_vec, c("NEURO", "SYNAP"))
deg_aso_p_go_sel = deg_aso_p_go[rownames(perm_top_go), ]
#deg_aso_go_sel = go_extract_keywords(deg_aso_go[which(deg_aso_go$pvalue < 0.05), ], "term", keyword_vec)
#deg_aso_dn_p_go_sel = go_extract_keywords(deg_aso_dn_p_go[which(deg_aso_dn_p_go$pvalue < 0.05), ], "term", keyword_vec)
#deg_aso_dn_go_sel = go_extract_keywords(deg_aso_dn_go[which(deg_aso_dn_go$pvalue < 0.05), ], "term", keyword_vec)
#deg_aso_up_p_go_sel = go_extract_keywords(deg_aso_up_p_go[which(deg_aso_up_p_go$pvalue < 0.05), ], "term", keyword_vec)
#deg_aso_up_go_sel = go_extract_keywords(deg_aso_up_go[which(deg_aso_up_go$pvalue < 0.05), ], "term", keyword_vec)


#non-directional DEGs for GO
#go_plot(deg_aso_p_go_sel[order(deg_aso_p_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
#ggsave(paste0("pdf/GOplot_", aso0, "DEG_","pval_all",".jpg"))
#go_plot(deg_aso_go_sel[order(deg_aso_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
#ggsave(paste0("pdf/GOplot_", aso0, "DEG_","FDR_all",".jpg"))
#intersect(deg_aso_p_go_sel[order(deg_aso_p_go_sel$pvalue), ][1:20, "term"], deg_aso_go_sel[order(deg_aso_go_sel$pvalue), ][1:20, "term"])


#Downregulated DEGs for GO
#go_plot(deg_aso_dn_p_go_sel[order(deg_aso_dn_p_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6)
#ggsave(paste0("pdf/GOplot_", aso0, "DEG_","pval_dn",".jpg"))
#go_plot(deg_aso_dn_go_sel[order(deg_aso_dn_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6)
#ggsave(paste0("pdf/GOplot_", aso0, "DEG_","FDR_dn",".jpg"))
#intersect(deg_aso_dn_p_go_sel[order(deg_aso_dn_p_go_sel$pvalue), ][, "term"], deg_aso_dn_go_sel[order(deg_aso_dn_go_sel$pvalue), ][, "term"])


#Upregulated DEGs for GO
#go_plot(deg_aso_up_p_go_sel[order(deg_aso_up_p_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
#ggsave(paste0("pdf/GOplot_", aso0, "DEG_","pval_up",".jpg"))
#go_plot(deg_aso_up_go_sel[order(deg_aso_up_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
#ggsave(paste0("pdf/GOplot_", aso0, "DEG_","FDR_up",".jpg"))
#intersect(deg_aso_up_p_go_sel[order(deg_aso_up_p_go_sel$pvalue), ][1:20, "term"], deg_aso_up_go_sel[order(deg_aso_up_go_sel$pvalue), ][1:20, "term"])


#Read in GO results for XDP DEGs
#perm_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_all",".txt"), header = T, row.names = 1, sep = "\t")
#perm_up_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_up",".txt"), header = T, row.names = 1, sep = "\t")
#perm_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_dn",".txt"), header = T, row.names = 1, sep = "\t")
#repl_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_all",".txt"), header = T, row.names = 1, sep = "\t")
#repl_up_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_up",".txt"), header = T, row.names = 1, sep = "\t")
#repl_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_dn",".txt"), header = T, row.names = 1, sep = "\t")
#stri_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_all",".txt"), header = T, row.names = 1, sep = "\t")
#stri_up_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_up",".txt"), header = T, row.names = 1, sep = "\t")
#stri_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_dn",".txt"), header = T, row.names = 1, sep = "\t")

#extract terms of interest
#aloy_vec = scan("go_terms_roi_v2.txt", what = "character")
#keyword_vec = c(aloy_vec, c("NEURO", "SYNAP"))
#perm_go_sel = go_extract_keywords(perm_go[which(perm_go$pvalue < 0.05), ], "term", keyword_vec)
#perm_dn_go_sel = go_extract_keywords(perm_dn_go[which(perm_dn_go$pvalue < 0.05), ], "term", keyword_vec)
#perm_up_go_sel = go_extract_keywords(perm_up_go[which(perm_up_go$pvalue < 0.05), ], "term", keyword_vec)
#repl_go_sel = go_extract_keywords(repl_go[which(repl_go$pvalue < 0.05), ], "term", keyword_vec)
#repl_dn_go_sel = go_extract_keywords(repl_dn_go[which(repl_dn_go$pvalue < 0.05), ], "term", keyword_vec)
#repl_up_go_sel = go_extract_keywords(repl_up_go[which(repl_up_go$pvalue < 0.05), ], "term", keyword_vec)
#stri_go_sel = go_extract_keywords(stri_go[which(stri_go$pvalue < 0.05), ], "term", keyword_vec)
#stri_dn_go_sel = go_extract_keywords(stri_dn_go[which(stri_dn_go$pvalue < 0.05), ], "term", keyword_vec)
#stri_up_go_sel = go_extract_keywords(stri_up_go[which(stri_up_go$pvalue < 0.05), ], "term", keyword_vec)


#Comparison of GO Terms from DEGs before and after ASO
perm_top_go_sel = perm_top_go
p = go_change_plot(perm_top_go_sel, deg_aso_p_go, "Permissive XDP Gene Signatures", paste0("Nominal ", aso0, " Gene Signatures"), "query_annotated", "fdr", 0.1, ntop = 20, text_width = 30, text_size = 6, transform_val = F, show_hit = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c(paste0("Nominal ", aso0, " Gene Signatures"), "Permissive XDP Gene Signatures"), 15), values = c("darkgrey", "#0066FF")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Top_Permissive_Set_vs_", aso1, "vCON_Nominal_Set_v1.jpg"), width = 6.1, height = 4.24)

perm_top_go_sel = go_extract_keywords(perm_top_go[which(perm_top_go$pvalue < 0.05), ], "term", keyword_vec)
p = go_change_plot(perm_top_go_sel, deg_aso_p_go, "Permissive XDP Gene Signatures", paste0("Nominal ", aso0, " Gene Signatures"), "query_annotated", "pvalue", 0.05, ntop = 20, text_width = 30, text_size = 6, transform_val = F, show_hit = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c(paste0("Nominal ", aso0, " Gene Signatures"), "Permissive XDP Gene Signatures"), 15), values = c("darkgrey", "#0066FF")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Top_Permissive_Set_vs_", aso1, "vCON_Nominal_Set_v2.jpg"), width = 6.1, height = 4.24)




p = go_change_plot(repl_go_sel, deg_aso_p_go[, 1:15], "Replicated XDP Gene Signatures", paste0("Nominal ", aso0, " Gene Signatures"), "query_annotated", "pvalue", 0.05, ntop = 20, text_width = 30, text_size = 6, transform_val = F, show_hit = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c(paste0("Nominal ", aso0, " Gene Signatures"), "Replicated XDP Gene Signatures"), 15), values = c("darkgrey", "#EC8014")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Replication_Set_vs_", aso1, "vCON_Nominal_Set.jpg"), width = 6.1, height = 4.24)


p = go_change_plot(stri_go_sel, deg_aso_go[, 1:15], "Stringent XDP Gene Signatures", paste0("FDR ", aso0, " Gene Signatures"), "query_annotated", "pvalue", 0.05, ntop = 20, text_width = 30, text_size = 6, transform_val = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c(paste0("FDR ", aso0, " Gene Signatures"), "Stringent XDP Gene Signatures"), 15), values = c("darkgrey", "red")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Stringent_Set_vs_", aso1, "vCON_FDR_Set.jpg"), width = 6.1, height = 4.24)

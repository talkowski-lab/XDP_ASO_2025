library(ggplot2)
library(ggforce)
library(stringr)
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

#Permissive XDP Gene Signatures: XDP_vanilla vs CON at p-value significance
perm = deg_xc_p
perm_up = deg_xc_up_p
perm_dn = deg_xc_dn_p
perm_bg = rownames(res_xc)

#Replicated XDP Gene Signatures: (DEGs of XDP_vanilla vs CON at p-value significance) overlap with (DEGs of XDP_noEdit vs CON at p-value significance)
repl = deg_ov2_dir_p
repl_up = deg_ov2_up_p
repl_dn = deg_ov2_dn_p
repl_bg = bg_ov2

#Stringent XDP Gene Signatures: (DEGs of XDP_vanilla vs CON at FDR significance) overlap with (DEGs of XDP_noEdit vs CON at FDR significance)
stri = deg_ov2_dir
stri_up = deg_ov2_up
stri_dn = deg_ov2_dn
stri_bg = bg_ov2


####################################################################################################################################################################################
#dSVA Signature and Rescue Effects
####################################################################################################################################################################################
####################################################################################################################################################################################
#dSVA DEG list and rescued genes
#############################################
#Baseline: dSVA vs CON 
key_name = "ENSG00000147133"
grp = read.table(paste0("model_params_colData_glm_Exp_Genotype_dSVA_vs_CON.txt"), header = T, row.names = 1, sep = "\t", check.names = F)
grp$MIN = factor(grp$MIN)

mat1 = read.table("model_params_deSV_glm_Exp_Genotype_dSVA_vs_CON.txt", header = T, row.names = 1, sep = "\t")
rnaPCA(mat1, grp, c("Genotype"), ntop = 500)
ggsave("pdf/PCA_rnaseq_Exp_Genotype_CON_vs_dSVA_corrected_new.jpg")

new1 = cbind(Exp = unlist(mat1[key_name, ]), grp)
ggplot(new1, aes(Genotype, Exp)) + geom_boxplot(outlier.alpha = 0) + geom_noBG() + labs(x = "", y = "Expression (log-scale)", color = "Donor")
ggsave("pdf/boxplot_rnaseq_Exp_Genotype_dSVA_vs_CON_corrected_new.jpg")

#DEGs of dSVA vs CON
res_dsva = read.table("res_FDR_LFC_glm_Exp_Genotype_dSVA_vs_CON.txt", header = T, sep = "\t")
rownames(res_dsva) = as.vector(res_dsva$name)
res_dsva[key_name, ]
bg_dsva = rownames(res_dsva)
deg_dsva_up = rownames(res_dsva[which(res_dsva$FDR.GenotypedSVA < 0.1 & res_dsva$LFC.GenotypedSVA > 0), ])
deg_dsva_dn = rownames(res_dsva[which(res_dsva$FDR.GenotypedSVA < 0.1 & res_dsva$LFC.GenotypedSVA < 0), ])
deg_dsva = unique(c(deg_dsva_up, deg_dsva_dn))
deg_dsva_up_p = rownames(res_dsva[which(res_dsva$pval.GenotypedSVA < 0.05 & res_dsva$LFC.GenotypedSVA > 0), ])
deg_dsva_dn_p = rownames(res_dsva[which(res_dsva$pval.GenotypedSVA < 0.05 & res_dsva$LFC.GenotypedSVA < 0), ])
deg_dsva_p = unique(c(deg_dsva_up_p, deg_dsva_dn_p))
deg_dsva_nc = rownames(res_dsva[which(res_dsva$FDR.GenotypedSVA >= 0.1), ])
deg_dsva_nc_p = rownames(res_dsva[which(res_dsva$pval.GenotypedSVA >= 0.05), ])

#Number of genes rescued by dSVA
perm_rescue_dsva = intersect(perm, deg_dsva_nc_p)
length(perm_rescue_dsva)
repl_rescue_dsva = intersect(repl, deg_dsva_nc_p)
length(repl_rescue_dsva)
stri_rescue_dsva = intersect(stri, deg_dsva_nc)
length(stri_rescue_dsva)

saveRDS(perm_rescue_dsva, "RData/perm_rescue_dsva.rds")
saveRDS(repl_rescue_dsva, "RData/repl_rescue_dsva.rds")
saveRDS(stri_rescue_dsva, "RData/stri_rescue_dsva.rds")

#Significance of genes rescued by dSVA, compared to randomness
overlap_hgtest(perm, deg_dsva_nc_p, intersect(bg_ov2, bg_dsva))
overlap_hgtest(repl, deg_dsva_nc_p, intersect(bg_ov2, bg_dsva))
overlap_hgtest(stri, deg_dsva_nc_p, intersect(bg_ov2, bg_dsva))

#make plot to summarize rescue
trt = "dSVA"
mf = "Genotype"
df = data.frame(count = c(length(perm), length(repl), length(stri), length(perm_rescue_dsva), length(repl_rescue_dsva), length(stri_rescue_dsva)), 
                category = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", "Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures"),
                rescue = c("No", "No", "No", "Yes", "Yes", "Yes"),
                label = c(length(perm), length(repl), length(stri), round(100*length(perm_rescue_dsva)/length(perm), 2), round(100*length(repl_rescue_dsva)/length(repl), 2), round(100*length(stri_rescue_dsva)/length(stri), 2)))

df[df$rescue == "No", "label"] = as.character(as.integer(df[df$rescue == "No", "label"]))
df[df$rescue == "Yes", "label"] = paste0(df[df$rescue == "Yes", "label"], "%")
df$class = df$category
df[df$rescue == "Yes", "class"] = paste0(trt, "-rescued XDP Gene Signatures")
df$category = factor(df$category, levels = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures"))
df$class = factor(df$class, levels = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", paste0(trt, "-rescued XDP Gene Signatures")))

ggplot(df, aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 1.5, position = position_dodge(0)) + geom_text(aes(y = count + 50)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("pdf/barplot_XDPDEG_3Sets_and_dSVA_Rescued.jpg")


ggplot(df[1:3, ], aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 0.75) + 
  geom_text(aes(y = count + 50)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("pdf/barplot_XDPDEG_3Sets.jpg")


ggplot(df[c(1, 2, 4, 5), ], aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 1.5, position = position_dodge(0)) + geom_text(aes(y = count + 100)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("pdf/barplot_XDPDEG_2Sets_and_dSVA_Rescued.jpg")


ggplot(df[1:2, ], aes(x = category, y = count, fill = class, label = label)) + 
  geom_col(width = 0.75) + 
  geom_text(aes(y = count + 100)) + 
  labs(x = "", y = "Number of Genes") + geom_noBG() + 
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1))
ggsave("pdf/barplot_XDPDEG_2Sets.jpg")


####################################################################################################################################################################################
#XDP Gene Signatures in dSVA vs XDP
###################################
trt = "dSVA"
mf = "Genotype"
res_trt = read.table("res_FDR_LFC_glmm_Exp_Genotype_dSVA_vs_XDP.txt", header = T, sep = "\t")
rownames(res_trt) = as.vector(res_trt$name)
res_trt[key_name, ]
bg_trt = rownames(res_trt)
deg_trt_up = rownames(res_trt[which(res_trt[[paste0("FDR.", mf, trt)]] < 0.1 & res_trt[[paste0("LFC.", mf, trt)]] > 0), ])
deg_trt_dn = rownames(res_trt[which(res_trt[[paste0("FDR.", mf, trt)]] < 0.1 & res_trt[[paste0("LFC.", mf, trt)]] < 0), ])
deg_trt = unique(c(deg_trt_up, deg_trt_dn))
deg_trt_up_p = rownames(res_trt[which(res_trt[[paste0("pval.", mf, trt)]] < 0.05 & res_trt[[paste0("LFC.", mf, trt)]] > 0), ])
deg_trt_dn_p = rownames(res_trt[which(res_trt[[paste0("pval.", mf, trt)]] < 0.05 & res_trt[[paste0("LFC.", mf, trt)]] < 0), ])
deg_trt_p = unique(c(deg_trt_up_p, deg_trt_dn_p))


#Average change of XDP DEGs vs their Average change against dSVA
res_xv_perm_up = res_xv[perm_up, ]
res_trt_perm_up = res_trt[intersect(perm_up, rownames(res_trt)), ]
res_xc_repl_up = res_xc[repl_up, ]
res_trt_repl_up = res_trt[intersect(repl_up, rownames(res_trt)), ]
res_xc_stri_up = res_xc[stri_up, ]
res_trt_stri_up = res_trt[intersect(stri_up, rownames(res_trt)), ]

res_xv_perm_dn = res_xv[perm_dn, ]
res_trt_perm_dn = res_trt[intersect(perm_dn, rownames(res_trt)), ]
res_xc_repl_dn = res_xc[repl_dn, ]
res_trt_repl_dn = res_trt[intersect(repl_dn, rownames(res_trt)), ]
res_xc_stri_dn = res_xc[stri_dn, ]
res_trt_stri_dn = res_trt[intersect(stri_dn, rownames(res_trt)), ]

df_up = data.frame(FC = c(mean(exp(res_xv_perm_up$LFC.Genotype1XDP_vanilla)), 
                          mean(exp(res_xc_repl_up$LFC.GenotypeXDP)), 
                          mean(exp(res_xc_stri_up$LFC.GenotypeXDP)), 
                          mean(exp(res_trt_perm_up$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_repl_up$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_stri_up$LFC.GenotypedSVA))),
                   category = c("Permissive XDP Gene Signatures", 
                                "Replicated XDP Gene Signatures", 
                                "Stringent XDP Gene Signatures", 
                                "Permissive XDP Gene Signatures", 
                                "Replicated XDP Gene Signatures", 
                                "Stringent XDP Gene Signatures"),
                   class = c("Permissive XDP Gene Signatures", 
                             "Replicated XDP Gene Signatures", 
                             "Stringent XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures"))

df_up$class = factor(df_up$class, levels = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", "dSVA-rescued XDP Gene Signatures")) 
df_up$direction = "Upregulated"

df_dn = data.frame(FC = c(mean(exp(res_xv_perm_dn$LFC.Genotype1XDP_vanilla)), 
                          mean(exp(res_xc_repl_dn$LFC.GenotypeXDP)), 
                          mean(exp(res_xc_stri_dn$LFC.GenotypeXDP)), 
                          mean(exp(res_trt_perm_dn$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_repl_dn$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_stri_dn$LFC.GenotypedSVA))),
                   category = c("Permissive XDP Gene Signatures", 
                                "Replicated XDP Gene Signatures", 
                                "Stringent XDP Gene Signatures", 
                                "Permissive XDP Gene Signatures", 
                                "Replicated XDP Gene Signatures", 
                                "Stringent XDP Gene Signatures"),
                   class = c("Permissive XDP Gene Signatures", 
                             "Replicated XDP Gene Signatures", 
                             "Stringent XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures"))

df_dn$class = factor(df_dn$class, levels = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", "dSVA-rescued XDP Gene Signatures")) 
df_dn$direction = "Downregulated"
df = rbind(df_up, df_dn)

ggplot(df, aes(x = category, y = 100* (FC -1), fill = class)) + geom_col(width = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  ylim(c(-40, 60)) + 
  xlab("") + ylab("Average Change (%)") +
  geom_noBG() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  facet_wrap(~direction)

ggsave("pdf/barplot_FC_XDPDEG_3Sets_against_dSVA.jpg")


#Average change of dSVA-rescued XDP DEGs vs their Average change against dSVA
perm_rescue_dsva_up = intersect(perm_rescue_dsva, perm_up)
perm_rescue_dsva_dn = intersect(perm_rescue_dsva, perm_dn)
repl_rescue_dsva_up = intersect(repl_rescue_dsva, repl_up)
repl_rescue_dsva_dn = intersect(repl_rescue_dsva, repl_dn)
stri_rescue_dsva_up = intersect(stri_rescue_dsva, stri_up)
stri_rescue_dsva_dn = intersect(stri_rescue_dsva, stri_dn)

res_xv_perm_rescue_dsva_up = res_xv[perm_rescue_dsva_up, ]
res_trt_perm_rescue_dsva_up = res_trt[intersect(perm_rescue_dsva_up, rownames(res_trt)), ]
res_xc_repl_rescue_dsva_up = res_xc[repl_rescue_dsva_up, ]
res_trt_repl_rescue_dsva_up = res_trt[intersect(repl_rescue_dsva_up, rownames(res_trt)), ]
res_xc_stri_rescue_dsva_up = res_xc[stri_rescue_dsva_up, ]
res_trt_stri_rescue_dsva_up = res_trt[intersect(stri_rescue_dsva_up, rownames(res_trt)), ]

res_xv_perm_rescue_dsva_dn = res_xv[perm_rescue_dsva_dn, ]
res_trt_perm_rescue_dsva_dn = res_trt[intersect(perm_rescue_dsva_dn, rownames(res_trt)), ]
res_xc_repl_rescue_dsva_dn = res_xc[repl_rescue_dsva_dn, ]
res_trt_repl_rescue_dsva_dn = res_trt[intersect(repl_rescue_dsva_dn, rownames(res_trt)), ]
res_xc_stri_rescue_dsva_dn = res_xc[stri_rescue_dsva_dn, ]
res_trt_stri_rescue_dsva_dn = res_trt[intersect(stri_rescue_dsva_dn, rownames(res_trt)), ]

df_up = data.frame(FC = c(mean(exp(res_xv_perm_rescue_dsva_up$LFC.Genotype1XDP_vanilla)), 
                          mean(exp(res_xc_repl_rescue_dsva_up$LFC.GenotypeXDP)), 
                          mean(exp(res_xc_stri_rescue_dsva_up$LFC.GenotypeXDP)), 
                          mean(exp(res_trt_perm_rescue_dsva_up$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_repl_rescue_dsva_up$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_stri_rescue_dsva_up$LFC.GenotypedSVA))),
                   category = c("Permissive XDP Gene Signatures Rescued", 
                                "Replicated XDP Gene Signatures Rescued", 
                                "Stringent XDP Gene Signatures Rescued", 
                                "Permissive XDP Gene Signatures Rescued", 
                                "Replicated XDP Gene Signatures Rescued", 
                                "Stringent XDP Gene Signatures Rescued"),
                   class = c("Permissive XDP Gene Signatures", 
                             "Replicated XDP Gene Signatures", 
                             "Stringent XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures"))

df_up$class = factor(df_up$class, levels = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", "dSVA-rescued XDP Gene Signatures")) 
df_up$direction = "Upregulated"

df_dn = data.frame(FC = c(mean(exp(res_xv_perm_rescue_dsva_dn$LFC.Genotype1XDP_vanilla)), 
                          mean(exp(res_xc_repl_rescue_dsva_dn$LFC.GenotypeXDP)), 
                          mean(exp(res_xc_stri_rescue_dsva_dn$LFC.GenotypeXDP)), 
                          mean(exp(res_trt_perm_rescue_dsva_dn$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_repl_rescue_dsva_dn$LFC.GenotypedSVA)), 
                          mean(exp(res_trt_stri_rescue_dsva_dn$LFC.GenotypedSVA))),
                   category = c("Permissive XDP Gene Signatures Rescued", 
                                "Replicated XDP Gene Signatures Rescued", 
                                "Stringent XDP Gene Signatures Rescued", 
                                "Permissive XDP Gene Signatures Rescued", 
                                "Replicated XDP Gene Signatures Rescued", 
                                "Stringent XDP Gene Signatures Rescued"),
                   class = c("Permissive XDP Gene Signatures", 
                             "Replicated XDP Gene Signatures", 
                             "Stringent XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures", 
                             "dSVA-rescued XDP Gene Signatures"))

df_dn$class = factor(df_dn$class, levels = c("Permissive XDP Gene Signatures", "Replicated XDP Gene Signatures", "Stringent XDP Gene Signatures", "dSVA-rescued XDP Gene Signatures")) 
df_dn$direction = "Downregulated"
df = rbind(df_up, df_dn)

write.table(df, paste0("dSVA","_rescued_changes_before_and_after_treatment.txt" ), row.names = F, col.names = T, quote = F, sep = "\t")

ggplot(df, aes(x = category, y = 100* (FC -1), fill = class)) + geom_col(width = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "red", "darkgrey")) +
  ylim(c(-40, 60)) + 
  xlab("") + ylab("Average Change (%)") +
  geom_noBG() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  facet_wrap(~direction)

ggsave("pdf/barplot_FC_XDPDEG_rescue_3Sets_against_dSVA.jpg")






df = data.frame(FC = c(c(exp(res_xv_perm_dn$LFC.Genotype1XDP_vanilla), exp(res_xv_perm_up$LFC.Genotype1XDP_vanilla)),
                       c(exp(res_xc_repl_dn$LFC.GenotypeXDP), exp(res_xc_repl_up$LFC.GenotypeXDP)),
                       c(exp(res_xc_stri_dn$LFC.GenotypeXDP), exp(res_xc_stri_up$LFC.GenotypeXDP)),
                       
                       c(exp(res_xv_perm_rescue_dsva_dn$LFC.Genotype1XDP_vanilla), exp(res_xv_perm_rescue_dsva_up$LFC.Genotype1XDP_vanilla)),
                       c(exp(res_xc_repl_rescue_dsva_dn$LFC.GenotypeXDP), exp(res_xc_repl_rescue_dsva_up$LFC.GenotypeXDP)),
                       c(exp(res_xc_stri_rescue_dsva_dn$LFC.GenotypeXDP), exp(res_xc_stri_rescue_dsva_up$LFC.GenotypeXDP)),
                       
                       c(exp(res_trt_perm_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_perm_rescue_dsva_up$LFC.GenotypedSVA)),
                       c(exp(res_trt_repl_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_repl_rescue_dsva_up$LFC.GenotypedSVA)),
                       c(exp(res_trt_stri_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_stri_rescue_dsva_up$LFC.GenotypedSVA))),
                
                category = c(rep("Permissive XDP Gene Signatures", length(c(exp(res_xv_perm_dn$LFC.Genotype1XDP_vanilla), exp(res_xv_perm_up$LFC.Genotype1XDP_vanilla)))),
                             rep("Replicated XDP Gene Signatures", length(c(exp(res_xc_repl_dn$LFC.GenotypeXDP), exp(res_xc_repl_up$LFC.GenotypeXDP)))),
                             rep("Stringent XDP Gene Signatures", length(c(exp(res_xc_stri_dn$LFC.GenotypeXDP), exp(res_xc_stri_up$LFC.GenotypeXDP)))),
                             
                             rep("Permissive XDP Gene Signatures", length(c(exp(res_xv_perm_rescue_dsva_dn$LFC.Genotype1XDP_vanilla), exp(res_xv_perm_rescue_dsva_up$LFC.Genotype1XDP_vanilla)))), 
                             rep("Replicated XDP Gene Signatures", length(c(exp(res_xc_repl_rescue_dsva_dn$LFC.GenotypeXDP), exp(res_xc_repl_rescue_dsva_up$LFC.GenotypeXDP)))), 
                             rep("Stringent XDP Gene Signatures", length(c(exp(res_xc_stri_rescue_dsva_dn$LFC.GenotypeXDP), exp(res_xc_stri_rescue_dsva_up$LFC.GenotypeXDP)))), 
                             
                             rep("Permissive XDP Gene Signatures", length(c(exp(res_trt_perm_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_perm_rescue_dsva_up$LFC.GenotypedSVA)))), 
                             rep("Replicated XDP Gene Signatures", length(c(exp(res_trt_repl_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_repl_rescue_dsva_up$LFC.GenotypedSVA)))), 
                             rep("Stringent XDP Gene Signatures", length(c(exp(res_trt_stri_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_stri_rescue_dsva_up$LFC.GenotypedSVA))))),
                
                class = c(rep("Permissive all", length(c(exp(res_xv_perm_dn$LFC.Genotype1XDP_vanilla), exp(res_xv_perm_up$LFC.Genotype1XDP_vanilla)))),
                          rep("Replicated all", length(c(exp(res_xc_repl_dn$LFC.GenotypeXDP), exp(res_xc_repl_up$LFC.GenotypeXDP)))),
                          rep("Stringent all", length(c(exp(res_xc_stri_dn$LFC.GenotypeXDP), exp(res_xc_stri_up$LFC.GenotypeXDP)))),
                          
                          rep("dSVA-rescued", length(c(exp(res_xv_perm_rescue_dsva_dn$LFC.Genotype1XDP_vanilla), exp(res_xv_perm_rescue_dsva_up$LFC.Genotype1XDP_vanilla)))), 
                          rep("dSVA-rescued", length(c(exp(res_xc_repl_rescue_dsva_dn$LFC.GenotypeXDP), exp(res_xc_repl_rescue_dsva_up$LFC.GenotypeXDP)))), 
                          rep("dSVA-rescued", length(c(exp(res_xc_stri_rescue_dsva_dn$LFC.GenotypeXDP), exp(res_xc_stri_rescue_dsva_up$LFC.GenotypeXDP)))), 
                          
                          rep("dSVA-rescued in treated vs untreated XDP", length(c(exp(res_trt_perm_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_perm_rescue_dsva_up$LFC.GenotypedSVA)))), 
                          rep("dSVA-rescued in treated vs untreated XDP", length(c(exp(res_trt_repl_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_repl_rescue_dsva_up$LFC.GenotypedSVA)))), 
                          rep("dSVA-rescued in treated vs untreated XDP", length(c(exp(res_trt_stri_rescue_dsva_dn$LFC.GenotypedSVA), exp(res_trt_stri_rescue_dsva_up$LFC.GenotypedSVA))))))

df$class = factor(df$class, levels = c("Permissive all", "Replicated all", "Stringent all", "dSVA-rescued", "dSVA-rescued in treated vs untreated XDP")) 
df = df[which(!grepl("Stringent", df$category)), ]

p = ggplot(df, aes(x = class, y = 100* (FC - 0))) + 
      geom_sina(aes(color = class), size = 0.1, position = position_dodge(1), maxwidth = 0.5, scale = "count") + 
      geom_hline(yintercept = 100, linetype = "dashed") +
      scale_color_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "darkgrey", "white")) +
      xlab("") + 
      ylab("Fold Change (%)") +
      geom_noBG() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
      scale_y_continuous(limits = c(0, 200), breaks=seq(0, 200, 20)) +
      theme(axis.text.x = element_text(angle = 45, hjust=1), strip.text = element_text(size = 8), legend.position="none") +
      facet_wrap(~category, scales = "free_x")

p
ggsave("pdf/sinaplot_FC_XDPDEG_rescue_2Sets_against_dSVA_white.jpg", width = 6.1, height = 4.24)

p + scale_color_manual(name = "class", labels = function(x) str_wrap(x, width = 15), values = c("#0066FF", "#EC8014", "darkgrey", "pink"))
ggsave("pdf/sinaplot_FC_XDPDEG_rescue_2Sets_against_dSVA.jpg", width = 6.1, height = 4.24)


####################################################################################################################################################################################
#This study
#dSVA vs CON Pathway analyses
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

deg_dsva_gs = ref[deg_dsva, 1]
deg_dsva_up_gs = ref[deg_dsva_up, 1]
deg_dsva_dn_gs = ref[deg_dsva_dn, 1]
bg_dsva_gs = ref[bg_dsva, 1]

deg_dsva_p_gs = ref[deg_dsva_p, 1]
deg_dsva_up_p_gs = ref[deg_dsva_up_p, 1]
deg_dsva_dn_p_gs = ref[deg_dsva_dn_p, 1]

#GO for DEGs at FDR level
deg_dsva_go = go_test(deg_dsva_gs, bg_dsva_gs, godb_list)
deg_dsva_up_go = go_test(deg_dsva_up_gs, bg_dsva_gs, godb_list)
deg_dsva_dn_go = go_test(deg_dsva_dn_gs, bg_dsva_gs, godb_list)
write.table(deg_dsva_go, file = paste0("GSEA_GO_dSVADEG_","FDR_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(deg_dsva_up_go, file = paste0("GSEA_GO_dSVADEG_","FDR_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(deg_dsva_dn_go, file = paste0("GSEA_GO_dSVADEG_","FDR_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#GO for DEGs at Nominal level
deg_dsva_p_go = go_test(deg_dsva_p_gs, bg_dsva_gs, godb_list)
deg_dsva_up_p_go = go_test(deg_dsva_up_p_gs, bg_dsva_gs, godb_list)
deg_dsva_dn_p_go = go_test(deg_dsva_dn_p_gs, bg_dsva_gs, godb_list)
write.table(deg_dsva_p_go, file = paste0("GSEA_GO_dSVADEG_","pval_all",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(deg_dsva_up_p_go, file = paste0("GSEA_GO_dSVADEG_","pval_up",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")
write.table(deg_dsva_dn_p_go, file = paste0("GSEA_GO_dSVADEG_","pval_dn",".txt"), row.names = T, col.names = T, quote = F, sep = "\t")

#Read in GO results
deg_dsva_p_go = read.table(paste0("GSEA_GO_dSVADEG_","pval_all",".txt"), header = T, row.names = 1, sep = "\t")
deg_dsva_up_p_go = read.table(paste0("GSEA_GO_dSVADEG_","pval_up",".txt"), header = T, row.names = 1, sep = "\t")
deg_dsva_dn_p_go = read.table(paste0("GSEA_GO_dSVADEG_","pval_dn",".txt"), header = T, row.names = 1, sep = "\t")
deg_dsva_go = read.table(paste0("GSEA_GO_dSVADEG_","FDR_all",".txt"), header = T, row.names = 1, sep = "\t")
deg_dsva_up_go = read.table(paste0("GSEA_GO_dSVADEG_","FDR_up",".txt"), header = T, row.names = 1, sep = "\t")
deg_dsva_dn_go = read.table(paste0("GSEA_GO_dSVADEG_","FDR_dn",".txt"), header = T, row.names = 1, sep = "\t")

#force the same order
deg_dsva_p_go = deg_dsva_p_go[rownames(deg_dsva_p_go), ]
deg_dsva_up_p_go = deg_dsva_up_p_go[rownames(deg_dsva_p_go), ]
deg_dsva_dn_p_go = deg_dsva_dn_p_go[rownames(deg_dsva_p_go), ]
deg_dsva_go = deg_dsva_go[rownames(deg_dsva_p_go), ]
deg_dsva_up_go = deg_dsva_up_go[rownames(deg_dsva_p_go), ]
deg_dsva_dn_go = deg_dsva_dn_go[rownames(deg_dsva_p_go), ]


#extract terms of interest
aloy_vec = scan("go_terms_roi_v2.txt", what = "character")
keyword_vec = c(aloy_vec, c("NEURO", "SYNAP"))
deg_dsva_p_go_sel = go_extract_keywords(deg_dsva_p_go[which(deg_dsva_p_go$pvalue < 0.05), ], "term", keyword_vec)
deg_dsva_go_sel = go_extract_keywords(deg_dsva_go[which(deg_dsva_go$pvalue < 0.05), ], "term", keyword_vec)
deg_dsva_dn_p_go_sel = go_extract_keywords(deg_dsva_dn_p_go[which(deg_dsva_dn_p_go$pvalue < 0.05), ], "term", keyword_vec)
deg_dsva_dn_go_sel = go_extract_keywords(deg_dsva_dn_go[which(deg_dsva_dn_go$pvalue < 0.05), ], "term", keyword_vec)
deg_dsva_up_p_go_sel = go_extract_keywords(deg_dsva_up_p_go[which(deg_dsva_up_p_go$pvalue < 0.05), ], "term", keyword_vec)
deg_dsva_up_go_sel = go_extract_keywords(deg_dsva_up_go[which(deg_dsva_up_go$pvalue < 0.05), ], "term", keyword_vec)


#non-directional DEGs for GO
go_plot(deg_dsva_p_go_sel[order(deg_dsva_p_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_dSVADEG_","pval_all",".jpg"))
go_plot(deg_dsva_go_sel[order(deg_dsva_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_dSVADEG_","FDR_all",".jpg"))
intersect(deg_dsva_p_go_sel[order(deg_dsva_p_go_sel$pvalue), ][1:20, "term"], deg_dsva_go_sel[order(deg_dsva_go_sel$pvalue), ][1:20, "term"])


#Downregulated DEGs for GO
go_plot(deg_dsva_dn_p_go_sel[order(deg_dsva_dn_p_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_dSVADEG_","pval_dn",".jpg"))
go_plot(deg_dsva_dn_go_sel[order(deg_dsva_dn_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_dSVADEG_","FDR_dn",".jpg"))
intersect(deg_dsva_dn_p_go_sel[order(deg_dsva_dn_p_go_sel$pvalue), ][, "term"], deg_dsva_dn_go_sel[order(deg_dsva_dn_go_sel$pvalue), ][, "term"])


#Upregulated DEGs for GO
go_plot(deg_dsva_up_p_go_sel[order(deg_dsva_up_p_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_dSVADEG_","pval_up",".jpg"))
go_plot(deg_dsva_up_go_sel[order(deg_dsva_up_go_sel$pvalue), ][1:20, ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_dSVADEG_","FDR_up",".jpg"))
intersect(deg_dsva_up_p_go_sel[order(deg_dsva_up_p_go_sel$pvalue), ][1:20, "term"], deg_dsva_up_go_sel[order(deg_dsva_up_go_sel$pvalue), ][1:20, "term"])


#Read in GO results for XDP DEGs
perm_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_all",".txt"), header = T, row.names = 1, sep = "\t")
perm_up_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_up",".txt"), header = T, row.names = 1, sep = "\t")
perm_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","permissive_dn",".txt"), header = T, row.names = 1, sep = "\t")
repl_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_all",".txt"), header = T, row.names = 1, sep = "\t")
repl_up_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_up",".txt"), header = T, row.names = 1, sep = "\t")
repl_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","replication_dn",".txt"), header = T, row.names = 1, sep = "\t")
stri_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_all",".txt"), header = T, row.names = 1, sep = "\t")
stri_up_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_up",".txt"), header = T, row.names = 1, sep = "\t")
stri_dn_go = read.table(paste0("GSEA_GO_XDPDEG_","stringent_dn",".txt"), header = T, row.names = 1, sep = "\t")

#extract terms of interest
aloy_vec = scan("go_terms_roi.txt", what = "character")
keyword_vec = c(aloy_vec, c("NEURO", "SYNAP"))
perm_go_sel = go_extract_keywords(perm_go[which(perm_go$pvalue < 0.05), ], "term", keyword_vec)
perm_dn_go_sel = go_extract_keywords(perm_dn_go[which(perm_dn_go$pvalue < 0.05), ], "term", keyword_vec)
perm_up_go_sel = go_extract_keywords(perm_up_go[which(perm_up_go$pvalue < 0.05), ], "term", keyword_vec)
repl_go_sel = go_extract_keywords(repl_go[which(repl_go$pvalue < 0.05), ], "term", keyword_vec)
repl_dn_go_sel = go_extract_keywords(repl_dn_go[which(repl_dn_go$pvalue < 0.05), ], "term", keyword_vec)
repl_up_go_sel = go_extract_keywords(repl_up_go[which(repl_up_go$pvalue < 0.05), ], "term", keyword_vec)
stri_go_sel = go_extract_keywords(stri_go[which(stri_go$pvalue < 0.05), ], "term", keyword_vec)
stri_dn_go_sel = go_extract_keywords(stri_dn_go[which(stri_dn_go$pvalue < 0.05), ], "term", keyword_vec)
stri_up_go_sel = go_extract_keywords(stri_up_go[which(stri_up_go$pvalue < 0.05), ], "term", keyword_vec)


#Comparison of GO Terms from DEGs before and after dSVA
p = go_change_plot(perm_go_sel, deg_dsva_p_go, "Permissive XDP Gene Signatures", "Nominal dSVA Gene Signatures", "query_annotated", "pvalue", 0.05, ntop = 20, text_width = 30, text_size = 6, transform_val = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c("Nominal dSVA Gene Signatures", "Permissive XDP Gene Signatures"), 15), values = c("darkgrey", "#0066FF")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Permissive_Set_vs_dSVAvCON_Nominal_Set.jpg"), width = 6.1, height = 4.24)

p = go_change_plot(repl_go_sel, deg_dsva_p_go[, 1:15], "Replicated XDP Gene Signatures", "dSVA Nominal Gene Signatures", "query_annotated", "pvalue", 0.05, ntop = 20, text_width = 30, text_size = 6, transform_val = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c("Nominal dSVA Gene Signatures", "Replicated XDP Gene Signatures"), 15), values = c("darkgrey", "#EC8014")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Replicated_Set_vs_dSVAvCON_Nominal_Set.jpg"), width = 6.1, height = 4.24)


p = go_change_plot(stri_go_sel, deg_dsva_go[, 1:15], "Stringent XDP Gene Signatures", "FDR dSVA Gene Signatures", "query_annotated", "pvalue", 0.05, ntop = 20, text_width = 30, text_size = 6, transform_val = F)
p + scale_fill_manual(name = "Set", labels = str_wrap(c("FDR dSVA Gene Signatures", "Stringent XDP Gene Signatures"), 15), values = c("darkgrey", "red")) + ylab("Number of Hit Gene")
ggsave(paste0("pdf/GOplot_changes_Stringent_Set_vs_dSVAvCON_FDR_Set.jpg"), width = 6.1, height = 4.24)

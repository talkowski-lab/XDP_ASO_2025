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


#######################################################################################
#Define grid-search functions 
#######################################################################################
godb_list = list(
  GO_BP = "~/ref/gsea/gsea_go_bp.rds",
  GO_CC = "~/ref/gsea/gsea_go_cc.rds",
  GO_MF = "~/ref/gsea/gsea_go_mf.rds",
  PATHWAYS = "~/ref/gsea/gsea_pathways.rds",
  REG = "~/ref/gsea/gsea_reg.rds",
  HPO = "~/ref/gsea/gsea_hpo.rds"
)

res_xc = read.table("res_FDR_LFC_glm_Exp_Genotype_XDP_vs_CON.txt", header = T, sep = "\t", check.names = F)
rownames(res_xc) = as.vector(res_xc$name)

p = search_fc_cutoff_go(res_xc, col_fc = "LFC.GenotypeXDP", fc_cutoffs = seq(1.1, 2.7, 0.1), col_p = "pval.GenotypeXDP", p_cutoff = 0.05, go_list = godb_list)
ggsave("pdf/barplot_FC_cutoff_vs_minimal_FDR_GO.jpg")

deg = rownames(res_xc[which((abs(res_xc[["LFC.GenotypeXDP"]]) >= log(1.2)) & (res_xc[["pval.GenotypeXDP"]] < 0.05)), ])
gs_ref = "GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt"
ref = read.table(gs_ref, header = F, row.names = 1, sep = "\t")
deg_gs = ref[deg, 1]
bg_gs = ref[rownames(res_xc), 1]
deg_go = go_test(deg_gs, bg_gs, godb_list)
head(deg_go[, 1:6], 3)
go_plot(deg_go[1:3, ], "fdr", 0.1, text_width = 30, text_size = 6) + geom_fill(c("#0066FF"))
ggsave("pdf/GOplot_permissive_set_at_FC120percent.jpg")

matrisome_genes = unique(c(strsplit(deg_go[1,ncol(deg_go)], ", ")[[1]], strsplit(deg_go[2,ncol(deg_go)], ", ")[[1]], strsplit(deg_go[3,ncol(deg_go)], ", ")[[1]]))
summary(res_xc[deg1, "FDR.GenotypeXDP"])

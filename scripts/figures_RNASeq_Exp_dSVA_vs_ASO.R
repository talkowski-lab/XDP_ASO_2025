library(ggplot2)
library(VennDiagram)
library(stringr)
source("/home/dadi/projects/xdp_aso/nsc/scripts/go.R")

setwd("~/projects/xdp_aso/nsc/output/")



godb_list = list(
  GO_BP = "~/ref/gsea/gsea_go_bp.rds",
  GO_CC = "~/ref/gsea/gsea_go_cc.rds",
  GO_MF = "~/ref/gsea/gsea_go_mf.rds",
  PATHWAYS = "~/ref/gsea/gsea_pathways.rds",
  REG = "~/ref/gsea/gsea_reg.rds",
  HPO = "~/ref/gsea/gsea_hpo.rds"
)
ref = read.table("GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", header = F, row.names = 1, sep = "\t")

aloy_vec = scan("go_terms_roi_v2.txt", what = "character")
keyword_vec = c(aloy_vec, c("NEURO", "SYNAP"))

res_dsva = read.table("res_FDR_LFC_glm_Exp_Trt_dSVA_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_880 = read.table("res_FDR_LFC_glm_Exp_Trt_XDP880_vs_CON.txt", header = T, sep = "\t", check.names = F)
res_131 = read.table("res_FDR_LFC_glm_Exp_Trt_XDP131_vs_CON.txt", header = T, sep = "\t", check.names = F)
bg = Reduce(intersect, list(as.vector(res_dsva$name), as.vector(res_880$name), as.vector(res_131$name)))
bg_gs = ref[bg, 1]


perm_rescue_dsva = readRDS("RData/perm_rescue_dsva.rds")
perm_rescue_131 = readRDS("RData/perm_rescue_ASO131.rds")
perm_rescue_880 = readRDS("RData/perm_rescue_ASO880.rds")

v2 <- venn.diagram(list(dSVA = perm_rescue_dsva, ASO131 = perm_rescue_131, ASO880 = perm_rescue_880),
                   fill = c("white", "pink", "purple"),
                   alpha = c(0.5, 0.5, 0.5),
                   print.mode="raw",
                   filename=NULL)

jpeg("pdf/venn_permissive_set_rescued.jpg")
grid.newpage()
grid.draw(v2)
dev.off()

perm_rescue_ov = Reduce(intersect, list(perm_rescue_dsva, perm_rescue_880, perm_rescue_131))
perm_rescue_ov_gs = ref[perm_rescue_ov, 1]
perm_rescue_ov_go = go_test(perm_rescue_ov_gs, bg_gs, godb_list)
perm_rescue_ov_go_sel = go_extract_keywords(perm_rescue_ov_go[which(perm_rescue_ov_go$pvalue < 0.05), ], "term", keyword_vec)

go_plot(perm_rescue_ov_go_sel[order(perm_rescue_ov_go_sel$pvalue), ][1:10,], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_ASOdSVA_sharedRescue_permissive_set_pval_all.jpg"))


repl_rescue_dsva = readRDS("RData/repl_rescue_dsva.rds")
repl_rescue_131 = readRDS("RData/repl_rescue_ASO131.rds")
repl_rescue_880 = readRDS("RData/repl_rescue_ASO880.rds")

v2 <- venn.diagram(list(dSVA = repl_rescue_dsva, ASO131 = repl_rescue_131, ASO880 = repl_rescue_880),
                   fill = c("white", "pink", "purple"),
                   alpha = c(0.5, 0.5, 0.5),
                   print.mode="raw",
                   filename=NULL)

jpeg("pdf/venn_replication_set_rescued.jpg")
grid.newpage()
grid.draw(v2)
dev.off()

repl_rescue_ov = Reduce(intersect, list(repl_rescue_dsva, repl_rescue_880, repl_rescue_131))
repl_rescue_ov_gs = ref[repl_rescue_ov, 1]
repl_rescue_ov_go = go_test(repl_rescue_ov_gs, bg_gs, godb_list)
repl_rescue_ov_go_sel = go_extract_keywords(repl_rescue_ov_go[which(repl_rescue_ov_go$pvalue < 0.05), ], "term", keyword_vec)

go_plot(repl_rescue_ov_go_sel[order(repl_rescue_ov_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6, show_hit = F)
ggsave(paste0("pdf/GOplot_ASOdSVA_sharedRescue_replication_set_pval_all.jpg"))


stri_rescue_dsva = readRDS("RData/stri_rescue_dsva.rds")
stri_rescue_131 = readRDS("RData/stri_rescue_aso131.rds")
stri_rescue_880 = readRDS("RData/stri_rescue_aso880.rds")

v2 <- venn.diagram(list(dSVA = stri_rescue_dsva, ASO131 = stri_rescue_131, ASO880 = stri_rescue_880),
                   fill = c("white", "pink", "purple"),
                   alpha = c(0.5, 0.5, 0.5),
                   print.mode="raw",
                   filename=NULL)

jpeg("pdf/venn_strigent_set_rescued.jpg")
grid.newpage()
grid.draw(v2)
dev.off()

stri_rescue_ov = Reduce(intersect, list(stri_rescue_dsva, stri_rescue_880, stri_rescue_131))
stri_rescue_ov_gs = ref[stri_rescue_ov, 1]
stri_rescue_ov_go = go_test(stri_rescue_ov_gs, bg_gs, godb_list)
stri_rescue_ov_go_sel = go_extract_keywords(stri_rescue_ov_go[which(stri_rescue_ov_go$pvalue < 0.05), ], "term", keyword_vec)

go_plot(stri_rescue_ov_go_sel[order(stri_rescue_ov_go_sel$pvalue), ], "pvalue", 0.05, text_width = 30, text_size = 6)
ggsave(paste0("pdf/GOplot_ASOdSVA_sharedRescue_stringent_set_pval_all.jpg"))

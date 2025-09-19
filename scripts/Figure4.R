## Transcriptome-side changes in CRISPR/ASO tratment ##

setwd("/Volumes//talkowski/Samples/XDP/ASO_dSVA/Results_draft1/2025-06-18_enrichment.results")

save.image(file='../RY_Fig4.session.RData')

DEGs_XDP_DF<-read_tsv("NO.CON+NO.XDP-LFC+Var.tsv")

DEGs_XDP_FDR<-DEGs_XDP_DF[DEGs_XDP_DF$fdr < 0.1,]

DEGs_XDP_Nom<-DEGs_XDP_DF[DEGs_XDP_DF$pvalue < 0.05 & abs(DEGs_XDP_DF$lfc) >= log(1.2),]

DEGs_XDP_FDR_Enrichment<-read_tsv("Enrichments-NO.CON.vs.NO.XDP-FDR.tsv")
DEGs_XDP_Nom_Enrichment<-read_tsv("Enrichments-NO.CON.vs.NO.XDP-Nominal.tsv")
 
RescuedGenes_DF<-read_tsv("RescueData-AllTreatments-AllCriteria.tsv")

#Number of treatments
unique(RescuedGenes_DF$Treatment)
# "dSVA" "131"  "307"  "876"  "877"  "879"  "880"  "881" 

RescuedGenes_FDR<-RescuedGenes_DF[RescuedGenes_DF$EnsemblID %in% DEGs_XDP_FDR$EnsemblID,]
RescuedGenes_Nom<-RescuedGenes_DF[RescuedGenes_DF$EnsemblID %in% DEGs_XDP_Nom$EnsemblID,]

RescuedGenes_DF %>% group_by(Treatment) %>% filter(isRescued.2ZTests == TRUE) %>% summarise(count = n())
# 1 131         279 (Fam: 131)
# 2 307         387 (Fam: 126)
# 3 876         348 (Fam: 131)
# 4 877         277 (Fam: 131)
# 5 879         257 (Fam: 134)
# 6 880         326 (Fam: 134)
# 7 881         203 (Fam: 134)
# 8 dSVA        139

RescuedGenes_FDR %>% group_by(Treatment) %>% filter(isRescued.Naive == TRUE) %>% summarise(count = n())
# 1 131         446
# 2 307         476
# 3 876         513
# 4 877         429
# 5 879         413
# 6 880         472
# 7 881         424
# 8 dSVA        456

RescuedGenes_Nom %>% group_by(Treatment) %>% filter(isRescued.2ZTests == TRUE) %>% summarise(count = n())
# 1 131         226
# 2 307         352
# 3 876         304
# 4 877         213
# 5 879         252
# 6 880         273
# 7 881         183
# 8 dSVA        155

RescuedGenes_Nom %>% group_by(Treatment) %>% filter(isRescued.Naive == TRUE) %>% summarise(count = n())
# 1 131         690
# 2 307         650
# 3 876         700
# 4 877         659
# 5 879         652
# 6 880         702
# 7 881         670
# 8 dSVA        705
## 877 doesn't rescue AS

FDR.enrichment.results.df <- load_enrichment_data(file_pattern='Enrichments-RescuedGenes-.*-isRescued.2ZTests-XDP.FDR.tsv')
Nom.enrichment.results.df <- load_enrichment_data(file_pattern='Enrichments-RescuedGenes-.*-isRescued.Both-XDP.Nominal.tsv')


FDR.enrichment.results.df_naive <- load_enrichment_data(file_pattern='Enrichments-RescuedGenes-.*isRescued.Naive-XDP.FDR.tsv')


rescuedGenes_list_FDR<-list(
 dSVA=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="dSVA" & RescuedGenes_FDR$isRescued.Naive == TRUE],
  ASO131=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="131" & RescuedGenes_FDR$isRescued.Naive == TRUE],
  ASO307=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="307" & RescuedGenes_FDR$isRescued.Naive == TRUE],
  ASO876=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="876" & RescuedGenes_FDR$isRescued.Naive == TRUE],
#  ASO877=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="877" & RescuedGenes_FDR$isRescued.Naive == TRUE],
  ASO879=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="879" & RescuedGenes_FDR$isRescued.Naive == TRUE],
  ASO880=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="880" & RescuedGenes_FDR$isRescued.Naive == TRUE],
  ASO881=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="881" & RescuedGenes_FDR$isRescued.Naive == TRUE]
)

rescuedGenes_list_Nom<-list(
# dSVA=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="dSVA" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO131=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="131" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO307=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="307" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO876=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="876" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO877=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="877" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO879=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="879" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO880=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="880" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO881=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="881" & RescuedGenes_Nom$isRescued.Naive == TRUE]
)


rescuedGenes_list_Nom<-list(
  # dSVA=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="dSVA" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO131=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="131" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO307=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="307" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO876=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="876" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO877=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="877" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO879=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="879" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO880=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="880" & RescuedGenes_Nom$isRescued.Naive == TRUE],
  ASO881=RescuedGenes_Nom$gene_name[RescuedGenes_Nom$Treatment=="881" & RescuedGenes_Nom$isRescued.Naive == TRUE]
)

XDPSignaturePathways<-DEGs_XDP_FDR_Enrichment[DEGs_XDP_FDR_Enrichment$p.value < 0.05,]

XDPSignaturePathways$RescuedBydSVA<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,productdSVA)
XDPSignaturePathways$RescuedBy880<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product880)
XDPSignaturePathways$RescuedBy131<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product131)
XDPSignaturePathways$RescuedBy307<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product307)
XDPSignaturePathways$RescuedBy879<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product879)
XDPSignaturePathways$RescuedBy876<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product876)
XDPSignaturePathways$RescuedBy881<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product881)
XDPSignaturePathways$RescuedBy877<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,product877)
XDPSignaturePathways$RescuedByAnyTreatment<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,productAnyTreatment)
XDPSignaturePathways$RescuedByAllTreatment<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,productAllTreatment)
XDPSignaturePathways$RescuedByAnyASO<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,productAnyASO)
XDPSignaturePathways$RescuedByAllASO<-apply(as.data.frame(XDPSignaturePathways$is_term_is_hit),1,productAllASO)

XDPSignaturePathways$RescuedBydSVA_numbers<-listLen(XDPSignaturePathways$RescuedBydSVA)
XDPSignaturePathways$RescuedBydSVA_perc<-(XDPSignaturePathways$RescuedBydSVA_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy880_numbers<-listLen(XDPSignaturePathways$RescuedBy880)
XDPSignaturePathways$RescuedBy880_perc<-(XDPSignaturePathways$RescuedBy880_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy131_numbers<-listLen(XDPSignaturePathways$RescuedBy131)
XDPSignaturePathways$RescuedBy131_perc<-(XDPSignaturePathways$RescuedBy131_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy307_numbers<-listLen(XDPSignaturePathways$RescuedBy307)
XDPSignaturePathways$RescuedBy307_perc<-(XDPSignaturePathways$RescuedBy307_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy879_numbers<-listLen(XDPSignaturePathways$RescuedBy879)
XDPSignaturePathways$RescuedBy879_perc<-(XDPSignaturePathways$RescuedBy879_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy876_numbers<-listLen(XDPSignaturePathways$RescuedBy876)
XDPSignaturePathways$RescuedBy876_perc<-(XDPSignaturePathways$RescuedBy876_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy881_numbers<-listLen(XDPSignaturePathways$RescuedBy881)
XDPSignaturePathways$RescuedBy881_perc<-(XDPSignaturePathways$RescuedBy881_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedBy877_numbers<-listLen(XDPSignaturePathways$RescuedBy877)
XDPSignaturePathways$RescuedBy877_perc<-(XDPSignaturePathways$RescuedBy877_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedByAnyTreatment_numbers<-listLen(XDPSignaturePathways$RescuedByAnyTreatment)
XDPSignaturePathways$RescuedByAnyTreatment_perc<-(XDPSignaturePathways$RescuedByAnyTreatment_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedByAllTreatment_numbers<-listLen(XDPSignaturePathways$RescuedByAllTreatment)
XDPSignaturePathways$RescuedByAllTreatment_perc<-(XDPSignaturePathways$RescuedByAllTreatment_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedByAnyASO_numbers<-listLen(XDPSignaturePathways$RescuedByAnyASO)
XDPSignaturePathways$RescuedByAnyASO_perc<-(XDPSignaturePathways$RescuedByAnyASO_numbers/XDPSignaturePathways$is_term_is_hit)*100
XDPSignaturePathways$RescuedByAllASO_numbers<-listLen(XDPSignaturePathways$RescuedByAllASO)
XDPSignaturePathways$RescuedByAllASO_perc<-(XDPSignaturePathways$RescuedByAllASO_numbers/XDPSignaturePathways$is_term_is_hit)*100

## Nominal DEG pathway rescue ##
XDPSignaturePathways_Nom<-DEGs_XDP_Nom_Enrichment[DEGs_XDP_Nom_Enrichment$p.value < 0.05,]

#XDPSignaturePathways_Nom$RescuedBydSVA<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,productdSVA)
XDPSignaturePathways_Nom$RescuedBy880<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product880)
XDPSignaturePathways_Nom$RescuedBy131<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product131)
#XDPSignaturePathways_Nom$RescuedBy307<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product307)
XDPSignaturePathways_Nom$RescuedBy879<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product879)
XDPSignaturePathways_Nom$RescuedBy876<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product876)
XDPSignaturePathways_Nom$RescuedBy881<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product881)
#XDPSignaturePathways_Nom$RescuedBy877<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,product877)
XDPSignaturePathways_Nom$RescuedByAnyTreatment<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,productAnyTreatment)
#XDPSignaturePathways_Nom$RescuedByAllTreatment<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,productAllTreatment)
XDPSignaturePathways_Nom$RescuedByAnyASO<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,productAnyASO)
#XDPSignaturePathways_Nom$RescuedByAllASO<-apply(as.data.frame(XDPSignaturePathways_Nom$is_term_is_hit),1,productAllASO)

# XDPSignaturePathways_Nom$RescuedBydSVA_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBydSVA)
# XDPSignaturePathways_Nom$RescuedBydSVA_perc<-(XDPSignaturePathways_Nom$RescuedBydSVA_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedBy880_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy880)
XDPSignaturePathways_Nom$RescuedBy880_perc<-(XDPSignaturePathways_Nom$RescuedBy880_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedBy131_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy131)
XDPSignaturePathways_Nom$RescuedBy131_perc<-(XDPSignaturePathways_Nom$RescuedBy131_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
# XDPSignaturePathways_Nom$RescuedBy307_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy307)
# XDPSignaturePathways_Nom$RescuedBy307_perc<-(XDPSignaturePathways_Nom$RescuedBy307_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedBy879_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy879)
XDPSignaturePathways_Nom$RescuedBy879_perc<-(XDPSignaturePathways_Nom$RescuedBy879_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedBy876_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy876)
XDPSignaturePathways_Nom$RescuedBy876_perc<-(XDPSignaturePathways_Nom$RescuedBy876_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedBy881_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy881)
XDPSignaturePathways_Nom$RescuedBy881_perc<-(XDPSignaturePathways_Nom$RescuedBy881_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
# XDPSignaturePathways_Nom$RescuedBy877_numbers<-listLen(XDPSignaturePathways_Nom$RescuedBy877)
# XDPSignaturePathways_Nom$RescuedBy877_perc<-(XDPSignaturePathways_Nom$RescuedBy877_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedByAnyTreatment_numbers<-listLen(XDPSignaturePathways_Nom$RescuedByAnyTreatment)
XDPSignaturePathways_Nom$RescuedByAnyTreatment_perc<-(XDPSignaturePathways_Nom$RescuedByAnyTreatment_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
# XDPSignaturePathways_Nom$RescuedByAllTreatment_numbers<-listLen(XDPSignaturePathways_Nom$RescuedByAllTreatment)
# XDPSignaturePathways_Nom$RescuedByAllTreatment_perc<-(XDPSignaturePathways_Nom$RescuedByAllTreatment_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
XDPSignaturePathways_Nom$RescuedByAnyASO_numbers<-listLen(XDPSignaturePathways_Nom$RescuedByAnyASO)
XDPSignaturePathways_Nom$RescuedByAnyASO_perc<-(XDPSignaturePathways_Nom$RescuedByAnyASO_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100
# XDPSignaturePathways_Nom$RescuedByAllASO_numbers<-listLen(XDPSignaturePathways_Nom$RescuedByAllASO)
# XDPSignaturePathways_Nom$RescuedByAllASO_perc<-(XDPSignaturePathways_Nom$RescuedByAllASO_numbers/XDPSignaturePathways_Nom$is_term_is_hit)*100

### 1490 rescue by 2Ztest ##
# 1 131         279 (Fam: 131)
# 2 307         387 (Fam: 126)
# 3 876         348 (Fam: 131)
# 4 877         277 (Fam: 131)
# 5 879         257 (Fam: 134)
# 6 880         326 (Fam: 134)
# 7 881         203 (Fam: 134)

rescuedGenes_list_FDR_Ztest<-list(
  DEGs_XDP_FDR=DEGs_XDP_FDR$gene_name
 dSVA=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="dSVA" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 ASO307=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="307" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 ASO131=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="131" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 ASO876=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="876" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 # ASO877=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="877" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 ASO880=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="880" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 ASO881=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="881" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
 ASO879=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="879" & RescuedGenes_FDR$isRescued.2ZTests == TRUE]
 )
AnyTreatmentRescue<-Reduce(union,rescuedGenes_list_FDR_Ztest)
AnyTreatmentRescue_3Fam
dSVA_rescued<-RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="dSVA" & RescuedGenes_FDR$isRescued.2ZTests == TRUE]
pdf("/Volumes/Macintosh HD/Users/ry077/Dropbox (Partners HealthCare)/XDP/ASOpaper2023/V1Figures/Temp_ASO2Ztets_rescue_overlaps_heatmap.pdf",height = 4, width = 4)
Univ=as.vector(unlist(Reduce(union,list(rescuedGenes_list_FDR_Ztest))))
plotEnrichments(rescuedGenes_list_FDR_Ztest,rescuedGenes_list_FDR_Ztest,Univ)
dev.off()

onlydSVA_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="dSVA" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)
only307_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="307" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)
only131_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="131" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)
only876_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="876" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)
only880_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="880" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)
only881_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="881" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)
only879_rescue<-setdiff(RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="879" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],AnyTreatmentRescue)

signgleRescues<-list(
  DEGs_XDP_FDR=DEGs_XDP_FDR$gene_name,
  onlydSVA_rescue=onlydSVA_rescue,
  only307_rescue=only307_rescue,
  only131_rescue=only131_rescue,
  only876_rescue=only876_rescue,
  only880_rescue=only880_rescue,
  only881_rescue=only881_rescue,
  only879_rescue=only879_rescue
)
Univ=as.vector(unlist(Reduce(union,list(signgleRescues))))
plotEnrichments(signgleRescues,signgleRescues,Univ)

XDPSignaturePathways_1<-DEGs_XDP_FDR_Enrichment[DEGs_XDP_FDR_Enrichment$p.value < 0.05,]
write_tsv(XDPSignaturePathways_1,file = "../XDPSignaturePathways_1.txt")

XDPSignaturePathways_1$RescuedBydSVA<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productdSVA)
XDPSignaturePathways_1$RescuedBy880<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product880)
XDPSignaturePathways_1$RescuedBy131<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product131)
XDPSignaturePathways_1$RescuedBy307<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product307)
XDPSignaturePathways_1$RescuedBy879<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product879)
XDPSignaturePathways_1$RescuedBy876<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product876)
XDPSignaturePathways_1$RescuedBy881<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product881)
XDPSignaturePathways_1$RescuedBy877<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product877)
XDPSignaturePathways_1$RescuedByAnyTreatment<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productAnyTreatment)
#XDPSignaturePathways_1$RescuedByAllTreatment<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productAllTreatment)
XDPSignaturePathways_1$RescuedByAnyASO<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productAnyASO)
XDPSignaturePathways_1$RescuedByAllASO<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productAllASO)

XDPSignaturePathways_1$RescuedBydSVA_numbers<-listLen(XDPSignaturePathways_1$RescuedBydSVA)
XDPSignaturePathways_1$RescuedBydSVA_perc<-(XDPSignaturePathways_1$RescuedBydSVA_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy880_numbers<-listLen(XDPSignaturePathways_1$RescuedBy880)
XDPSignaturePathways_1$RescuedBy880_perc<-(XDPSignaturePathways_1$RescuedBy880_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy131_numbers<-listLen(XDPSignaturePathways_1$RescuedBy131)
XDPSignaturePathways_1$RescuedBy131_perc<-(XDPSignaturePathways_1$RescuedBy131_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy307_numbers<-listLen(XDPSignaturePathways_1$RescuedBy307)
XDPSignaturePathways_1$RescuedBy307_perc<-(XDPSignaturePathways_1$RescuedBy307_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy879_numbers<-listLen(XDPSignaturePathways_1$RescuedBy879)
XDPSignaturePathways_1$RescuedBy879_perc<-(XDPSignaturePathways_1$RescuedBy879_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy876_numbers<-listLen(XDPSignaturePathways_1$RescuedBy876)
XDPSignaturePathways_1$RescuedBy876_perc<-(XDPSignaturePathways_1$RescuedBy876_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy881_numbers<-listLen(XDPSignaturePathways_1$RescuedBy881)
XDPSignaturePathways_1$RescuedBy881_perc<-(XDPSignaturePathways_1$RescuedBy881_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedBy877_numbers<-listLen(XDPSignaturePathways_1$RescuedBy877)
XDPSignaturePathways_1$RescuedBy877_perc<-(XDPSignaturePathways_1$RescuedBy877_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByAnyTreatment_numbers<-listLen(XDPSignaturePathways_1$RescuedByAnyTreatment)
XDPSignaturePathways_1$RescuedByAnyTreatment_perc<-(XDPSignaturePathways_1$RescuedByAnyTreatment_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
#XDPSignaturePathways_1$RescuedByAllTreatment_numbers<-listLen(XDPSignaturePathways_1$RescuedByAllTreatment)
#XDPSignaturePathways_1$RescuedByAllTreatment_perc<-(XDPSignaturePathways_1$RescuedByAllTreatment_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByAnyASO_numbers<-listLen(XDPSignaturePathways_1$RescuedByAnyASO)
XDPSignaturePathways_1$RescuedByAnyASO_perc<-(XDPSignaturePathways_1$RescuedByAnyASO_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByAllASO_numbers<-listLen(XDPSignaturePathways_1$RescuedByAllASO)
XDPSignaturePathways_1$RescuedByAllASO_perc<-(XDPSignaturePathways_1$RescuedByAllASO_numbers/XDPSignaturePathways_1$is_term_is_hit)*100

write_tsv(XDPSignaturePathways_1,file = "../XDPSignaturePathways_1_rescuedGenes_perTerm.txt")

XDPSignaturePathways_1$RescuedByOnlydSVA<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnlydSVA)
XDPSignaturePathways_1$RescuedByOnly880<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnly880)
XDPSignaturePathways_1$RescuedByOnly131<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnly131)
XDPSignaturePathways_1$RescuedByOnly307<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnly307)
XDPSignaturePathways_1$RescuedByOnly879<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnly879)
XDPSignaturePathways_1$RescuedByOnly876<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnly876)
XDPSignaturePathways_1$RescuedByOnly881<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,productOnly881)

XDPSignaturePathways_1$RescuedByOnlydSVA_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnlydSVA)
XDPSignaturePathways_1$RescuedByOnlydSVA_perc<-(XDPSignaturePathways_1$RescuedByOnlydSVA_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByOnly880_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnly880)
XDPSignaturePathways_1$RescuedByOnly880_perc<-(XDPSignaturePathways_1$RescuedByOnly880_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByOnly131_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnly131)
XDPSignaturePathways_1$RescuedByOnly131_perc<-(XDPSignaturePathways_1$RescuedByOnly131_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByOnly307_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnly307)
XDPSignaturePathways_1$RescuedByOnly307_perc<-(XDPSignaturePathways_1$RescuedByOnly307_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByOnly879_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnly879)
XDPSignaturePathways_1$RescuedByOnly879_perc<-(XDPSignaturePathways_1$RescuedByOnly879_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByOnly876_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnly876)
XDPSignaturePathways_1$RescuedByOnly876_perc<-(XDPSignaturePathways_1$RescuedByOnly876_numbers/XDPSignaturePathways_1$is_term_is_hit)*100
XDPSignaturePathways_1$RescuedByOnly881_numbers<-listLen(XDPSignaturePathways_1$RescuedByOnly881)
XDPSignaturePathways_1$RescuedByOnly881_perc<-(XDPSignaturePathways_1$RescuedByOnly881_numbers/XDPSignaturePathways_1$is_term_is_hit)*100

XDPSignaturePathways_1$RescuedBy3ASO<-apply(as.data.frame(XDPSignaturePathways_1$term_and_hit_genes),1,product3FAMASO)
XDPSignaturePathways_1$RescuedBy3ASO_numbers<-listLen(XDPSignaturePathways_1$RescuedBy3ASO)
XDPSignaturePathways_1$RescuedBy3ASO_perc<-(XDPSignaturePathways_1$RescuedBy3ASO_numbers/XDPSignaturePathways_1$is_term_is_hit)*100


XDPSignaturePathways_2<-XDPSignaturePathways_1 %>% slice_max(order_by = RescuedBy3ASO_numbers, n=20) 
write_tsv(XDPSignaturePathways_1,file = "../XDPSignaturePathways_3FAM_rescuedGenes_perTerm.txt")

Rescued_Sig_filtered <- XDPSignaturePathways_2 %>% filter(term_genes > 50, term_genes < 1000) %>% filter(RescuedBydSVA_perc > 30)
Rescued_Sig_filtered2 <- Rescued_Sig_filtered %>% filter(RescuedByAnyASO_perc > 30) %>% arrange(desc(RescuedByAnyASO_numbers))

plotEnrichments(rescuedGenes_list_FDR_Ztest,rescuedGenes_list_FDR_Ztest,Univ)
# dSVA        ASO307       ASO131        ASO876       ASO880        ASO881       ASO879
# dSVA   1.207404e-216  7.079454e-08 4.403838e-07  2.798676e-07 1.700709e-04  6.440633e-06 2.319843e-02
# ASO307  7.079454e-08  0.000000e+00 2.178331e-44 3.117766e-103 5.777564e-44  2.578402e-21 1.080337e-29
# ASO131  4.403838e-07  2.178331e-44 0.000000e+00  5.213547e-54 1.002226e-44  2.476470e-20 3.668458e-24
# ASO876  2.798676e-07 3.117766e-103 5.213547e-54  0.000000e+00 2.391345e-41  1.997090e-25 2.249128e-29
# ASO880  1.700709e-04  5.777564e-44 1.002226e-44  2.391345e-41 0.000000e+00  5.000350e-14 2.455626e-30
# ASO881  6.440633e-06  2.578402e-21 2.476470e-20  1.997090e-25 5.000350e-14 1.182698e-282 5.305191e-36
# ASO879  2.319843e-02  2.073410e-30 3.668458e-24  2.249128e-29 2.455626e-30  6.246631e-37 0.000000e+00

common_elements <- list()
for (i in 1:nrow(XDPSignaturePathways_1)) {
  common_elements[[i]] <- intersect(XDPSignaturePathways_1$RescuedBydSVA[[i]], XDPSignaturePathways_1$RescuedBy307[[i]])
}
names(common_elements)<-XDPSignaturePathways_1$term_name

common_elements_1 <- list()
for (i in 1:nrow(XDPSignaturePathways_1)) {
  common_elements_1[[i]] <- intersect(XDPSignaturePathways_1$RescuedBydSVA[[i]], XDPSignaturePathways_1$RescuedBy131[[i]])
}
names(common_elements_1)<-XDPSignaturePathways_1$term_name

# names(common_elements[lengths(common_elements) >= 10])
# [1] "GOCC_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE"                 
# [2] "Mouse_Mef2cKO_DEGs"                                          
# [3] "Mouse_Mef2cKO_UpRegDEGs"                                     
# [4] "GOCC_ANCHORING_JUNCTION"                                     
# [5] "NABA_MATRISOME"                                              
# [6] "GOBP_POSITIVE_REGULATION_OF_MULTICELLULAR_ORGANISMAL_PROCESS"
# [7] "ZIM3_TARGET_GENES"                                           
# [8] "GOBP_NEGATIVE_REGULATION_OF_SIGNALING"                       
# [9] "GOBP_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY"       
# [10] "GOBP_POSITIVE_REGULATION_OF_MOLECULAR_FUNCTION"              

View(DEGs_XDP_FDR_Enrichment[DEGs_XDP_FDR_Enrichment$p.value < 0.05 & DEGs_XDP_FDR_Enrichment$enrichment_zscore > 3 & (DEGs_XDP_FDR_Enrichment$term_genes > 50 & DEGs_XDP_FDR_Enrichment$term_genes < 1000),c("term_name","p.value")])

#XDP_Path<-c("GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE","GOBP_REGULATION_OF_MORPHOGENESIS_OF_AN_EPITHELIUM","GOCC_SYNAPTIC_VESICLE_MEMBRANE","GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES","GOBP_POSITIVE_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY","GOBP_REGULATION_OF_PEPTIDYL_SERINE_PHOSPHORYLATION","GOCC_EXOCYTIC_VESICLE","GOMF_BETA_CATENIN_BINDING","GOBP_POSITIVE_REGULATION_OF_FIBROBLAST_PROLIFERATION","HP_POOR_FINE_MOTOR_COORDINATION","GOBP_NEURAL_CREST_CELL_DIFFERENTIATION","GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT")
XDP_Path<-c("GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE","GOBP_REGULATION_OF_MORPHOGENESIS_OF_AN_EPITHELIUM","GOCC_SYNAPTIC_VESICLE_MEMBRANE","GOBP_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_ADHESION_MOLECULES","GOBP_POSITIVE_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY","GOBP_REGULATION_OF_PEPTIDYL_SERINE_PHOSPHORYLATION","GOCC_EXOCYTIC_VESICLE","HP_POOR_FINE_MOTOR_COORDINATION","GOBP_NEURAL_CREST_CELL_DIFFERENTIATION","GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT")

XDP_signature_pathways_finaldata<-XDPSignaturePathways_1[XDPSignaturePathways_1$term_name %in% XDP_Path,]
saveRDS(XDP_signature_pathways_finaldata, "../XDP_signature_pathways_finaldata.rds")

XDP_signature_pathways_finaldata<-readRDS("../XDP_signature_pathways_finaldata.rds")
## Figure 4A ##
res<-read.table("NO.CON+NO.XDP-LFC+Var.tsv",sep = "\t",header = T)
lab_italics <- paste0("italic('", res$gene_name, "')")

res1<-res[order(res$fdr),]

res<-res[res$lfc <4,]

pdf("/Users/ry077/Dropbox (Partners HealthCare)/XDP/ASOpaper2023/FigureJune2025/Fig4A.pdf",height = 6,width = 6)
EnhancedVolcano(res,
                lab = res$gene_name,
                x = 'lfc',
                y = 'fdr',
                selectLab = c("CAT","TAF1","TNFRSF11A","ZIC2","ALK","PRKCQ","ZNF300","CHCHD2","NLRP2"),
                xlab = bquote(~Log~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'FDR'),
                pCutoff = 10e-06,
                #  FCcutoff = 2.0,
                pointSize = 2.0,
                labSize = 5.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                col=c('black', 'black','red3','red3'),
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'None',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "",
                subtitle = "")
dev.off()

XDP_signature_pathways_finaldata<-XDP_signature_pathways_finaldata[order(XDP_signature_pathways_finaldata$is_term_is_hit,decreasing = F),]
XDP_signature_pathways_finaldata$logP<-round(-log10(XDP_signature_pathways_finaldata$p.value))

XDP_signature_pathways_finaldata$Temp_print<-gsub("stat","STAT",gsub("dna","DNA",gsub("erk","ERK",str_to_sentence(gsub("_"," ",gsub("HP_","",gsub("GO[BP|MF|CC]*_","",XDP_signature_pathways_finaldata$term_name)))))))

XDP_signature_pathways_finaldata$Temp_print<-factor(XDP_signature_pathways_finaldata$Temp_print,level=XDP_signature_pathways_finaldata$Temp_print)


pdf("/Users/ry077/Dropbox (Partners HealthCare)/XDP/ASOpaper2023/FigureJune2025/Fig4B.pdf",height = 6,width = 6)
ggplot(data=XDP_signature_pathways_finaldata, aes(x=as.factor(Temp_print), y=is_term_is_hit,fill=logP)) +
  geom_bar(stat="identity", width=0.8,position = "dodge") +
  scale_fill_continuous(low = "#1B9E77", high = "#105e47") +
  theme_bw() +
  # remove gridlines. panel.grid.major is for verical lines, panel.grid.minor is for horizontal lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # remove borders
        panel.border = element_blank(),
        # removing borders also removes x and y axes. Add them back
        axis.line = element_line(),strip.text.y = element_text(angle = 180), 
        #legend.position = "bottom",
        axis.text.y = element_text(size = 12)) +
  coord_flip() +  scale_x_discrete(labels = label_wrap(width = 40))+
  # remove x-axis label
  xlab("XDP signature genes functional enrichment (p < 0.05)") +
  # change the y-axis label
  ylab("Gene counts") + labs(fill = bquote(~-log[10]~ 'p-value'))
dev.off()

## Figure new 5B (old 4C) ##
Plot_rescueDF<-XDP_signature_pathways_finaldata[XDP_signature_pathways_finaldata$term_name %in% c("GOBP_NEURAL_CREST_CELL_DIFFERENTIATION","GOBP_POSITIVE_REGULATION_OF_NF_KAPPAB_TRANSCRIPTION_FACTOR_ACTIVITY","GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE","GOBP_GLIOGENESIS","GOCC_PRESYNAPSE","GOBP_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY"),c(2,8,30,36,34,32)]
colnames(Plot_rescueDF)<-c("Pathway","XDP Signatures","dSVA_rescued","ASO307_rescue","ASO131_rescued","ASO880_rescued")
Plot_rescueDF_M<-melt(Plot_rescueDF)


#Plot_rescueDF_M$Pathway<-factor(Plot_rescueDF_M$Pathway,levels = c("GOBP_POSITIVE_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY","GOBP_REGULATION_OF_PEPTIDYL_SERINE_PHOSPHORYLATION","GOBP_NEURAL_CREST_CELL_DIFFERENTIATION"
                                                         ,"GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT") )


pdf("/Users/ry077/Dropbox (Partners HealthCare)/XDP/ASOpaper2023/FigureJune2025/Fig4D_4.pdf",height = 8,width = 6)
ggplot(Plot_rescueDF_M, aes(x = reorder(variable,value), y = value, fill = variable)) + 
  ylab("Gene counts") +
  xlab("XDP associated functions") +
  geom_col() + 
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~ Pathway, strip.position = "left", nrow = 4,
             labeller = as_labeller(c(GOBP_NEURAL_CREST_CELL_DIFFERENTIATION = "Neural crest\ncell differentiation",
                                      GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT = "Receptor signaling\npathway via STAT",
                                      GOBP_REGULATION_OF_PEPTIDYL_SERINE_PHOSPHORYLATION = "Regulation of peptidyl\nserine phosphorylation",
                                      GOBP_POSITIVE_REGULATION_OF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY = "Regulation of DNA\nbinding TF activity"))) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))
dev.off()

# new 4A
library(ggplot2)
brewer.pal(8, "Dark2")
colours<-c("#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666")
# Create a data frame in "long" format
x = c("dSVA", "ASO880", "ASO131")
rescued_genes = c(456,472,446)
XDP_signature_genes = c(1490-456,1490-472,1490-446)

#rank_data<-as.data.frame(Treatment = c("dSVA", "ASO880", "ASO131"),
 #            XDP_signature_genes = c(1490,1490,1490), 
#             p-value=c())
to_plot <- data.frame(Treatment=x,rescued_genes=rescued_genes,XDP_signature_genes=XDP_signature_genes)
melted<-melt(to_plot, id="Treatment")

# Create a new column to define the fill color for the 'Incorrect' portions
# All 'Correct' portions will have the same fill color.
# All 'Incorrect' portions will have a unique fill color based on the student.
melted <- melted %>%
  mutate(Gene_Type = ifelse(
    variable == "XDP_signature_genes",
    "XDP_signature_genes",
    paste0("Rescued by ", Treatment)
  ))
melted$Gene_Type<-factor(melted$Gene_Type,levels=c("XDP_signature_genes","Rescued by dSVA","Rescued by ASO880","Rescued by ASO131"))
melted$Treatment<-factor(melted$Treatment,levels = c("dSVA","ASO880","ASO131"))
melted$variable<-factor(melted$variable,levels = c("XDP_signature_genes","rescued_genes"))

pdf("/Users/ry077/Dropbox (Partners HealthCare)/XDP/ASOpaper2023/FigureSep2025/Fig4A_piecharts.pdf",height = 6,width = 8)
# Create the stacked bar chart with different 'Incorrect' fill colors
#ggplot(melted, aes(x = Treatment, y = value, fill = variable)) +
  # Create the stacked bars
#  geom_bar(position = "fill", stat = "identity", width = 0.7) +  
  # Remove position = "fill" for counts
ggplot(melted, aes(x = "", y = value, fill = variable)) +
  geom_bar( stat = "identity", width = 1,color="black") +
 # scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = c( "grey","blue")) +
  coord_polar("y", start = 0) +
  facet_grid(.~Treatment)+
  labs(
    title = "",
    x = "",
   y = "Treatments",
    fill = "Gene type"
  ) +
  # Add percentage labels to the segments
  theme_minimal() +
  # Add percentage labels only to the "Correct" segments
  geom_text(
    aes(
      label = ifelse(variable == "rescued_genes",
                     scales::percent(value/1490),
                     "")),
  #  position = position_fill(vjust = 500), # Center labels vertically
  #  size = 4) +
  #  geom_text(aes(label = paste0(percentage, "%")),
              position = position_stack(vjust = 0.5),
              color = "black", size = 4) +
 # geom_text(
  #  data = rank_data,
  #  aes(x = student, y = total_count, label = paste("p-value:", rank), fill = NULL),
  #  vjust = -0.5,
  #  size = 4,
  #  color = "black"
  #) +
  # Remove the legend for the text color
  guides(color = "none") +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))
dev.off()



## data ##
rescuedGenes_list_Ztest<-list(
  DEGs_XDP_FDR=DEGs_XDP_FDR$gene_name,
  dSVA=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="dSVA" & RescuedGenes_DF$isRescued.2ZTests == TRUE],
  ASO307=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="307" & RescuedGenes_DF$isRescued.2ZTests == TRUE],
  ASO131=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="131" & RescuedGenes_DF$isRescued.2ZTests == TRUE],
  ASO876=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="876" & RescuedGenes_DF$isRescued.2ZTests == TRUE],
  # ASO877=RescuedGenes_FDR$gene_name[RescuedGenes_FDR$Treatment=="877" & RescuedGenes_FDR$isRescued.2ZTests == TRUE],
  ASO880=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="880" & RescuedGenes_DF$isRescued.2ZTests == TRUE],
  ASO881=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="881" & RescuedGenes_DF$isRescued.2ZTests == TRUE],
  ASO879=RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment=="879" & RescuedGenes_DF$isRescued.2ZTests == TRUE]
)
pdf("/Volumes/Macintosh HD/Users/ry077/Dropbox (Partners HealthCare)/XDP/ASOpaper2023/V1Figures/XDR_FDR_ASO2Ztets_rescue_overlaps_heatmap.pdf",height = 4, width = 4)
Univ=as.vector(unlist(Reduce(union,list(rescuedGenes_list_Ztest))))
plotEnrichments(rescuedGenes_list_Ztest,rescuedGenes_list_Ztest,Univ)
dev.off()

# Functions ##
load_enrichment_data <- function(file_pattern){
  '/Volumes//talkowski/Samples/XDP/ASO_dSVA/Results_draft1/2025-06-18_enrichment.results' %>%
    list.files(
      pattern=file_pattern,
      full.names=TRUE
    ) %>%
    tibble(filepath=.) %>%
    mutate(
      fileinfo=
        filepath %>%
        basename() %>% 
        str_remove('.tsv')
    ) %>%
    separate_wider_delim(
      fileinfo,
      delim='-',
      names=
        c(
          NA,
          NA,
          NA,
          'Rescue.Criteria',
          'OverlapList'
        )
    ) %>% 
    rowwise() %>% 
    mutate(enrichment.results=list(read_tsv(filepath))) %>% 
    ungroup() %>% 
    unnest(enrichment.results) %>% 
    dplyr::select(-c(filepath))
}

productdSVA<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "dSVA" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product131<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "131" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product307<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "307" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product876<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "876" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product877<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "877" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product879<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "879" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product880<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "880" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}
product881<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(RescuedGenes_DF$gene_name[RescuedGenes_DF$Treatment == "880" & RescuedGenes_DF$isRescued.2ZTests == TRUE],y)
  return(z)
}

AnyTreatmentRescue<-Reduce(union,rescuedGenes_list_FDR_Ztest)
productAnyTreatment<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(AnyTreatmentRescue,y)
  return(z)
}
AllTreatmentRescue<-Reduce(intersect,rescuedGenes_list_FDR_Ztest)
productAllTreatment<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(AllTreatmentRescue,y)
  return(z)
}
AnyASORescue<-Reduce(union,rescuedGenes_list_FDR_Ztest)
productAnyASO<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(AnyASORescue,y)
  return(z)
}
AllASORescue<-Reduce(intersect,rescuedGenes_list_FDR_Ztest)
productAllASO<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(AllASORescue,y)
  return(z)
}


productOnlydSVA<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(onlydSVA_rescue,y)
  return(z)
}
productOnly131<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(only131_rescue,y)
  return(z)
}
productOnly307<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(only307_rescue,y)
  return(z)
}
productOnly876<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(only876_rescue,y)
  return(z)
}
productOnly879<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(only879_rescue,y)
  return(z)
}
productOnly880<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(only880_rescue,y)
  return(z)
}
productOnly881<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(only881_rescue,y)
  return(z)
}

AnyTreatmentRescue_3Fam<-Reduce(union,rescuedGenes_list_FDR_Ztest)
product3FAMASO<-function(x){
  y=strsplit(x,",")[[1]]
  z=intersect(AnyTreatmentRescue_3Fam,y)
  return(z)
}
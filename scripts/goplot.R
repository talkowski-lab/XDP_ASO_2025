library(ggplot2)
library(ggrepel)
library(GOstats)
library(GO.db)
library(annotate)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

source("/data/talkowski/dg520/app/lib/R/geom_noBG.R")


gotest = function(target_gene_symbols, background_gene_symbols, species = "human", ont = "BP") {
  if (species=="human"){
    annot="org.Hs.eg.db"
    db = toTable(org.Hs.egALIAS2EG)
  }else if (species=="mouse"){
    annot="org.Mm.eg.db"
    db = toTable(org.Mm.egALIAS2EG)
  }else{
    stop("'species' has to be either 'human' or 'mouse'.")
  }

  genes = db[which(db$alias_symbol %in% target_gene_symbols), 1]
  target_gene_found = which(target_gene_symbols %in% db$alias_symbol)
  message(paste0("Found ",length(target_gene_found)," out of ",length(target_gene_symbols)," (",round(100*length(target_gene_found)/length(target_gene_symbols),2),"%) target genes in the GO database."))
  universeGenes = db[which(db$alias_symbol %in% background_gene_symbols), 1]
  background_gene_found = which(background_gene_symbols %in% db$alias_symbol)
  message(paste0("Found ",length(background_gene_found)," out of ",length(background_gene_symbols)," (",round(100*length(background_gene_found)/length(background_gene_symbols),2),"%) background genes in the GO database."))

  params = new("GOHyperGParams",
                geneIds = genes,
                universeGeneIds = universeGenes,
                annotation = annot,
                ontology = ont,
                pvalueCutoff = 0.05,
                testDirection = "over",
                conditional = TRUE)
  hgOver = hyperGTest(params)
  s = summary(hgOver, pvalue=1)
  pvals.adj = p.adjust(s[,2], method = "BH")
  s2 = cbind(ont,s, pvals.adj)
  colnames(s2) = c("Class","GO.ID","p.value","OR","ExpCount","Count","Size","Terms","FDR")
  s2 = s2[,c("Class","Terms","p.value","FDR","GO.ID","OR","ExpCount","Count","Size")]
  # sort results by adjusted p-value
  s2 = s2[order(s2[,"FDR"]),]
  return(s2)
}


golist = function(tb,goclass="BF"){
  gt=read.table(tb,header=T,sep="\t")
  gt=cbind(goclass,gt[,c("Description","P.value","FDR.q.value")])
  colnames(gt)=c("Class","Terms","p.value","FDR")
  return(gt)
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}


goplot = function(pos,cutoff=0.05,neg=NULL,method = "FDR"){
  if (!(method %in% c("FDR","p.value"))){stop("'method' has to be either 'FDR' or 'p.value'")}
  if (method=="FDR"){
    ylabel="FDR (-log10-transformed)"
  }else if (method=="p.value"){
    ylabel="Raw p value (-log10-transformed)"
  }
  go=pos
  colnames(go)=c("Class","Terms","p.value","FDR")
  go=go[order(go[,method],decreasing = T),]
  #go=go[go$FDR<cutoff,]
  go$transformedFDR=-log10(go[,method])
  vv=as.vector(go$Terms)
  vv=firstup(vv)
  go$Terms=vv
  go$Terms=factor(go$Terms,levels=vv)
  #go$direct="Nondirectional"
  ggcut=geom_hline(yintercept = -log10(cutoff),linetype="dashed",color="black")
  if (is.null(neg)){
    p=ggplot(go,aes(x=Terms,y=transformedFDR,fill="#FF6666"))+geom_bar(stat="identity",position='dodge')+xlab("GO")+ylab(ylabel)+coord_flip()+ggcut+geom_noBG()
  }else{
    go$direct="Upregulation"
    go1=neg
    colnames(go1)=c("Class","Terms","FDR")
    go1=go1[order(go[,method],decreasing = T),]
    #go1=go1[go1$FDR<cutoff,]
    go1$transformedFDR=log10(go1[,method])
    vv1=as.vector(go1$Terms)
    vv1=firstup(vv1)
    go1$Terms=vv1
    go1$Terms=factor(go1$Terms,levels=vv1)
    go1$direct="Downregulation"
    go=rbind(go,go1)
    ggcut=geom_hline(yintercept = c(log10(cutoff),-log10(cutoff)),linetype="dashed",color="black")
    p=ggplot(go,aes(x=Terms,y=transformedFDR,fill=direct))+geom_bar(stat="identity",position='dodge')+geom_hline(yintercept = 0,linetype="solid",color="black",size=1)+geom_vline(xintercept = 0,linetype="solid",color="black",size=1)+xlab("GO")+ylab(ylabel)+coord_flip()+ggcut+geom_noBG()
  }
  
  return(p)
}

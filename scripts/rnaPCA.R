library(matrixStats)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(corrplot)


rnaPCA <- function (object, sampleTable, intgroup, textLabel = NULL, xPC = 1,yPC = 2, ntop = 500, pointSize = 3, labelSize = 3, labelOverlaps = nrow(sampleTable), returnData = FALSE) 
{
  mat = as.matrix(object[,rownames(sampleTable)])
  rv <- rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(mat[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(sampleTable))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(sampleTable[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = "_"))
  }
  else {
    sampleTable[[intgroup]]
  }
  d <- data.frame(pca$x, group = group, sampleTable)
  colnames(d)[1:dim(pca$x)[2]]=paste0("PC",c(1:dim(pca$x)[2]))
  attr(d, "percentVar") <- percentVar
  #d <- data.frame(xPC = pca$x[, 1], yPC = pca$x[, 2], group = group, intgroup.df, name = colData(rld)[,1])
  #plot(p)
  if (returnData) {
    return(d)
  } else {
    p0 = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = axis_line, axis.ticks=axis_line, axis.ticks.length=tick_length, axis.title=axis_text, axis.text=tick_text)
    if (is.null(textLabel)){
	  p = ggplot(data = d, aes(x = .data[[paste0("PC", xPC)]], y = .data[[paste0("PC", yPC)]], color = .data[["group"]] )) + geom_point(size = pointSize) + xlab(paste0("PC", xPC,": ", round(attr(d, "percentVar")[xPC] * 100), "% variance")) + ylab(paste0("PC", yPC,": ", round(attr(d, "percentVar")[yPC] * 100), "% variance")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + p0
	} else {
	  p = ggplot(data = d, aes(x = .data[[paste0("PC", xPC)]], y = .data[[paste0("PC", yPC)]], color = .data[["group"]], label = .data[[textLabel]])) + geom_point(size = pointSize) + xlab(paste0("PC", xPC,": ", round(attr(d, "percentVar")[xPC] * 100), "% variance")) + ylab(paste0("PC", yPC,": ", round(attr(d, "percentVar")[yPC] * 100), "% variance")) + geom_text_repel(size = labelSize, max.overlaps = labelOverlaps) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + p0
    }
	#plot(p)
    return(p)
  }
}


corPCA <- function(rnaPCA_data, cols = NULL, pcs = 1:10, returnData = FALSE){
  if (is.null(cols)){
    new = cbind(rnaPCA_data[, grep("^PC",colnames(rnaPCA_data))][, pcs], rnaPCA_data[, which(!grepl("^PC",colnames(rnaPCA_data)))])
  } else {
    new = cbind(rnaPCA_data[, grep("^PC",colnames(rnaPCA_data))][, pcs], rnaPCA_data[, cols])
  }

  new = as.data.frame(unclass(new),stringsAsFactors = T)
  new[, sapply(new, is.factor)] = sapply(new[, sapply(new, is.factor)], unclass)
  cm = rcorr(as.matrix(new))
  cm_p = cor.mtest(as.matrix(new))
  corrplot(cm$r, type="upper", order="original",p.mat = cm_p$p, sig.level = 0.05, insig = "blank")
  if (returnData){
    return(new)
  }
}

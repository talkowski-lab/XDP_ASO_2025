library(ggplot2)
source("~/lib/R/geom_noBG.R")


stat_corplot = function(mat1,
                        mat2, 
                        val1,    # "LFC.GenotypeXDP"
                        val2,    # "LFC.GenotypeXDP"
                        filter1, # "FDR.GenotypeXDP"
                        filter2, # "FDR.GenotypeXDP"
                        name1,   # "Current"
                        name2,   # "CellPaper"
                        name3 = "Both",
                        cutoff1 = 0.1,
                        cutoff2 = 0.1,
                        point_size = 0.5,
                        highlight_size = 2,
                        xr = c(-3, 3),
                        yr = c(-3, 3),
                        content = "DEG Fold Changes"){
  
  shared = intersect(rownames(mat1), rownames(mat2))
  df = data.frame(val1 = mat1[shared, val1] , val2 = mat2[shared, val2], filter1 = mat1[shared, filter1], filter2 = mat2[shared, filter2])
  rownames(df) = shared
  df[is.na(df)] = 1
  df$Group = "None"
  df[(df$filter1 < cutoff1) & (df$filter2 >= cutoff2), "Group"] = name1
  df[(df$filter1 >= cutoff1) & (df$filter2 < cutoff2), "Group"] = name2
  df[(df$filter1 < cutoff1) & (df$filter2 < cutoff2), "Group"] = name3
  df = df[df$Group != "None", ]
  df$Group = factor(df$Group, levels = c(name1, name2, name3))
  df$size = point_size
  df[df$Group == name3, "size"] = highlight_size
  plt = ggplot(df, aes(x = val1, y = val2, color = Group)) + geom_point(size = df$size) + 
    geom_color(c("orange", "black", "red")) + 
    xlim(xr) + 
    ylim(yr) + 
    geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_abline(slope = 1, color = "blue") + 
    xlab(paste0(content, " in ", name1)) + 
    ylab(paste0(content, " in ", name2)) + 
    geom_noBG()
  
  plot(plt)
  
  return(plt)
}

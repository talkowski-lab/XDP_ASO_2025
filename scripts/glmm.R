library(lme4)
library(pbkrtest)
library(DESeq2)
library(Hmisc)
library(BSDA)
library(aod)


get_ci = function(res, beta = NULL, stderr = NULL, pval = NULL, fdr = NULL, term, class,  ci = 0.95){
  if (!is.null(beta)){
    beta = read.table(beta, header = T, row.names = 1, sep = "\t", check.names = F)
    std  = read.table(stderr, header = T, row.names = 1, sep = "\t", check.names = F)
    pv  = read.table(pval, header = T, row.names = 1, sep = "\t", check.names = F)
    bh = read.table(fdr, header = T, row.names = 1, sep = "\t", check.names = F)
  } else {
    beta = as.data.frame(res@beta)
    std = as.data.frame(res@stderr)
    pv  = as.data.frame(res@pval)
    bh = as.data.frame(res@fdr)
  }
  
  alpha = 1 - ci
  k = length(which(bh[[term]] < alpha))
  o = length(which(pv[[term]] < alpha))
  cl = 1 - (k / o) * alpha / 2
  margin = apply(std[,term, drop = F], 1, FUN = function(x){return(qnorm(cl) * x)})
  df = as.data.frame(cbind(LFC = beta[,term], margin = margin))
  df[["name"]] = as.vector(rownames(df))
  df[["class"]] = class
  rownames(df) = NULL
  
  return(df)
}


get_fdr = function(res, beta = NULL, pval = NULL, fdr = NULL, term, class, keyword = NULL){
  if (!is.null(beta)){
    fc = read.table(beta, header = T, row.names = 1, sep = "\t", check.names = F)
    pv  = read.table(pval, header = T, row.names = 1, sep = "\t", check.names = F)
    bh = read.table(fdr, header = T, row.names = 1, sep = "\t", check.names = F)
  } else {
    fc = as.data.frame(res@beta)
    pv  = as.data.frame(res@pval)
    bh = as.data.frame(res@fdr)
  }
  
  new0= fc[,grep(term, colnames(fc)), drop = F]
  colnames(new0) = paste0("LFC.",colnames(fc)[grep(term, colnames(fc))])
  new1= pv[,grep(term, colnames(pv)), drop = F]
  colnames(new1) = paste0("pval.",colnames(pv)[grep(term, colnames(pv))])
  new2 = bh[,grep(term, colnames(bh)), drop = F]
  colnames(new2) = paste0("FDR.",colnames(bh)[grep(term, colnames(bh))])
  if (!is.null(keyword)){
    new0 = new0[keyword, , drop = F]
    new1 = new1[keyword, , drop = F]
    new2 = new2[keyword, , drop = F]
  }
  new = as.data.frame(new2)
  new = cbind(new, new0, new1)
  #colnames(new) = paste0("V",1:ncol(new))
  new$name = as.vector(rownames(new))
  new$class = class
  rownames(new) = NULL
  
  return(new)
}


get_intron_length = function(irnames){
  start = as.numeric(sapply(strsplit(irnames, "/"), "[[", 2))
  end   = as.numeric(sapply(strsplit(irnames, "/"), "[[", 3))
  return(end-start+1)
}


cross_model_z_test = function(term1, term2 = NULL, n1 = NULL, cond1 = NULL, lvl1, beta1, stderr1, n2 = NULL, cond2 = NULL, lvl2, beta2, stderr2, class = NULL){
  if (is.null(term2)){term2 = term1}
  if (is.null(n1)){
    if (is.null(cond2)){stop("Need to provide either 'n1' (integer) or 'cond1' (colData file name).")}
    grp1 = read.table(cond1, header = T, row.names = 1, sep = "\t", check.names = F)
    n1 = nrow(grp1[grp1[[term1]] == lvl1, ])
  }
  if (is.null(n2)){
    if (is.null(cond2)){stop("Need to provide either 'n2' (integer) or 'cond2' (colData file name).")}
    grp2 = read.table(cond2, header = T, row.names = 1, sep = "\t", check.names = F)
    n2 = nrow(grp2[grp2[[term2]] == lvl2, ])
  }
  b1 = read.table(beta1, header = T, row.names = 1, sep = "\t", check.names = F)
  b2 = read.table(beta2, header = T, row.names = 1, sep = "\t", check.names = F)
  shared_genes = intersect(rownames(b1), rownames(b2))
  bt1 = b1[shared_genes, paste0(term1, lvl1)]
  bt2 = b2[shared_genes, paste0(term2, lvl2)]
  e1 = read.table(stderr1, header = T, row.names = 1, sep = "\t", check.names = F)
  e2 = read.table(stderr2, header = T, row.names = 1, sep = "\t", check.names = F)
  er1 = e1[shared_genes, paste0(term1, lvl1)]
  er2 = e2[shared_genes, paste0(term2, lvl2)]
  if (is.null(bt1)){bt1 = rep(0, length(shared_genes))}
  if (is.null(bt2)){bt2 = rep(0, length(shared_genes))}
  if (is.null(er1)){er1 = rep(0, length(shared_genes))}
  if (is.null(er2)){er2 = rep(0, length(shared_genes))}
  sd1 = er1 * sqrt(n1)
  sd2 = er2 * sqrt(n2)
  tt = t(sapply(1:length(shared_genes), FUN = function(x){tmp = zsum.test(mean.x = bt1[x], mean.y = bt2[x], sigma.x = sd1[x], sigma.y = sd2[x], n.x = n1, n.y = n2, mu = 0); return(c(tmp$p.value, mean(tmp$conf.int)))}))
  tt = as.data.frame(tt)
  colnames(tt) = c(paste0("pval.", term1, lvl1), paste0("LFC.", term1, lvl1))
  tt[[paste0("FDR.", term1, lvl1)]] = p.adjust(tt[[paste0("pval.", term1, lvl1)]], method = "BH")
  tt = tt[, c(paste0("FDR.", term1, lvl1), paste0("LFC.", term1, lvl1), paste0("pval.", term1, lvl1))]
  tt$name = shared_genes
  if (is.null(class)){
    tt$class = paste0(term1, lvl1, "_vs_", term2, lvl2)
  } else {
    tt$class = class
  }
  
  return(tt)
}



get_se_from_contrast = function(vcov_mat, contrast_vector){
  new_var = t(contrast_vector) %*% tcrossprod(as.matrix(vcov_mat), t(contrast_vector))
  new_se = sqrt(new_var)

  return(new_se)
}


wald_test_from_contrast = function(object, beta = NULL, var_cov = NULL, contrast = NULL, term1, term2 = NULL, lvl1, lvl2, vcov_name = "name"){
  if (!is.null(beta)){
    b = read.table(beta_file, header = T, row.names = 1, sep = "\t", check.names = F)
    v = read.table(vcov_file, header = T, sep = "\t", check.names = F)
  } else {
    b = object@beta
    v = object@vcov
  }

  if (is.null(contrast)){
    if (is.null(term2)){term2 = term1}
    contrast_vec = rep(0, ncol(b))
    idx1 = which(colnames(b) == paste0(term1, lvl1))
    idx2 = which(colnames(b) == paste0(term2, lvl2))
    contrast_vec[idx1] = 1
    contrast_vec[idx2] = -1
  } else {
    contrast_vec = contrast
  }
  
  new_betas = as.matrix(b) %*% contrast_vec     # `new_betas` is a 1-column matrix

  vv0 = split(v, v[[vcov_name]])                 
  vv  = lapply(vv0, function(x){return(x[, -ncol(x)])})            #`vv` is a list, each element is a vcov data frame of a gene's model
  new_ses0 = lapply(vv, get_se_from_contrast, contrast_vector = contrast_vec) 
  new_ses = do.call(rbind, new_ses0)           # `new_ses` is a 1-column matrix

  new_z = new_betas / new_ses                  # `new_z` is a 1-column matrix
  new_p = apply(new_z, 1, FUN = function(x){return(2*pnorm(q = abs(x), lower.tail = FALSE))})    # `new_p` is a vector
  new_q = p.adjust(new_p, method = "BH")                                                         # `new_q` is a vector

  tt = data.frame(new_q, new_betas[, 1], new_p)
  colnames(tt) = c(paste0("FDR.", term1, lvl1), paste0("LFC.", term1, lvl1), paste0("pval.", term1, lvl1))
  tt$name = rownames(b)
  tt$class = paste0(term1, lvl1, "_vs_", term2, lvl2)
  
  return(tt)
}


overlap_hgtest = function(group1, group2, background, overlaps = NULL,  as_is = F){
  background = unique(background)
  if (as_is){
    g1 = unique(group1)
    g2 = unique(group2)
  } else {
    g1 = unique(group1[which(group1 %in% background)])
    g2 = unique(group2[which(group2 %in% background)])
  }
  if (is.null(overlaps)){
    ov = intersect(g1, g2)
  } else {
    ov = unique(overlaps)
    if (!as_is){
      ov = ov[which(ov %in% background)]
    }
  }
  
  p = phyper(length(ov) - 1, length(g1), length(background) - length(g1), length(g2), lower.tail=FALSE)
  
  return(p)
}


meta_p = function(sample_list, col_feature = NULL, col_pval, col_fc, padj_method = "bonferroni"){
  extract_func = function(x, i){return(unname(as.vector(x[, i])))}
  sel_func = function(x, i, ref){if (i > 0) {ft = unname(as.vector(x[, i]))} else {ft = rownames(x)}; rownames(x) = ft; return(x[ref, ])}
  metap_func = function(x){return(pchisq(-2 * sum(log(x)), df = length(x), lower.tail = F))}
  dir_func = function(x){up = all(x > 0); dn = all(x < 0); dir_ch = "unknown"; if ((!is.na(up)) & up){dir_ch = "up"}; if ((!is.na(dn)) & dn){dir_ch = "down"}; return(dir_ch)}
  
  if (!is.numeric(col_feature)){
    stop("'col.feature' must be provided as a single numeric value, indicating which column to look for feature names. Provide 0 to use row names.")
  } else if (col_feature == 0){
    feature_list = lapply(sample_list, rownames)
  } else {
    feature_list = lapply(sample_list, extract_func, col_feature)
  }
  shared = Reduce(intersect, feature_list)
  if (length(shared) == 0){stop("Error: no overlapped feature found across input.")}
  sample_list1 = lapply(sample_list, sel_func, col_feature, shared)
  
  pval_mat = sapply(sample_list1, extract_func, col_pval)
  pv_vec = apply(pval_mat, 1, metap_func)
  fc_mat   = sapply(sample_list1, extract_func, col_fc)
  fc_vec = apply(fc_mat, 1, dir_func)
  
  df = data.frame(meta_p = pv_vec, direction = fc_vec)
  rownames(df) = shared
  df$meta_fdr = p.adjust(df$meta_p, method = padj_method)
  
  return(df)
}


search_fc_cutoff_go = function(mat, col_fc, fc_cutoffs, col_p, p_cutoff, go_list, transform = T, link_func = log, keyword = NULL, gs_ref = "GRCh37.75_GRCm38.83_ensemblID_geneSymbol_pairs.txt", return_data = F){
  if (is.null(rownames(mat))){stop("'mat' must have row names as ensembl IDs.")}
  if (!(col_fc %in% colnames(mat))){stop("'col_fc' cannot be found in the column names of 'mat'.")}
  
  ref = read.table(gs_ref, header = F, row.names = 1, sep = "\t")
  
  df = c()
  for (x in fc_cutoffs){
    message("Using FC cutoff: ", x)
    x_t = x
    if (transform){x_t = link_func(x)}
    deg = rownames(mat[which((abs(mat[[col_fc]]) >= x_t) & (mat[[col_p]] < p_cutoff)), ])
    deg_gs = ref[deg, 1]
    bg_gs = ref[rownames(mat), 1]
    deg_go = go_test(deg_gs, bg_gs, go_list)
    min_fdr = min(deg_go$fdr)
    n_sig = length(which(deg_go$fdr < 0.1))
    df = rbind(df, c(x, min_fdr, n_sig))
  }
  colnames(df) = c("cutoff", "min_fdr", "n_sig")
  df = as.data.frame(df)
  
  p = ggplot(df, aes(x = cutoff, y = -log10(min_fdr) )) + 
    geom_col(color = "grey") + 
    ylab("Minimal FDR for GO Enrichment (-log10-transformed)") + 
    xlab("FC Cutoff") +
    geom_text(data = df, aes(x = cutoff, y = -log10(min_fdr) + 0.1, label = n_sig)) +
    scale_x_continuous(breaks = fc_cutoffs) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_noBG()
  
  plot(p)
  
  if (return_data){
    return(df)
  } else {
    return(p)
  }
}


rescue_fisher = function(mat1, mat2, genes){
  if (is.null(rownames(mat1)) || is.null(rownames(mat2))){
    stop("Error: both 'mat1' and 'mat2' must have row names as gene symbol/ID.")
  }
  if (!all(genes %in% rownames(mat1))){
    stop("Error: all 'genes' must be found in the row names 'mat1'.")
  }
  genes = unique(genes)
  df = data.frame(ref = sign(mat1[genes, 1]))
  df$query = 0
  rownames(df) = genes
  shared = intersect(genes, rownames(mat2))
  df[shared, "query"] = sign(mat2[shared, 1])
  df$res = df$ref * df$query
  df_up = df[which(df$ref > 0), ]
  df_dn = df[which(df$ref < 0), ]
  if (length(which(df$ref == 0))){
    stop("Error: some 'genes' in 'mat1' have zero fold changes in Column 1.")
  }
  #upreg - diff direction
  v1 = length(which(df_up$res == -1))
  #upreg - same direction
  v2 = length(which(df_up$res != -1))
  #downreg - diff direction
  v3 = length(which(df_dn$res == -1))
  #downreg - same direction
  v4 = length(which(df_dn$res != -1))
  tbt = data.frame("down_query" = c(v1, v4), "up_query" = c(v2, v3), row.names = c("up_ref", "down_ref"))
  pval = fisher.test(tbt)$p.value

  final = list(pval = pval, data = tbt, true_rescue_rate = (v1 + v3)/(v1 + v2 + v3 + v4))

  return(final)
}


trt_errorbar = function(baseline_ci_file, trt_ci_files, keyword){
  bl = read.table(baseline_ci_file, header = T, sep = "\t", check.names = F)
  bl = bl[which(bl[["name"]] == keyword), ]
  trt = c()
  for (v in trt_ci_files){
    tmp = read.table(v, header = T, sep = "\t", check.names = F)
    tmp = tmp[which(tmp[["name"]] == keyword), ]
    trt = rbind(trt, tmp)
  }
  trt = trt[order(trt[["LFC"]]), ]
  df = rbind(bl,trt)
  df = as.data.frame(df)
  #df$class = gsub("_vs_CON", "", df$class)
  df$class = factor(df$class, levels = as.vector(df[["class"]]))

  return(df)
}


find_paired_delta = function(l, col_paired, col_test, ref, target, col_val){
  tmp = split(l, l[[col_paired]])
  res = lapply(tmp, function(x){if ((nrow(x) == 2) & all(c(target, ref) %in% as.vector(x[[col_test]]))){return(x[which(x[[col_test]] == target), col_val] - x[which(x[[col_test]] == ref), col_val])}})
  
  return(unname(unlist(res)))
}


find_all_delta = function(df, col_test, ref, target, col_val, col_paired = NULL){
  res = c()
  if (is.null(col_paired)){
    ll = list(df)
  } else {
    df[[col_paired]] = factor(as.vector(df[[col_paired]]))
    ll = split(df, df[[col_paired]])
  }
  for (x in ll){
    tmp = split(x, x[[col_test]])
    if (length(tmp) == 2){
      val1 = tmp[[target]][[col_val]]
      val0 = tmp[[ref]][[col_val]]
      
      for (vv in val1){res = c(res, vv - val0)}
      
    } 
  }
  return(res)
}


trt_within_individual = function(df, col_nest, col_val, col_test, ref, target, col_paired = NULL, return_stats = T){
  subm = df[which(df[[col_test]] == ref | df[[col_test]] == target), ]
  subm[[col_test]] = factor(subm[[col_test]], levels = c(ref, target))
  tt = split(subm, subm[[col_nest]])
  dff = lapply(tt, FUN = find_all_delta, col_val = col_val, col_test = col_test, ref = ref, target = target, col_paired = col_paired)
  dff = unname(unlist(dff))
  dff = dff[is.finite(dff)]
  if (return_stats){
    u = mean(dff)
    if (is.null(col_paired)){
      n_target = nrow(subm[which(subm[[col_test]] == target), ])
      n_ref    = nrow(subm[which(subm[[col_test]] == ref), ])
      dof = max(n_target, n_ref)
    } else {
      dof = length(dff)
    }
    mgn = 1.96 * sd(dff)/sqrt(dof)
    
    return(data.frame(LFC = u, margin = mgn, n = dof, name = "X", class = target))
  } else {
    return(dff)
  }
}

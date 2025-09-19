library(stringr)
source("~/lib/R/geom_noBG.R")


enrichment_test <- function(term, go_db, query_genes, background_genes){
  gene_list <- toupper(go_db[[term]])
  
  term_gene <- length(gene_list)
  query <- length(query_genes)
  background <- length(background_genes)
  
  background_annotated_genes <- intersect(background_genes, gene_list)
  background_annotated <- length(background_annotated_genes)
  query_annotated_genes <- intersect(query_genes, gene_list)
  query_annotated_genes_print = paste(c(query_annotated_genes), collapse=', ' )
  query_annotated <- length(query_annotated_genes)
  
  not_query_annotated <- background_annotated  - query_annotated
  query_unannotated <-  query - query_annotated
  not_query_unannotated <- background - query - not_query_annotated
  
  expected <- query * background_annotated / background
  variance_term_1 <- (background - background_annotated) / background
  variance_term_2 <- (background - query) / (background - 1)
  stand_deviation <- sqrt(expected * variance_term_1 * variance_term_2)
  enrichment_zscore <- (query_annotated - expected) / stand_deviation
  tbt_table <- matrix(c(query_annotated, query_unannotated, not_query_annotated, not_query_unannotated), nrow = 2)
  fisher_result <- fisher.test(tbt_table, alternative = "greater")
  pvalue <- fisher_result$p.value
  odds_ratio <- unname(fisher_result$estimate)
  ci_lower <- fisher_result$conf.int[1]
  ci_upper <- fisher_result$conf.int[2]

  result_line <- data.frame(term, pvalue, odds_ratio, ci_lower, ci_upper, query_annotated, query_unannotated, not_query_annotated, not_query_unannotated, term_gene, enrichment_zscore, query_annotated_genes_print)
  
  return(result_line)
}


go_test <- function(target, background, db) {
  db_name = names(db)
  if (is.null(db_name)){db_name = paste0("GO_DB", 1:length(db))}
  
  res_final = c()
  for (i in 1:length(db)){
    ref_db = readRDS(db[[i]])
    message("GO database: ", db_name[i])
    
    db_genes = toupper(unique(unname(do.call(c, ref_db))))
    target <- toupper(target)
    background <- toupper(background)
    target_gene_found = which(target %in% db_genes)
    background_gene_found = which(background %in% db_genes)
    
    message(paste0("Found ",length(target_gene_found)," out of ",length(target)," (",round(100*length(target_gene_found)/length(target),2),"%) target genes in the GO database."))
    message(paste0("Found ",length(background_gene_found)," out of ",length(background)," (",round(100*length(background_gene_found)/length(background),2),"%) background genes in the GO database."))
    
    res_list = lapply(names(ref_db), enrichment_test, go_db = ref_db, query_genes = target[target_gene_found], background_genes = background[background_gene_found])
    res = do.call(rbind, res_list)
    
    res$fdr <- p.adjust(res$pvalue, method="BH")
    res$bonferroni <- p.adjust(res$pvalue, method="bonferroni")
    res$ontology <- db_name[i]
    
    res = res[,c(1, ncol(res), (ncol(res) - 1), (ncol(res) - 2), 2:(ncol(res) - 3))]
    
    res_final = rbind(res_final, res)
  }
  
  res_final = res_final[with(res_final, order(fdr, pvalue)), ]
  
  return(res_final)
}


go_extract_keywords = function(dat, col_term, keywords){
  idx_list = lapply(keywords, FUN = function(x){return(grep(x, as.vector(dat[[col_term]])))})
  idx = do.call(c, idx_list)
  idx = unique(idx)
  
  return(dat[idx, ])
}


go_plot = function(dat, col_val, cutoff, ontology = NULL, ntop = 5, show_hit = T, text_width = 20, text_size = 8){
  col_transformed = paste0(col_val, "_transformed")
  dat[["term"]] = gsub("_", " ", dat[["term"]])
  dat[["term"]] = factor(dat[["term"]], levels = rev(as.vector(dat[["term"]])))
  dat[[col_transformed]] = -log10(dat[[col_val]])
  dat[["plot_label"]] = as.vector(paste0("(", dat[["query_annotated"]], "/" , dat[["term_gene"]], ")"))
  if ((!is.null(ontology)) & (ontology %in% as.vector(dat[["ontology"]])) ){dat = dat[which(dat[["ontology"]] %in% ontology), ]}
  if ((!is.null(ntop)) | ntop > 0){dat = dat[1:ntop, ]}
  p = ggplot(dat, aes_string(x = "term", y = col_transformed, fill = "ontology")) + 
    geom_col() + 
    ylab("Significance (-log10-transformed)") + xlab("") + 
    coord_flip() + 
    geom_hline(yintercept = -log10(cutoff),linetype = "dashed", color = "black") + 
    geom_noBG()
  
  if (show_hit){
    p = p + geom_text(aes_string(x = "term", y = 0.05, label = "plot_label"), data = dat, hjust = "left")
  }
  
  p = p + scale_x_discrete(labels = function(x) str_wrap(x, width = text_width)) + theme(axis.text.y = element_text(size = text_size))
  plot(p)
  
  return(p)
}


go_change_plot = function(dat1, dat2, name1 = "ref", name2 = "query", col_val, col_filter, cutoff, show_hit = T, ntop = 20, text_width = 20, text_size = 8, y_label = "Significance (-log10-transformed)", transform_val = T){
  dat1[["term"]] = gsub("_", " ", dat1[["term"]])
  dat2[["term"]] = gsub("_", " ", dat2[["term"]])
  col_transformed = paste0(col_val, "_transformed")
  if (transform_val){
    dat1[[col_transformed]] = -log10(dat1[[col_val]])
    dat2[[col_transformed]] = -log10(dat2[[col_val]])
  } else {
    dat1[[col_transformed]] = dat1[[col_val]]
    dat2[[col_transformed]] = dat2[[col_val]]
  }
  
  shared = intersect(rownames(dat1), rownames(dat2))
  dat1 = dat1[shared, ][order(dat1[[col_filter]]), ]
  dat1 = dat1[which(dat1[[col_filter]] < cutoff), ]
  dat1 = dat1[1:min(ntop, nrow(dat1)), ]
  dat2 = dat2[rownames(dat1), ]
  dat1[["Set_name"]] = name1
  dat2[["Set_name"]] = name2
  dat = rbind(dat1, dat2)
  dat[["term"]] = factor(dat[["term"]], levels = rev(as.vector(dat1[["term"]])))
  dat[["plot_label"]] = as.vector(paste0("(", dat[["query_annotated"]], "/" , dat[["term_gene"]], ")"))
  p = ggplot(dat, aes_string(x = "term", y = col_transformed, fill = "Set_name", label = "plot_label")) + 
    geom_col(position = "dodge") + 
    ylab(y_label) + xlab("") + 
    coord_flip() + 
    geom_noBG()
  
  if (transform_val){
    p = p + geom_hline(yintercept = -log10(cutoff),linetype = "dashed", color = "black")
  }
  
  if (show_hit){
    p = p + geom_text(aes(y = 0.05), hjust = "left", position = position_dodge(width = 0.9), size = text_size/4)
  }
  
  p = p + scale_x_discrete(labels = function(x) str_wrap(x, width = text_width)) + theme(axis.text.y = element_text(size = text_size))
  plot(p)
  
  return(p)
}

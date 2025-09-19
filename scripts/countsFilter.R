createConditionList = function(meta, conditions){
  meta_new = meta[, conditions]
  new_condition = apply(meta_new, 1, paste, collapse = "_")
  conditionList = lapply(unique(new_condition), function(x){return(unname(which(new_condition == x)))})
  return(conditionList)
}


countsFilter = function(mat, conditionList, cut = 0.1, pct = 0.5){
  condition_valid_pct = lapply(conditionList, function(x){length(which( mat[x] >= cut )) / length(mat[x]) })
  if (max(unlist(condition_valid_pct)) >= pct){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


psiFilter = function(mat, conditionList, min_valid_sample = 3){
  condition_valid = lapply(conditionList, function(x){length(which(   !is.na(mat[x])    )) })
  if (min(unlist(condition_valid)) >= min_valid_sample){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


getTotalReads = function(pathList){
  sizeVec=c()
  sampleVec = sapply(strsplit(pathList,"\\."),"[[",2)
  for (i in 1:length(pathList)){
    temp = scan(paste0(pathList[i],"/",sampleVec[i],".Log.final.out"),what="character",sep="\n",quiet=TRUE)
    for (x in temp){
      if (grepl("Uniquely mapped reads number", x)){
        sizeVec = c(sizeVec,as.numeric(strsplit(x,"\t")[[1]][2]))
      }
    }
  } 
  names(sizeVec) = sampleVec
  return(sizeVec)
}

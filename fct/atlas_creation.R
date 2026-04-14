
atlas_creation <- function(X, PATH_MetaboAnalyst, output_dir = dir()){
  
  metaboAnalyst_results <- read.table(PATH_MetaboAnalyst, sep = ',', header = T)
  
  ## PATHWAYS
  rownames(metaboAnalyst_results) = metaboAnalyst_results$Query
  metaboAnalyst_results = metaboAnalyst_results[colnames(X),]
  
  # a/ rename
  Q = metaboAnalyst_results$Query ; K = metaboAnalyst_results$KEGG
  miss_matched = Q[is.na(K)]
  
  K[is.na(K)] = seq(1, length(which(is.na(K))), by = 1)
  
  # b/ get pathways from KEGG
  mat = matrix(NA, nrow = length(K), ncol = 50)
  
  for(i in 1:length(K)){
    
    if(is.na(as.numeric(K[i]))==TRUE){
      pw = KEGGREST::keggLink("pathway", K[i])
      
      if(length(pw)>0){
        res = KEGGREST::keggGet(KEGGREST::keggLink("pathway", K[i]))
        N = unlist(lapply(res, function(x) x$NAME))
        mat[i,1:(5+length(N))] = c(as.character(metaboAnalyst_results[i,1:5]), N)
        
      }else{mat[i,1:5] = as.character(metaboAnalyst_results[i,1:5])}
    }else{mat[i,1:5] = as.character(metaboAnalyst_results[i,1:5])}
  }
  
  colnames(mat) = c(colnames(metaboAnalyst_results)[1:5], paste("pathways", 1:45))
  write.csv(mat, file = paste0(dir(),"/atlas.csv"))
  
}

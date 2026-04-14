library(KEGGREST)
library(tidyverse)
library(KODAMA)

hierarchical_creation <- function(X, normalisation = TRUE, PATH_MetaboAnalyst, cutoff = 5){
  
  ## NORMALISATION
  if(normalisation == TRUE){
    u <- X
    result <- KODAMA::normalization(u, method = "pqn")
    data <- result$newXtrain
  }else{data = X}
  
  ## B/ PATHWAYS
  metaboAnalyst_results <- read.table(PATH_MetaboAnalyst, sep = ',', header = T)
  
  rownames(metaboAnalyst_results) = metaboAnalyst_results$Query
  metaboAnalyst_results = metaboAnalyst_results[colnames(X),]
  
  # a/ rename
  Q = metaboAnalyst_results$Query ; K = metaboAnalyst_results$KEGG
  miss_matched = Q[is.na(K)]
  
  # b/ get pathways from KEGG
  pathways <- keggLink("pathway", K)
  
  pw = data.frame(KEGG_ID = sub("cpd:", "", names(pathways)), 
                  PTW = sub("path:", "",as.character(pathways)))
  
  K[is.na(K)] = seq(1, length(which(is.na(K))), by = 1)
  K[which(duplicated(K))] = 
    paste0(K[which(duplicated(K))], "_1")
  
  rownames(metaboAnalyst_results) = K
  
  MB_pathways = list()
  class_metabolites = NULL
  
  for(i in unique(pw$PTW)){
    # Metabolites inside the Block
    metabolites_inside = t(data[, metaboAnalyst_results[pw$KEGG_ID[which(pw$PTW==i)],]$Query])
    metabolites_inside = metabolites_inside[duplicated(rownames(metabolites_inside))==FALSE,]
    
    if(nrow(metabolites_inside) > cutoff){ # We select only pathways with more than 5 metabolites.
      # Metadata of the pathway
      CLASS = keggGet(i)
      if(is.character(CLASS[[1]]$CLASS)){class_metabolites = c(class_metabolites, CLASS[[1]]$CLASS)
      }else{class_metabolites = c(class_metabolites, NA)}
      
      # Construction of the list
      if(length(CLASS[[1]]$NAME)>1){
        MB_pathways[[CLASS[[1]]$NAME$NAME]] = metabolites_inside
      }else{MB_pathways[[CLASS[[1]]$NAME]] = metabolites_inside}
    }
  }
  
  return(list("X_list" = MB_pathways, 
              "Functions of MB" = class_metabolites,
              "Metabolties_unrecognised" = miss_matched))
}








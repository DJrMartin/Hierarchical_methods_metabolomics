library(FactoMineR)

rda_manual <- function(Y, X, scale.Y = TRUE, scale.X = TRUE, col_sites, pch_sites, graph = TRUE, nperm = 999,
                       arrows_scale = 2, sel.var = 0.2) {
  
  # Y = Y[, groups == "O2/CO2"]
  # X = X$`Metabolic pathways`
  # pch_sites = pch.CMS
  # col_sites = as.character(col.Y)
  # graph = F
  
  ## --- 1. Préparation des matrices ---
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  
  w_y <- which(rowSums(is.na(Y)) == 0 & rowSums(is.na(X)) == 0)
  
  if(scale.Y){Y <- scale(Y[w_y,])}else{Y <- Y[w_y,]}
  if(scale.X){X <- scale(X[w_y,])}else{X <- X[w_y,]}
  
  col_sites = col_sites[w_y]
  pch_sites = pch_sites[w_y]
  
  n <- nrow(X)
  p <- ncol(X)
  
  if(p > 20){
    w = which(rowMeans(abs(cor(X, Y)))>0.3)
    if(length(w)>1){X = X[,w]}else{X = X}
  }
  
  ## --- 2. PCA sur Y ---
  pca_Y <- PCA(Y, graph = F)
  
  Y_pca <- pca_Y$ind$coord                    # scores individus
  eig_Y <- pca_Y$eig[,1]                      # valeurs propres PCA
  
  ## --- 3. Projection RDA (régression multiple) ---
  Y_hat <- apply(Y_pca, 2, function(z)
    lm(z ~ X)$fitted.values
  )
  
  Y_res <- Y_pca - Y_hat
  
  ## --- 4. PCA sur la partie contrainte ---
  pca_RDA <- PCA(Y_hat, graph = F)
  eig_RDA <- pca_RDA$eig[,1]                      
  
  ## --- 5. Partition de variance ---
  var_total <- sum(apply(Y, 2, var))
  var_constrained <- sum(apply(Y_hat, 2, var))
  var_residual <- var_total - var_constrained
  
  ## --- 5. Test de permutation (global RDA) ---
  R2_perm <- numeric(nperm)
  for (i in seq_len(nperm)) {
    
    # 1. Permutation des lignes
    Y_perm <- Y[sample(1:n), ]
    
    # 2. PCA sur Y permuté
    pca_perm <- FactoMineR::PCA(Y_perm, graph = FALSE)
    Y_pca_perm <- pca_perm$ind$coord
    
    # 3. Régression (RDA)
    Y_hat_perm <- apply(Y_pca_perm, 2, function(z)
      lm(z ~ X)$fitted.values
    )
    
    # 4. Variance expliquée
    var_constrained_perm <- sum(apply(Y_hat_perm, 2, var))
    
    # 5. R²
    R2_perm[i] <- var_constrained_perm / var_total
  }
  
  # R² ajusté observé
  R2_obs <- var_constrained / var_total
  R2_adj <- 1 - (1 - R2_obs) * (n - 1) / (n - ncol(X) - 1)
  
  # p-value
  p_value <- (sum(R2_perm >= R2_obs) + 1) / (nperm + 1)
  
  partition <- c(
    constrained = var_constrained / var_total,
    residual    = var_residual    / var_total,
    p.value = p_value,
    R2_adj = R2_adj
  )
  
  ## --- 6. Scores des sites (scaling 1) ---
  sites <- sweep(
    pca_RDA$ind$coord,
    2,
    sqrt(eig_RDA),
    "/"
  )
  
  ## --- 7. Flèches Y (variables réponses) ---
  var_Y <- cor(Y, pca_RDA$ind$coord, use = "complete.obs")
  var_Y <- sweep(var_Y, 2, sqrt(eig_RDA), "/")
  
  ## --- 8. Flèches X (variables explicatives) ---
  var_X <- cor(X, pca_RDA$ind$coord, use = "complete.obs")
  var_X <- sweep(var_X, 2, sqrt(eig_RDA), "/")
  
  select.var.Y = var_Y[which(rowMeans(abs(var_Y[,1:2]))>sel.var),]
  select.var.X = var_X[which(rowMeans(abs(var_X[,1:2]))>sel.var),]
  
  if(graph == TRUE){
    plot(
      scale(sites),
      xlab = paste0("RDA 1 (", round(eig_RDA[1]/sum(eig_RDA), 2)*100, "%)"),
      ylab = paste0("RDA 2 (", round(eig_RDA[2]/sum(eig_RDA), 2)*100, "%)"),
      asp = 1,
      pch = pch_sites,
      col = col_sites
    )
    abline(h = 0, v = 0, lty = 2, col = "grey")

    ## Flèches Y
    arrows(
      0, 0,
      select.var.Y[, 1]*arrows_scale, select.var.Y[, 2]*arrows_scale,
      length = 0.08,
      col = "darkgreen"
    )
    text(
      select.var.Y[, 1]*arrows_scale, select.var.Y[, 2]*arrows_scale,
      labels = rownames(select.var.Y),
      col = "darkgreen",
      cex = 0.7,
      pos = 3
    )
    
    ## Flèches X
    arrows(
      0, 0,
      select.var.X[, 1]*arrows_scale, select.var.X[, 2]*arrows_scale,
      length = 0.08,
      col = "cornflowerblue"
    )
    text(
      select.var.X[, 1]*arrows_scale, select.var.X[, 2]*arrows_scale,
      labels = rownames(select.var.X),
      col = "cornflowerblue",
      cex = 0.8,
      pos = 4
    )
  }
  
  ## --- 9. Sortie ---
  return(list(
    eig_RDA   = eig_RDA,
    partition = data.frame(partition),
    sites     = sites,
    var_Y     = var_Y,
    var_X     = var_X
  ))

}

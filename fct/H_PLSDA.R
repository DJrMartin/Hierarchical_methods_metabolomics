library(mixOmics)

H_PLSDA <- function(X_TRAIN, Y_TRAIN, new.data = NULL, new.Y = NULL, 
                    ncomp_lower.level = 2, ncomp_upper.level = 5){
  
  if(is.null(new.data)==TRUE){new.data = X_TRAIN ; new.Y = Y_TRAIN}

  #################
  ## LOWER LAYER ##
  #################
  res_pls = list()
  LOADINGS_variables = list()
  SCORES_train_list = list()
  SCORES_test_list = list()
  
  for(j in 1:length(X_TRAIN)){
    
    block_name <- names(X_TRAIN)[j]
    
    # PLS-DA
    res_pls[[block_name]] <- mixOmics::plsda(
      X = X_TRAIN[[j]], 
      Y = Y_TRAIN, 
      ncomp = ncomp_lower.level
    )
    
    LOADINGS_variables[[block_name]] <- res_pls[[block_name]]$loadings$X
    
    ## TRAIN SCORES
    scores_train <- res_pls[[block_name]]$variates$X
    
    ## HMFA normalization (by first eigenvalue)
    lambda1 <- var(scores_train[,1])
    scores_train <- scores_train / sqrt(lambda1)
    
    SCORES_train_list[[block_name]] <- scores_train
    
    ## TEST SCORES
    scores_test <- predict(res_pls[[block_name]], newdata = new.data[[j]])$variates
    
    scores_test <- scores_test / sqrt(lambda1)
    
    SCORES_test_list[[block_name]] <- scores_test
  }
  
  ## CONCATENATION
  SCORES_train <- do.call(cbind, SCORES_train_list)
  SCORES_test  <- do.call(cbind, SCORES_test_list)
  
  colnames(SCORES_train) <- colnames(SCORES_test) <-
    paste(rep(names(X_TRAIN), each = ncomp_lower.level), 1:ncomp_lower.level, sep = "_")
  
  #################
  ## UPPER LAYER ##
  #################
  ## TRAINING model. Input : PLS scores of the train set from the lower level PLS models
  ### MODEL
  res.plsda <- mixOmics::plsda(X = SCORES_train, Y = Y_TRAIN, ncomp = ncomp_upper.level)
  ### SCORES
  scores_train_set <- data.frame(res.plsda$variates$X)
  ### EIGENVALUES
  Eigenvalue = res.plsda$prop_expl_var$X
  ### PIP
  LOADINGS_pathways = res.plsda$loadings$X
  
  ## PREDICTION. Input : PLS scores of the test set.
  ### SCORES
  scores_test_set <- predict(res.plsda, SCORES_test)$variates
  ### MODEL ACCURACY
  Accuracy <- caret::confusionMatrix(factor(predict(res.plsda, SCORES_test)$class$max.dist[,ncomp_upper.level], levels = levels(new.Y)), new.Y)
  
  ####################
  ## MODEL ACCURACY ##
  ####################
  ## Q2
  Y_test_dummy <- model.matrix(~ new.Y - 1)
  Y_test_pred <- predict(res.plsda, SCORES_test)$predict[, , ncomp_upper.level]
  
  PRESS <- sum((Y_test_dummy - Y_test_pred)^2)
  TSS <- sum((Y_test_dummy - mean(Y_test_dummy))^2)
  
  Q2 <- 1 - PRESS / TSS
  
  ## R2
  Yb_dummy <- model.matrix(~ Y_TRAIN - 1)
  Yb_pred <- predict(res.plsda, SCORES_train)$predict[, , ncomp_upper.level]
  
  SS_res <- sum((Yb_dummy - Yb_pred)^2)
  SS_tot <- sum((Yb_dummy - mean(Yb_dummy))^2)
  
  R2 <- 1 - SS_res / SS_tot
  
  return(list(Upper_loadings = LOADINGS_pathways,
              Lower_loadings = LOADINGS_variables,
              Scores = list("Train" = scores_train_set, "Test" = scores_test_set),
              Accuracy = Accuracy,
              Q2 = Q2, 
              R2 = R2,
              Eigenvalue = Eigenvalue))
}

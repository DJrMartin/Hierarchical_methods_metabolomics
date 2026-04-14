simulate_multiblock_data <- function(
    n = 200,
    p = c(80, 5, 100),
    rho = c(0.3, 0.5, 0.8),
    signal_strength = list(
      block1 = c(0, 0.4),   # weak signal
      block2 = c(0.2, 0.6)  # moderate signal
    ),
    noise_sd = 3,
    train_prop = 0.60,
    seed = NULL
){
  
  if(!is.null(seed)) set.seed(seed)
  
  library(MASS)
  library(caret)
  
  ## =========================
  ## 1. Generate blocks
  ## =========================
  
  X_list <- list()
  Sigma_list <- list()
  
  for(i in seq_along(p)){
    Sigma_list[[i]] <- rho[i]^abs(outer(1:p[i], 1:p[i], "-"))
    X_list[[i]] <- MASS::mvrnorm(n, mu = rep(0, p[i]), Sigma = Sigma_list[[i]])
  }
  
  names(X_list) <- paste0("block", seq_along(p))
  
  ## =========================
  ## 2. Generate signal
  ## =========================
  
  beta_list <- list()
  
  # Block 1 = weak signal
  beta_list[[1]] <- runif(p[1], signal_strength$block1[1], signal_strength$block1[2])
  
  # Block 2 = moderate signal
  beta_list[[2]] <- runif(p[2], signal_strength$block2[1], signal_strength$block2[2])
  
  # Block 3 = pure noise (no signal)
  beta_list[[3]] <- rep(0, p[3])
  
  signal_X1 <- X_list[[1]] %*% beta_list[[1]]
  signal_X2 <- X_list[[2]] %*% beta_list[[2]]
  
  epsilon <- rnorm(n, sd = noise_sd)
  
  ## =========================
  ## 3. Response variable
  ## =========================
  
  Y_cont <- signal_X2 + 0.5 * (X_list[[2]][,1]^2) + signal_X1 + epsilon
  
  Y <- ifelse(Y_cont > median(Y_cont), 1, 0)
  Y <- factor(Y)
  
  ## =========================
  ## 4. Normalization
  ## =========================
  
  X_list <- lapply(X_list, scale)
  MB <- lapply(X_list, data.frame)
  
  ## =========================
  ## 5. Train/Test split
  ## =========================
  
  intraining <- caret::createDataPartition(Y, p = train_prop, list = FALSE)
  
  train <- lapply(MB, function(x) x[intraining, ])
  test  <- lapply(MB, function(x) x[-intraining, ])
  
  Y_train <- Y[intraining]
  Y_test  <- Y[-intraining]
  
  ## =========================
  ## Output
  ## =========================
  
  return(list(
    MB = MB,
    train = train,
    test = test,
    Y = Y,
    Y_train = Y_train,
    Y_test = Y_test,
    beta = beta_list
  ))
}

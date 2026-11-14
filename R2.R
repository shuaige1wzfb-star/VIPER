library(randomForest)
reparameterize <- function(X, s) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  L <- as.numeric(X[, s] < 1e-9)
  
  W <- X
  W[, s] <- 0
  
  row_sums <- rowSums(W)
  row_sums[row_sums < 1e-9] <- 1
  W <- W / row_sums
  
  return(list(
    W = W,           # n×p 
    L = L           # n×1 
  ))
}

# W -> Y
select_reg_model <- function(W_train, Y_train, W_test) {
  rf_shallow <- randomForest(
    x = W_train,
    y = as.vector(Y_train),
    ntree = 250
  )
  
  return(list(pred = predict(rf_shallow, W_test), type = "shallow_rf"))
}

# W -> L 
select_class_model <- function(W_train, L_train, W_test) {
  rf_shallow <- randomForest(
    x = W_train,
    y = factor(L_train, levels = c(0, 1)),
    ntree = 250
  )
  probs <- predict(rf_shallow, W_test, type = "prob")
  
  return(list(
    pred = probs[, "1"],
    type = "shallow_rf",
    probs = probs
  ))
}

select_class_model2 <- function(W_train, L_train, W_test) {
  m_hat <- mean(L_train)
  m_pred <- rep(m_hat, nrow(W_test))

  list(
    pred = m_pred,
    type = "intercept_only",
    probs = cbind(`0` = 1 - m_pred, `1` = m_pred)
  )
}

# algorithm
R2_partial_linear <- function(Y, L, W, s, nfolds = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  W <- as.matrix(W)
  n <- nrow(W)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  u_hat <- v_hat <- numeric(n)
  
  for (k in seq_len(nfolds)) {
    test <- fold_id == k
    train <- !test
    W_train <- W[train, , drop = FALSE]
    W_test <- W[test, , drop = FALSE]
    
    g_pred <- select_reg_model(W_train, Y[train], W_test)$pred
    m_pred <- select_class_model(W_train, L[train], W_test)$pred
    
    u_hat[test] <- Y[test] - g_pred
    v_hat[test] <- L[test] - m_pred
  }
  
  J_hat <- mean(v_hat^2)
  theta_hat <- mean(u_hat * v_hat) / J_hat
  
  q_i <- (1 - L)^2
  q_hat <- mean(q_i)
  Y_centered_sq <- (Y - mean(Y))^2
  sigma_Y_sq <- mean(Y_centered_sq)
  phi_hat <- theta_hat^2 * q_hat / sigma_Y_sq
  
  varphi_theta <- v_hat * (u_hat - theta_hat * v_hat) / J_hat
  varphi <- (
    2 * theta_hat * q_hat * varphi_theta +
      theta_hat^2 * (q_i - q_hat) -
      phi_hat * (Y_centered_sq - sigma_Y_sq)
  ) / sigma_Y_sq
  
  list(
    phi_hat = phi_hat,
    sigma_phi_sq = mean(varphi^2),
    theta_hat = theta_hat,
    varphi = varphi,
    varphi_theta = varphi_theta,
    n = n,
    p = ncol(W),
    s = s
  )
}

# algorithm with intercept_only
R2_partial_linear2 <- function(Y, L, W, s, nfolds = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  W <- as.matrix(W)
  n <- nrow(W)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  u_hat <- v_hat <- numeric(n)
  
  for (k in seq_len(nfolds)) {
    test <- fold_id == k
    train <- !test
    W_train <- W[train, , drop = FALSE]
    W_test <- W[test, , drop = FALSE]
    
    g_pred <- select_reg_model(W_train, Y[train], W_test)$pred
    m_pred <- select_class_model2(W_train, L[train], W_test)$pred
    
    u_hat[test] <- Y[test] - g_pred
    v_hat[test] <- L[test] - m_pred
  }
  
  J_hat <- mean(v_hat^2)
  theta_hat <- mean(u_hat * v_hat) / J_hat
  
  q_i <- (1 - L)^2
  q_hat <- mean(q_i)
  Y_centered_sq <- (Y - mean(Y))^2
  sigma_Y_sq <- mean(Y_centered_sq)
  phi_hat <- theta_hat^2 * q_hat / sigma_Y_sq
  
  varphi_theta <- v_hat * (u_hat - theta_hat * v_hat) / J_hat
  varphi <- (
    2 * theta_hat * q_hat * varphi_theta +
      theta_hat^2 * (q_i - q_hat) -
      phi_hat * (Y_centered_sq - sigma_Y_sq)
  ) / sigma_Y_sq
  
  list(
    phi_hat = phi_hat,
    sigma_phi_sq = mean(varphi^2),
    theta_hat = theta_hat,
    varphi = varphi,
    varphi_theta = varphi_theta,
    n = n,
    p = ncol(W),
    s = s
  )
}

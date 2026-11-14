# ============================================================
# VIPER positive-threshold test under correct and misspecified models
#
# Procedures:
#      H0: phi_j <= c versus H1: phi_j > c
#   3. Log-contrast plug-in with full-refit pairs bootstrap
#   4. Random-forest plug-in with full-refit pairs bootstrap
#
# VIPER uses its influence-function standard error. The two plug-in methods
# use full-refit pairs bootstrap intervals. Each threshold test rejects when
# its one-sided lower confidence bound exceeds c.
#
# Data:
#   Correlated logistic-normal compositions
#   with controlled structural zeros
# ============================================================

rm(list = ls())

library(MASS)
library(randomForest)
library(doParallel)
library(foreach)

set.seed(2026)


# ============================================================
# 0. User settings
# ============================================================

# Run from a working directory containing R2.R.
algorithm_file <- normalizePath("R2.R", mustWork = TRUE)
source(algorithm_file)

n_iter <- 300        
n <- 100             
p <- 50            
nfolds <- 3
alpha <- 0.05
rf_ntree <- 300
n_boot <- 200
n_oracle <- 100000

# Choose "correct_logcontrast", "nonlinear_plm", or "both":
#   "correct_logcontrast": the log-contrast model is correctly specified.
#   "nonlinear_plm": the log-contrast model is misspecified, while the
#                    partially linear model for VIPER remains correct.
scenario_to_run <- "both"

# Ridge penalties used inside the log-contrast plug-in learner
ridge_lambda_grid <- 10^seq(-4, 8, by = 1)

sig <- 0.7
pseud <- 1e-10

# Positive VIPER threshold
c_threshold <- 0.01

# Fix one target component throughout the power curve.
target_j <- 1

# Used only under the correctly specified log-contrast scenario.
reference_j <- 7

# Fix the target component presence probability across scenarios.
q_value <- 0.2

# Dense candidate grids are used only to locate effects whose oracle phi lies
# below, near, and above c. The simulation subsequently uses the selected grid.
effect_candidates_correct <- seq(0, 0.15, length.out = 30)
effect_candidates_nonlinear <- seq(0, 1.50, length.out = 30)

# Target oracle-importance levels. These guarantee points below, near, and
# above the threshold whenever the candidate grid spans the requested range.
oracle_level_multipliers <- c(0, 0.3, 0.6, 1, 3, 5, 8, 10)

if (!scenario_to_run %in% c("correct_logcontrast", "nonlinear_plm", "both")) {
  stop("scenario_to_run must be 'correct_logcontrast', 'nonlinear_plm', or 'both'.")
}

scenario_names <- if (scenario_to_run == "both") {
  c("correct_logcontrast", "nonlinear_plm")
} else {
  scenario_to_run
}

J_signal <- target_j

base_seed <- 2026

# ============================================================
# 1. Multivariate normal settings
# ============================================================
rho <- 0.5
mu <- rep(0, p) 
Sigma <- toeplitz(rho^(0:(p - 1)))

# ============================================================
# 2. Logistic-normal composition with structural zeros
# ============================================================

generate_logistic_normal <- function(n, p, s, rho, presence_prob) {
  mu <- rep(0, p) 
  Sigma <- toeplitz(rho^(0:(p - 1)))
  
  W <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  X_raw <- exp(W)
  X_raw <- X_raw / rowSums(X_raw)
  
  present <- matrix(
    rbinom(n * p, size = 1, prob = presence_prob),
    nrow = n,
    ncol = p
  )
  
  X <- X_raw * present
  
  zero_rows <- which(rowSums(X) <= 1e-14)
  
  if (length(zero_rows) > 0L) {
    for (i in zero_rows) {
      keep <- sample(seq_len(p), size = 2L)
      X[i, keep] <- X_raw[i, keep]
    }
  }
  
  X <- X / rowSums(X)
  return(X)
}

make_beta <- function(effect_value, p, target_j, reference_j) {
  beta <- rep(0, p)
  beta[target_j] <- effect_value
  beta[reference_j] <- -effect_value
  
  return(beta)
}

get_true_regression <- function(X, target_j, effect_value, scenario, reference_j, pseud) {
  rep_X <- reparameterize(X = X, s = target_j)
  L <- rep_X$L
  W <- rep_X$W
  
  if (scenario == "correct_logcontrast") {
    beta <- make_beta(
      effect_value = effect_value,
      p = ncol(X),
      target_j = target_j,
      reference_j = reference_j
    )
    
    # This is the conditional mean used to generate the response.
    conditional_mean <- as.vector(log(X + pseud) %*% beta)
    
    return(list(
      conditional_mean = conditional_mean,
      beta = beta,
      L = L,
      W = W,
      h = NA_real_
    ))
  }
  
  # ----------------------------------------------------------
  # Nonlinear PLM:
  #   Y = theta_j L_j + h_j(W_j) + epsilon.
  # The log-contrast model is misspecified, but the VIPER PLM
  # is correctly specified by construction.
  # ----------------------------------------------------------
  
  available_columns <- setdiff(seq_len(ncol(W)), target_j)
  
  if (length(available_columns) < 6L) {
    stop("At least six non-target components are required.")
  }
  
  selected_columns <- available_columns[1:6]
  ratio_offset <- 1e-10
  
  balance_1 <- log(
    (W[, selected_columns[1]] + ratio_offset) /
      (W[, selected_columns[2]] + ratio_offset)
  )
  
  balance_2 <- log(
    (W[, selected_columns[3]] + ratio_offset) /
      (W[, selected_columns[4]] + ratio_offset)
  )
  
  balance_3 <- log(
    (W[, selected_columns[5]] + ratio_offset) /
      (W[, selected_columns[6]] + ratio_offset)
  )
  
  h <- 0.8 * sin(balance_1) + 0.5 * tanh(balance_2)^2 + 0.4 * tanh(balance_1) * tanh(balance_3)
  
  conditional_mean_Z <- effect_value * L + h
  
  # When L=0 the perturbation is effective; replacing Z by W changes
  # the conditional mean by theta_j. When L=1, Z already equals W.
  perturbation_contrast <- effect_value * (1 - L)
  
  return(list(
    conditional_mean = conditional_mean_Z,
    beta = NULL,
    L = L,
    W = W,
    h = h,
    perturbation_contrast = perturbation_contrast
  ))
}

generate_response <- function(X, target_j, effect_value, scenario, reference_j, sig, pseud) {
  truth <- get_true_regression(
    X = X,
    target_j = target_j,
    effect_value = effect_value,
    scenario = scenario,
    reference_j = reference_j,
    pseud = pseud
  )
  
  truth$y <- truth$conditional_mean + rnorm(nrow(X), sd = sig)
  truth
}

build_oracle_design <- function(scenario, q_value, effect_candidates, n_oracle, seed) {
  set.seed(seed)
  
  X_oracle <- generate_logistic_normal(
    n = n_oracle,
    p = p,
    s = s,
    rho = rho,
    presence_prob = q_value
  )
  oracle_error <- rnorm(n_oracle, mean = 0, sd = sig)
  
  oracle_rows <- lapply(seq_along(effect_candidates), function(effect_id) {
    effect_value <- effect_candidates[effect_id]
    
    cat(sprintf(
      "Oracle: scenario=%s, q=%.2f, effect=%d/%d (%.4f)\n",
      scenario,
      q_value,
      effect_id,
      length(effect_candidates),
      effect_value
    ))
    flush.console()
    
    truth_now <- get_true_regression(
      X = X_oracle,
      target_j = target_j,
      effect_value = effect_value,
      scenario = scenario,
      reference_j = reference_j,
      pseud = pseud
    )
    y_oracle <- truth_now$conditional_mean + oracle_error
    sigma_y_sq <- var(y_oracle)
    
    if (scenario == "correct_logcontrast") {
      # Evaluate the same known population regression function at Z and W.
      f_Z <- truth_now$conditional_mean
      f_W <- as.numeric(log(truth_now$W + pseud) %*% truth_now$beta)
      phi_value <- mean((f_W - f_Z)^2) / sigma_y_sq
      theta_value <- NA_real_
    } else {
      theta_value <- effect_value
      phi_value <- theta_value^2 * q_value / sigma_y_sq
    }
    
    data.frame(
      scenario = scenario,
      q = q_value,
      effect_value = effect_value,
      oracle_theta = theta_value,
      oracle_phi = phi_value,
      oracle_sigma_y_sq = sigma_y_sq
    )
  })
  
  do.call(rbind, oracle_rows)
}

select_effect_grid <- function(oracle_table, threshold, multipliers) {
  target_levels <- threshold * multipliers
  
  selected_indices <- vapply(target_levels, function(target_level) {
    which.min(abs(oracle_table$oracle_phi - target_level))
  }, integer(1))
  
  selected <- oracle_table[unique(selected_indices), ]
  selected <- selected[order(selected$effect_value), ]
  
  if (!any(selected$oracle_phi < threshold) ||
      !any(selected$oracle_phi >= threshold) ||
      max(selected$oracle_phi) <= threshold) {
    stop(
      "The candidate effect grid does not adequately cover phi<c, phi=c, and phi>c."
    )
  }
  
  selected$threshold_region <- ifelse(
    selected$oracle_phi < 0.90 * threshold,
    "Below",
    ifelse(
      selected$oracle_phi > 1.10 * threshold,
      "Above",
      "Boundary"
    )
  )
  
  selected
}

interpolate_threshold_effect <- function(oracle_table, threshold) {
  ordered <- oracle_table[order(oracle_table$oracle_phi), ]
  
  if (threshold < min(ordered$oracle_phi) ||
      threshold > max(ordered$oracle_phi)) {
    return(NA_real_)
  }
  
  approx(
    x = ordered$oracle_phi,
    y = ordered$effect_value,
    xout = threshold,
    ties = "ordered"
  )$y
}

# ============================================================
# 4. Prediction learners for perturbation plug-in importance
#
# Log-contrast plug-in uses constrained ridge regression: sum(beta)=0.
# ============================================================
fit_constrained_ridge <- function(y, x, lambda) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  
  x_mean <- colMeans(x)
  y_mean <- mean(y)
  
  xc <- sweep(x, 2, x_mean, FUN = "-")
  yc <- y - y_mean
  
  p_now <- ncol(xc)
  constraint <- matrix(1, nrow = p_now, ncol = 1)
  
  A <- crossprod(xc) + lambda * diag(p_now)
  rhs <- crossprod(xc, yc)
  
  A_inv_rhs <- solve(A, rhs)
  A_inv_constraint <- solve(A, constraint)
  
  correction <- A_inv_constraint %*%
    solve(t(constraint) %*% A_inv_constraint) %*%
    t(constraint) %*% A_inv_rhs
  
  beta_hat <- as.numeric(A_inv_rhs - correction)
  intercept_hat <- y_mean - sum(x_mean * beta_hat)
  
  return(list(
    intercept = intercept_hat,
    beta = beta_hat,
    lambda = lambda
  ))
}

predict_constrained_ridge <- function(fit, x_new) {
  as.numeric(fit$intercept + as.matrix(x_new) %*% fit$beta)
}

fit_logcontrast_cv <- function(y, x, lambda_grid, nfolds = 3, seed = 1) {
  n_now <- length(y)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n_now))
  
  cv_error <- rep(0, length(lambda_grid))
  
  for (lambda_id in seq_along(lambda_grid)) {
    lambda_now <- lambda_grid[lambda_id]
    fold_error <- rep(NA_real_, nfolds)
    
    for (fold_now in seq_len(nfolds)) {
      test_id <- which(fold_id == fold_now)
      train_id <- setdiff(seq_len(n_now), test_id)
      
      fit_now <- fit_constrained_ridge(
        y = y[train_id],
        x = x[train_id, , drop = FALSE],
        lambda = lambda_now
      )
      
      pred_now <- predict_constrained_ridge(
        fit_now,
        x[test_id, , drop = FALSE]
      )
      
      fold_error[fold_now] <- mean((y[test_id] - pred_now)^2)
    }
    
    cv_error[lambda_id] <- mean(fold_error)
  }
  
  selected_lambda <- lambda_grid[which.min(cv_error)]
  
  fit_final <- fit_constrained_ridge(
    y = y,
    x = x,
    lambda = selected_lambda
  )
  
  fit_final$cv_error <- cv_error
  return(fit_final)
}

# Cross-fitted original and perturbed predictions for plug-in importance.
compute_plugin_predictions <- function(X, y, J_signal, nfolds, seed, pseud, ridge_lambda_grid, rf_ntree) {
  n_now <- nrow(X)
  m <- length(J_signal)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n_now))
  
  full_lc_pred <- rep(NA_real_, n_now)
  full_rf_pred <- rep(NA_real_, n_now)
  
  plugin_lc_perturbed_pred <- matrix(NA_real_, nrow = n_now, ncol = m)
  plugin_rf_perturbed_pred <- matrix(NA_real_, nrow = n_now, ncol = m)
  
  # W_full_list retains all p coordinates. For each target component, the
  # regression fitted on Z is evaluated at both Z and W.
  W_full_list <- lapply(
    J_signal,
    function(j) {
      W <- X
      W[, j] <- 0
      remaining_mass <- rowSums(W)
      
      if (any(remaining_mass <= 1e-14)) {
        stop("Perturbation produced a row with no remaining mass.")
      }
      
      W / remaining_mass
    }
  )
  
  log_X <- log(X + pseud)
  log_W_full_list <- lapply(W_full_list, function(W) log(W + pseud))
  
  for (fold_now in seq_len(nfolds)) {
    test_id <- which(fold_id == fold_now)
    train_id <- setdiff(seq_len(n_now), test_id)
    
    # --------------------------------------------------------
    # Full log-contrast predictor
    # --------------------------------------------------------
    full_lc_fit <- fit_logcontrast_cv(
      y = y[train_id],
      x = log_X[train_id, , drop = FALSE],
      lambda_grid = ridge_lambda_grid,
      nfolds = nfolds,
      seed = seed + 1000L + fold_now
    )
    
    full_lc_pred[test_id] <- predict_constrained_ridge(
      full_lc_fit,
      log_X[test_id, , drop = FALSE]
    )
    
    # --------------------------------------------------------
    # Full random-forest predictor
    # --------------------------------------------------------
    set.seed(seed + 2000L + fold_now)
    full_rf_fit <- randomForest(
      x = X[train_id, , drop = FALSE],
      y = y[train_id],
      ntree = rf_ntree
    )
    
    full_rf_pred[test_id] <- predict(
      full_rf_fit,
      X[test_id, , drop = FALSE]
    )
    
    # --------------------------------------------------------
    # Evaluate the same fitted regressions at the perturbed compositions
    # --------------------------------------------------------
    for (k in seq_along(J_signal)) {
      W_full <- W_full_list[[k]]
      log_W_full <- log_W_full_list[[k]]
      
      plugin_lc_perturbed_pred[test_id, k] <-
        predict_constrained_ridge(
          full_lc_fit,
          log_W_full[test_id, , drop = FALSE]
        )
      
      plugin_rf_perturbed_pred[test_id, k] <- predict(
        full_rf_fit,
        W_full[test_id, , drop = FALSE]
      )
    }
  }
  
  return(list(
    full_lc = full_lc_pred,
    full_rf = full_rf_pred,
    plugin_lc_perturbed = plugin_lc_perturbed_pred,
    plugin_rf_perturbed = plugin_rf_perturbed_pred
  ))
}

compute_logcontrast_plugin_once <- function(X, y, target_j, nfolds, seed, pseud, ridge_lambda_grid) {
  n_now <- nrow(X)
  
  W <- X
  W[, target_j] <- 0
  remaining_mass <- rowSums(W)
  
  if (any(remaining_mass <= 1e-14)) {
    return(NA_real_)
  }
  
  W <- W / remaining_mass
  
  log_X <- log(X + pseud)
  log_W <- log(W + pseud)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n_now))
  
  original_pred <- rep(NA_real_, n_now)
  perturbed_pred <- rep(NA_real_, n_now)
  
  for (fold_now in seq_len(nfolds)) {
    test_id <- which(fold_id == fold_now)
    train_id <- setdiff(seq_len(n_now), test_id)
    
    fit_Z <- fit_logcontrast_cv(
      y = y[train_id],
      x = log_X[train_id, , drop = FALSE],
      lambda_grid = ridge_lambda_grid,
      nfolds = nfolds,
      seed = seed + 1000L + fold_now
    )
    
    original_pred[test_id] <- predict_constrained_ridge(
      fit_Z,
      log_X[test_id, , drop = FALSE]
    )
    
    perturbed_pred[test_id] <- predict_constrained_ridge(
      fit_Z,
      log_W[test_id, , drop = FALSE]
    )
  }
  
  mean((original_pred - perturbed_pred)^2) / var(y)
}

compute_random_forest_plugin_once <- function(X, y, target_j, nfolds, seed, rf_ntree) {
  n_now <- nrow(X)
  
  W <- X
  W[, target_j] <- 0
  remaining_mass <- rowSums(W)
  
  if (any(remaining_mass <= 1e-14)) {
    return(NA_real_)
  }
  
  W <- W / remaining_mass
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n_now))
  
  original_pred <- rep(NA_real_, n_now)
  perturbed_pred <- rep(NA_real_, n_now)
  
  for (fold_now in seq_len(nfolds)) {
    test_id <- which(fold_id == fold_now)
    train_id <- setdiff(seq_len(n_now), test_id)
    
    set.seed(seed + 1000L + fold_now)
    fit_Z <- randomForest(
      x = X[train_id, , drop = FALSE],
      y = y[train_id],
      ntree = rf_ntree
    )
    
    original_pred[test_id] <- predict(
      fit_Z,
      X[test_id, , drop = FALSE]
    )
    
    perturbed_pred[test_id] <- predict(
      fit_Z,
      W[test_id, , drop = FALSE]
    )
  }
  
  mean((original_pred - perturbed_pred)^2) / var(y)
}

bootstrap_logcontrast_plugin_test <- function(X, y, target_j, plugin_hat, c_threshold, alpha, n_boot, nfolds, seed, pseud, ridge_lambda_grid) {
  n_now <- nrow(X)
  bootstrap_estimates <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    set.seed(seed + b)
    bootstrap_id <- sample(seq_len(n_now), size = n_now, replace = TRUE)
    
    bootstrap_estimates[b] <- tryCatch(
      compute_logcontrast_plugin_once(
        X = X[bootstrap_id, , drop = FALSE],
        y = y[bootstrap_id],
        target_j = target_j,
        nfolds = nfolds,
        seed = seed + 10000L + b,
        pseud = pseud,
        ridge_lambda_grid = ridge_lambda_grid
      ),
      error = function(e) NA_real_
    )
  }
  
  bootstrap_estimates <- bootstrap_estimates[
    is.finite(bootstrap_estimates)
  ]
  
  if (length(bootstrap_estimates) < max(10L, ceiling(0.8 * n_boot))) {
    return(list(
      se = NA_real_,
      lower_bound = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      reject = NA_integer_,
      valid_bootstrap = length(bootstrap_estimates)
    ))
  }
  
  bootstrap_se <- sd(bootstrap_estimates)
  
  # One-sided basic bootstrap lower confidence bound:
  # hat(phi) - quantile_{1-alpha}{hat(phi)^* - hat(phi)}.
  bootstrap_error <- bootstrap_estimates - plugin_hat
  lower_bound <- plugin_hat - as.numeric(
    quantile(
      bootstrap_error,
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )
  )
  
  ci_lower <- plugin_hat - as.numeric(
    quantile(
      bootstrap_error,
      probs = 1 - alpha / 2,
      names = FALSE,
      type = 8
    )
  )
  
  ci_upper <- plugin_hat - as.numeric(
    quantile(
      bootstrap_error,
      probs = alpha / 2,
      names = FALSE,
      type = 8
    )
  )
  
  p_value <- mean(bootstrap_error >= plugin_hat - c_threshold)
  
  return(list(
    se = bootstrap_se,
    lower_bound = lower_bound,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    reject = as.integer(lower_bound > c_threshold),
    valid_bootstrap = length(bootstrap_estimates)
  ))
}

bootstrap_random_forest_plugin_test <- function(X, y, target_j, plugin_hat, c_threshold, alpha, n_boot, nfolds, seed, rf_ntree) {
  n_now <- nrow(X)
  bootstrap_estimates <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    set.seed(seed + b)
    bootstrap_id <- sample(seq_len(n_now), size = n_now, replace = TRUE)
    
    bootstrap_estimates[b] <- tryCatch(
      compute_random_forest_plugin_once(
        X = X[bootstrap_id, , drop = FALSE],
        y = y[bootstrap_id],
        target_j = target_j,
        nfolds = nfolds,
        seed = seed + 10000L + b,
        rf_ntree = rf_ntree
      ),
      error = function(e) NA_real_
    )
  }
  
  bootstrap_estimates <- bootstrap_estimates[
    is.finite(bootstrap_estimates)
  ]
  
  if (length(bootstrap_estimates) < max(10L, ceiling(0.8 * n_boot))) {
    return(list(
      se = NA_real_,
      lower_bound = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      reject = NA_integer_,
      valid_bootstrap = length(bootstrap_estimates)
    ))
  }
  
  bootstrap_error <- bootstrap_estimates - plugin_hat
  
  lower_bound <- plugin_hat - as.numeric(
    quantile(
      bootstrap_error,
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )
  )
  
  ci_lower <- plugin_hat - as.numeric(
    quantile(
      bootstrap_error,
      probs = 1 - alpha / 2,
      names = FALSE,
      type = 8
    )
  )
  
  ci_upper <- plugin_hat - as.numeric(
    quantile(
      bootstrap_error,
      probs = alpha / 2,
      names = FALSE,
      type = 8
    )
  )
  
  p_value <- mean(bootstrap_error >= plugin_hat - c_threshold)
  
  return(list(
    se = sd(bootstrap_estimates),
    lower_bound = lower_bound,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value,
    reject = as.integer(lower_bound > c_threshold),
    valid_bootstrap = length(bootstrap_estimates)
  ))
}

estimate_plugin_importance <- function(y, original_pred, perturbed_pred) {
  sigma_y_sq <- var(y)
  
  squared_prediction_change <- (original_pred - perturbed_pred)^2
  
  plugin_hat <- mean(squared_prediction_change) / sigma_y_sq
  
  return(list(
    estimate = plugin_hat
  ))
}

# ============================================================
# 5. VIPER positive-threshold tests
# ============================================================

get_threshold_test <- function(phi_hat, sigma_phi_sq, n, c_threshold, alpha) {
  if (!is.finite(phi_hat) ||
      !is.finite(sigma_phi_sq) ||
      sigma_phi_sq <= 0) {
    return(list(
      se = NA_real_,
      z = NA_real_,
      p_value = NA_real_,
      lower_bound = NA_real_,
      reject = NA_integer_
    ))
  }
  
  se_phi <- sqrt(sigma_phi_sq / n)
  
  z_stat <- (phi_hat - c_threshold) / se_phi
  
  # One-sided p-value for H0: phi <= c versus H1: phi > c
  p_value <- pnorm(z_stat, lower.tail = FALSE)
  
  # One-sided (1-alpha) lower confidence bound
  lower_bound <- phi_hat - qnorm(1 - alpha) * se_phi
  
  reject <- as.integer(lower_bound > c_threshold)
  
  return(list(
    se = se_phi,
    z = z_stat,
    p_value = p_value,
    lower_bound = lower_bound,
    reject = reject
  ))
}

# ============================================================
# 6. Fit the target component
# ============================================================

fit_threshold_tests <- function(X, y, J_signal, effect_value, iter_seed, c_threshold, alpha) {
  m <- length(J_signal)
  n_current <- nrow(X)
  
  output <- data.frame(
    component = J_signal,
    effect_value = effect_value,
    theta_hat = NA_real_,
    phi_hat = NA_real_,
    sigma_phi_sq = NA_real_,
    se_phi = NA_real_,
    q_hat = NA_real_,
    z_threshold = NA_real_,
    p_threshold = NA_real_,
    lower_bound = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    reject_threshold = NA_integer_
  )
  
  for (k in seq_along(J_signal)) {
    j <- J_signal[k]
    
    rep_X <- reparameterize(X = X, s = j)
    
    fit_j <- tryCatch(
      R2_partial_linear(
        Y = y,
        L = rep_X$L,
        W = rep_X$W,
        s = j,
        nfolds = nfolds,
        seed = iter_seed + j
      ),
      error = function(e) {
        message("Variable ", j, " failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(fit_j)) {
      test_j <- get_threshold_test(
        phi_hat = fit_j$phi_hat,
        sigma_phi_sq = fit_j$sigma_phi_sq,
        n = n_current,
        c_threshold = c_threshold,
        alpha = alpha
      )
      
      output$theta_hat[k] <- fit_j$theta_hat
      output$phi_hat[k] <- fit_j$phi_hat
      output$sigma_phi_sq[k] <- fit_j$sigma_phi_sq
      output$se_phi[k] <- test_j$se
      output$q_hat[k] <- mean((1 - rep_X$L)^2)
      output$z_threshold[k] <- test_j$z
      output$p_threshold[k] <- test_j$p_value
      output$lower_bound[k] <- test_j$lower_bound
      output$ci_lower[k] <- fit_j$phi_hat -
        qnorm(1 - alpha / 2) * test_j$se
      output$ci_upper[k] <- fit_j$phi_hat +
        qnorm(1 - alpha / 2) * test_j$se
      output$reject_threshold[k] <- test_j$reject
    }
  }
  
  return(output)
}

# ============================================================
# 7. Run simulation
# ============================================================

# Parallelize the outer Monte Carlo repetitions.  Each worker runs one
# complete simulation repetition, including both full-refit bootstraps.
# Keeping the bootstrap loops sequential avoids nested parallelism.
n_cores <- 15

# Prevent threaded BLAS/OpenMP code inside each worker from creating another
# layer of threads.  This is important when 15 R workers run simultaneously.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

cl <- parallel::makeCluster(n_cores)

parallel::clusterExport(cl, "algorithm_file")
parallel::clusterEvalQ(cl, {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1"
  )
  library(MASS)
  library(randomForest)
  source(algorithm_file)
  NULL
})

doParallel::registerDoParallel(cl)
cat("Parallel simulation using", n_cores, "workers.\n")

cat("Constructing independent large-sample oracle designs...\n")

oracle_dense_list <- list()
oracle_selected_list <- list()
oracle_id <- 1L

for (scenario_now in scenario_names) {
  effect_candidates <- if (scenario_now == "correct_logcontrast") {
    effect_candidates_correct
  } else {
    effect_candidates_nonlinear
  }
  
  for (q_now in q_value) {
    oracle_dense_now <- build_oracle_design(
      scenario = scenario_now,
      q_value = q_now,
      effect_candidates = effect_candidates,
      n_oracle = n_oracle,
      seed = base_seed + 9000000L + oracle_id * 10000L
    )
    
    oracle_selected_now <- select_effect_grid(
      oracle_table = oracle_dense_now,
      threshold = c_threshold,
      multipliers = oracle_level_multipliers
    )
    
    oracle_dense_list[[oracle_id]] <- oracle_dense_now
    oracle_selected_list[[oracle_id]] <- oracle_selected_now
    oracle_id <- oracle_id + 1L
  }
}

oracle_dense_results <- do.call(rbind, oracle_dense_list)
oracle_design <- do.call(rbind, oracle_selected_list)
rownames(oracle_dense_results) <- NULL
rownames(oracle_design) <- NULL

print(oracle_design, row.names = FALSE, digits = 4)

all_results <- list()

begin_time <- Sys.time()

for (scenario_now in scenario_names) {
  for (q_now in q_value) {
    design_now <- oracle_design[
      oracle_design$scenario == scenario_now &
        oracle_design$q == q_now,
    ]
    effect_grid <- design_now$effect_value
    
    for (effect_id in seq_along(effect_grid)) {
      effect_now <- effect_grid[effect_id]
      oracle_phi_now <- design_now$oracle_phi[effect_id]
      threshold_region_now <- design_now$threshold_region[effect_id]
      
      cat(sprintf(
        "Starting scenario=%s, n=%d, q=%.1f, effect=%.3f, oracle phi=%.4f\n",
        scenario_now,
        n,
        q_now,
        effect_now,
        oracle_phi_now
      ))
      
      effect_results <- foreach(
        iter = seq_len(n_iter),
        .packages = c("MASS", "randomForest"),
        .export = c(
          "generate_logistic_normal",
          "make_beta",
          "get_true_regression",
          "generate_response",
          "fit_constrained_ridge",
          "predict_constrained_ridge",
          "fit_logcontrast_cv",
          "compute_plugin_predictions",
          "compute_logcontrast_plugin_once",
          "compute_random_forest_plugin_once",
          "bootstrap_logcontrast_plugin_test",
          "bootstrap_random_forest_plugin_test",
          "estimate_plugin_importance",
          "get_threshold_test",
          "fit_threshold_tests"
        ),
        .errorhandling = "pass",
        .inorder = TRUE
      ) %dopar% {
        
        data_seed <-
          base_seed +
          n * 100L +
          round(q_now * 10000) +
          effect_id * 100000L +
          iter
        
        set.seed(data_seed)
        
        X <- generate_logistic_normal(
          n = n,
          p = p,
          s = s,
          rho = rho,
          presence_prob = q_now
        )
        
        response_now <- generate_response(
          X = X,
          target_j = target_j,
          effect_value = effect_now,
          scenario = scenario_now,
          reference_j = reference_j,
          sig = sig,
          pseud = pseud
        )
        
        y <- response_now$y
        
        plugin_predictions <- compute_plugin_predictions(
          X = X,
          y = y,
          J_signal = J_signal,
          nfolds = nfolds,
          seed = data_seed + 250000L,
          pseud = pseud,
          ridge_lambda_grid = ridge_lambda_grid,
          rf_ntree = rf_ntree
        )
        
        fit_now <- fit_threshold_tests(
          X = X,
          y = y,
          J_signal = J_signal,
          effect_value = effect_now,
          iter_seed = data_seed + 500000L,
          c_threshold = c_threshold,
          alpha = alpha
        )
        
        fit_now$plugin_logcontrast_hat <- NA_real_
        fit_now$plugin_rf_hat <- NA_real_
        
        fit_now$plugin_logcontrast_bootstrap_se <- NA_real_
        fit_now$plugin_logcontrast_bootstrap_lower <- NA_real_
        fit_now$plugin_logcontrast_bootstrap_ci_lower <- NA_real_
        fit_now$plugin_logcontrast_bootstrap_ci_upper <- NA_real_
        fit_now$plugin_logcontrast_bootstrap_p_threshold <- NA_real_
        fit_now$plugin_logcontrast_bootstrap_reject_threshold <- NA_integer_
        fit_now$plugin_logcontrast_valid_bootstrap <- NA_integer_
        
        fit_now$plugin_rf_bootstrap_se <- NA_real_
        fit_now$plugin_rf_bootstrap_lower <- NA_real_
        fit_now$plugin_rf_bootstrap_ci_lower <- NA_real_
        fit_now$plugin_rf_bootstrap_ci_upper <- NA_real_
        fit_now$plugin_rf_bootstrap_p_threshold <- NA_real_
        fit_now$plugin_rf_bootstrap_reject_threshold <- NA_integer_
        fit_now$plugin_rf_valid_bootstrap <- NA_integer_
        
        for (k in seq_along(J_signal)) {
          plugin_lc_k <- estimate_plugin_importance(
            y = y,
            original_pred = plugin_predictions$full_lc,
            perturbed_pred = plugin_predictions$plugin_lc_perturbed[, k]
          )
          
          plugin_rf_k <- estimate_plugin_importance(
            y = y,
            original_pred = plugin_predictions$full_rf,
            perturbed_pred = plugin_predictions$plugin_rf_perturbed[, k]
          )
          
          plugin_lc_bootstrap_k <- bootstrap_logcontrast_plugin_test(
            X = X,
            y = y,
            target_j = J_signal[k],
            plugin_hat = plugin_lc_k$estimate,
            c_threshold = c_threshold,
            alpha = alpha,
            n_boot = n_boot,
            nfolds = nfolds,
            seed = data_seed + 750000L + k * 1000L,
            pseud = pseud,
            ridge_lambda_grid = ridge_lambda_grid
          )
          
          plugin_rf_bootstrap_k <- bootstrap_random_forest_plugin_test(
            X = X,
            y = y,
            target_j = J_signal[k],
            plugin_hat = plugin_rf_k$estimate,
            c_threshold = c_threshold,
            alpha = alpha,
            n_boot = n_boot,
            nfolds = nfolds,
            seed = data_seed + 850000L + k * 1000L,
            rf_ntree = rf_ntree
          )
          
          fit_now$plugin_logcontrast_hat[k] <- plugin_lc_k$estimate
          fit_now$plugin_rf_hat[k] <- plugin_rf_k$estimate
          
          fit_now$plugin_logcontrast_bootstrap_se[k] <-
            plugin_lc_bootstrap_k$se
          fit_now$plugin_logcontrast_bootstrap_lower[k] <-
            plugin_lc_bootstrap_k$lower_bound
          fit_now$plugin_logcontrast_bootstrap_ci_lower[k] <-
            plugin_lc_bootstrap_k$ci_lower
          fit_now$plugin_logcontrast_bootstrap_ci_upper[k] <-
            plugin_lc_bootstrap_k$ci_upper
          fit_now$plugin_logcontrast_bootstrap_p_threshold[k] <-
            plugin_lc_bootstrap_k$p_value
          fit_now$plugin_logcontrast_bootstrap_reject_threshold[k] <-
            plugin_lc_bootstrap_k$reject
          fit_now$plugin_logcontrast_valid_bootstrap[k] <-
            plugin_lc_bootstrap_k$valid_bootstrap
          
          fit_now$plugin_rf_bootstrap_se[k] <-
            plugin_rf_bootstrap_k$se
          fit_now$plugin_rf_bootstrap_lower[k] <-
            plugin_rf_bootstrap_k$lower_bound
          fit_now$plugin_rf_bootstrap_ci_lower[k] <-
            plugin_rf_bootstrap_k$ci_lower
          fit_now$plugin_rf_bootstrap_ci_upper[k] <-
            plugin_rf_bootstrap_k$ci_upper
          fit_now$plugin_rf_bootstrap_p_threshold[k] <-
            plugin_rf_bootstrap_k$p_value
          fit_now$plugin_rf_bootstrap_reject_threshold[k] <-
            plugin_rf_bootstrap_k$reject
          fit_now$plugin_rf_valid_bootstrap[k] <-
            plugin_rf_bootstrap_k$valid_bootstrap
          
        }
        
        fit_now$q <- q_now
        fit_now$n <- n
        fit_now$scenario <- scenario_now
        fit_now$oracle_phi <- oracle_phi_now
        fit_now$threshold_region <- threshold_region_now
        fit_now$iteration <- iter
        
        fit_now
      }
      
      failed_parallel <- vapply(
        effect_results,
        inherits,
        logical(1),
        what = "error"
      )
      
      if (any(failed_parallel)) {
        warning(
          sum(failed_parallel),
          " parallel iterations failed for scenario=",
          scenario_now,
          ", effect=",
          effect_now
        )
      }
      
      effect_results <- effect_results[!failed_parallel]
      all_results <- c(all_results, effect_results)
      
      cat(sprintf(
        "Completed scenario=%s, effect=%.3f: %d/%d successful repetitions\n",
        scenario_now,
        effect_now,
        length(effect_results),
        n_iter
      ))
      cat("Total time:", as.character(Sys.time() - begin_time), "\n")
    }
  }
}

parallel::stopCluster(cl)
foreach::registerDoSEQ()

component_results <- do.call(rbind, all_results)
rownames(component_results) <- NULL


# ============================================================
# 8. Summarize threshold rejection probabilities
# ============================================================

summary_results <- do.call(rbind, lapply(
  split(component_results, list(component_results$scenario,
                                component_results$n, component_results$q,
                                component_results$effect_value), drop = TRUE),
  function(x) {
    oracle_phi <- unique(x$oracle_phi)
    coverage <- function(lower, upper) {
      if (oracle_phi <= 1e-10) return(NA_real_)
      mean(lower <= oracle_phi & oracle_phi <= upper, na.rm = TRUE)
    }
    
    data.frame(
      scenario = unique(x$scenario),
      q = unique(x$q),
      effect_value = unique(x$effect_value),
      oracle_phi = oracle_phi,
      mean_phi_hat = mean(x$phi_hat, na.rm = TRUE),
      mean_logcontrast_plugin = mean(x$plugin_logcontrast_hat, na.rm = TRUE),
      mean_random_forest_plugin = mean(x$plugin_rf_hat, na.rm = TRUE),
      viper_coverage = coverage(x$ci_lower, x$ci_upper),
      logcontrast_plugin_bootstrap_coverage = coverage(
        x$plugin_logcontrast_bootstrap_ci_lower,
        x$plugin_logcontrast_bootstrap_ci_upper
      ),
      random_forest_plugin_bootstrap_coverage = coverage(
        x$plugin_rf_bootstrap_ci_lower,
        x$plugin_rf_bootstrap_ci_upper
      ),
      threshold_rejection = mean(x$reject_threshold, na.rm = TRUE),
      logcontrast_plugin_bootstrap_rejection = mean(
        x$plugin_logcontrast_bootstrap_reject_threshold, na.rm = TRUE
      ),
      random_forest_plugin_bootstrap_rejection = mean(
        x$plugin_rf_bootstrap_reject_threshold, na.rm = TRUE
      )
    )
  }
))
# ============================================================
# 10. Main figure: mean importance, coverage, and rejection
# ============================================================
graphics.off()

method_col <- c(
  Oracle = "black",
  VIPER = "#0072B2",
  LC = "#D55E00",
  RF = "#009E73"
)

method_lty <- c(
  Oracle = 2,
  VIPER = 1,
  LC = 5,
  RF = 4
)

method_pch <- c(
  Oracle = 4,
  VIPER = 16,
  LC = 17,
  RF = 15
)

scenario_title <- c(
  correct_logcontrast = "Log-contrast model",
  nonlinear_plm = "Partial linear model"
)

get_plot_data <- function(scenario_now) {
  temp <- summary_results[
    summary_results$scenario == scenario_now &
      summary_results$q == q_value,
  ]
  temp <- temp[order(temp$effect_value), ]
  
  dense_now <- oracle_dense_results[
    oracle_dense_results$scenario == scenario_now &
      oracle_dense_results$q == q_value,
  ]
  
  list(
    temp = temp,
    threshold_effect = interpolate_threshold_effect(
      dense_now,
      c_threshold
    ),
    x_axis_label = if (scenario_now == "correct_logcontrast") {
      expression(beta[j])
    } else {
      expression(theta[j])
    }
  )
}

draw_mean_panel <- function(scenario_now, show_legend = FALSE, show_title = TRUE) {
  plot_data <- get_plot_data(scenario_now)
  temp <- plot_data$temp
  
  ylim_phi <- range(
    c(
      temp$oracle_phi,
      temp$mean_phi_hat,
      temp$mean_logcontrast_plugin,
      temp$mean_random_forest_plugin,
      c_threshold
    ),
    finite = TRUE
  )
  
  plot(
    temp$effect_value,
    temp$oracle_phi,
    type = "b",
    pch = method_pch["Oracle"],
    lty = method_lty["Oracle"],
    lwd = 2,
    col = method_col["Oracle"],
    ylim = ylim_phi,
    xlab = plot_data$x_axis_label,
    ylab = "Mean importance",
    main = if (show_title) {
      scenario_title[[scenario_now]]
    } else {
      ""
    }
  )
  
  lines(
    temp$effect_value,
    temp$mean_phi_hat,
    type = "b",
    pch = method_pch["VIPER"],
    lty = method_lty["VIPER"],
    lwd = 2,
    col = method_col["VIPER"]
  )
  
  lines(
    temp$effect_value,
    temp$mean_logcontrast_plugin,
    type = "b",
    pch = method_pch["LC"],
    lty = method_lty["LC"],
    lwd = 2,
    col = method_col["LC"]
  )
  
  lines(
    temp$effect_value,
    temp$mean_random_forest_plugin,
    type = "b",
    pch = method_pch["RF"],
    lty = method_lty["RF"],
    lwd = 2,
    col = method_col["RF"]
  )
  
  # abline(
  #   h = c_threshold,
  #   lty = 3,
  #   lwd = 1.5,
  #   col = "gray35"
  # )
  
  # if (is.finite(plot_data$threshold_effect)) {
  #   abline(
  #     v = plot_data$threshold_effect,
  #     lty = 3,
  #     lwd = 1.5,
  #     col = "gray55"
  #   )
  # }
  
  grid()
  
  if (show_legend) {
    legend(
      "topleft",
      legend = c(
        "Reference",
        "VIPER",
        "LC plug-in",
        "RF plug-in"
      ),
      col = c(
        method_col["Oracle"],
        method_col["VIPER"],
        method_col["LC"],
        method_col["RF"]
      ),
      pch = c(
        method_pch["Oracle"],
        method_pch["VIPER"],
        method_pch["LC"],
        method_pch["RF"]
      ),
      lty = c(
        method_lty["Oracle"],
        method_lty["VIPER"],
        method_lty["LC"],
        method_lty["RF"]
      ),
      lwd = c(2, 2, 2, 2),
      bty = "n",
      cex = 1
    )
  }
}

draw_coverage_panel <- function(scenario_now, show_legend = FALSE, show_title = TRUE) {
  plot_data <- get_plot_data(scenario_now)
  temp <- plot_data$temp
  
  plot(
    temp$effect_value,
    temp$viper_coverage,
    type = "b",
    pch = method_pch["VIPER"],
    lty = method_lty["VIPER"],
    lwd = 2,
    col = method_col["VIPER"],
    ylim = c(0, 1),
    xlab = plot_data$x_axis_label,
    ylab = "Empirical coverage",
    main = if (show_title) {
      scenario_title[[scenario_now]]
    } else {
      ""
    }
  )
  
  lines(
    temp$effect_value,
    temp$logcontrast_plugin_bootstrap_coverage,
    type = "b",
    pch = method_pch["LC"],
    lty = method_lty["LC"],
    lwd = 2,
    col = method_col["LC"]
  )
  
  lines(
    temp$effect_value,
    temp$random_forest_plugin_bootstrap_coverage,
    type = "b",
    pch = method_pch["RF"],
    lty = method_lty["RF"],
    lwd = 2,
    col = method_col["RF"]
  )
  
  abline(
    h = 1 - alpha,
    lty = 3,
    lwd = 1.5,
    col = "gray35"
  )
  
  # if (is.finite(plot_data$threshold_effect)) {
  #   abline(
  #     v = plot_data$threshold_effect,
  #     lty = 3,
  #     lwd = 1.5,
  #     col = "gray55"
  #   )
  # }
  
  grid()
  
  if (show_legend) {
    legend(
      "bottomleft",
      legend = c(
        "VIPER",
        "LC bootstrap",
        "RF bootstrap"
      ),
      col = c(
        method_col["VIPER"],
        method_col["LC"],
        method_col["RF"]
      ),
      pch = c(
        method_pch["VIPER"],
        method_pch["LC"],
        method_pch["RF"]
      ),
      lty = c(
        method_lty["VIPER"],
        method_lty["LC"],
        method_lty["RF"]
      ),
      lwd = c(2, 2, 2),
      bty = "n",
      cex = 1
    )
  }
}

draw_rejection_panel <- function(scenario_now, show_legend = FALSE) {
  plot_data <- get_plot_data(scenario_now)
  temp <- plot_data$temp
  
  plot(
    temp$effect_value,
    temp$threshold_rejection,
    type = "b",
    pch = method_pch["VIPER"],
    lty = method_lty["VIPER"],
    lwd = 2,
    col = method_col["VIPER"],
    ylim = c(0, 1),
    xlab = plot_data$x_axis_label,
    ylab = "Empirical rejection probability",
    main = scenario_title[[scenario_now]]
  )
  
  lines(
    temp$effect_value,
    temp$logcontrast_plugin_bootstrap_rejection,
    type = "b",
    pch = method_pch["LC"],
    lty = method_lty["LC"],
    lwd = 2,
    col = method_col["LC"]
  )
  
  lines(
    temp$effect_value,
    temp$random_forest_plugin_bootstrap_rejection,
    type = "b",
    pch = method_pch["RF"],
    lty = method_lty["RF"],
    lwd = 2,
    col = method_col["RF"]
  )
  
  abline(
    h = alpha,
    lty = 3,
    lwd = 1.5,
    col = "gray35"
  )
  
  if (is.finite(plot_data$threshold_effect)) {
    abline(
      v = plot_data$threshold_effect,
      lty = 3,
      lwd = 1.5,
      col = "gray55"
    )
  }
  
  grid()
  
  if (show_legend) {
    legend(
      "topleft",
      legend = c(
        "VIPER",
        "LC bootstrap",
        "RF bootstrap"
      ),
      col = c(
        method_col["VIPER"],
        method_col["LC"],
        method_col["RF"]
      ),
      pch = c(
        method_pch["VIPER"],
        method_pch["LC"],
        method_pch["RF"]
      ),
      lty = c(
        method_lty["VIPER"],
        method_lty["LC"],
        method_lty["RF"]
      ),
      lwd = c(2, 2, 2),
      bty = "n",
      cex = 1
    )
  }
}

# ============================================================
# Figure 1 (2 x 2):
# Columns: log-contrast model and partial linear model
# Rows: mean importance and empirical coverage
# ============================================================

draw_figure_1 <- function() {
  par(
    mfrow = c(2, 2),
    mar = c(4.5, 4.8, 3.0, 1.0),
    oma = c(0, 0, 0, 0)
  )
  
  draw_mean_panel(
    "correct_logcontrast",
    show_legend = TRUE,
    show_title = TRUE
  )
  
  draw_mean_panel(
    "nonlinear_plm",
    show_legend = FALSE,
    show_title = TRUE
  )
  
  draw_coverage_panel(
    "correct_logcontrast",
    show_legend = TRUE,
    show_title = TRUE
  )
  
  draw_coverage_panel(
    "nonlinear_plm",
    show_legend = FALSE,
    show_title = TRUE
  )
}

# ============================================================
# Figure 2:
# Positive-threshold rejection probabilities
# Left: log-contrast model
# Right: partial linear model
# ============================================================

draw_figure_2 <- function() {
  par(
    mfrow = c(1, length(scenario_names)),
    mar = c(4.5, 4.8, 3.0, 1.0),
    oma = c(0, 0, 0, 0)
  )
  
  for (scenario_now in scenario_names) {
    draw_rejection_panel(
      scenario_now,
      show_legend = scenario_now == scenario_names[1]
    )
  }
}

# ============================================================
# Save Figure 1
# ============================================================

figure_1_png <- sprintf(
  "VIPER_two_models_estimation_coverage_n=%d_q=%.1f_c=%.3f.png",
  n,
  q_value,
  c_threshold
)

png(
  filename = figure_1_png,
  width = 1500,
  height = 1350,
  res = 150
)

draw_figure_1()

dev.off()

# ============================================================
# Save Figure 2
# ============================================================

figure_2_png <- sprintf("VIPER_threshold_rejection_n=%d_q=%.1f_c=%.3f.png", n, q_value, c_threshold)

png(filename = figure_2_png,width = 1500, height = 750, res = 150)

draw_figure_2()

dev.off()

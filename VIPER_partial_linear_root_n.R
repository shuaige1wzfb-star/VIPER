library(MASS)
library(randomForest)
library(parallel)
library(doParallel)
library(foreach)

# Run from a working directory containing R2.R.
algorithm_file <- normalizePath("R2.R", mustWork = TRUE)
source(algorithm_file)

set.seed(123)

## ------------------------------------------------------------------
## Simulation settings
## ------------------------------------------------------------------
p <- 50
rho <- 0.5
presence_prob <- 0.2
sig <- 0.7
pseud <- 1e-10
target_dim <- 1

n_folds <- 3
n_reps <- 1000
n_cores <- 15
n_large <- 1e5 
sample_sizes <- c(50, 1000, 2000, 3000, 4000)
 
theta_plm <- 1

mu <- rep(0, p)
Sigma <- toeplitz(rho^(0:(p - 1)))

## ------------------------------------------------------------------
## Data-generating mechanism
## ------------------------------------------------------------------
generate_logistic_normal <- function(n, p, s, presence_prob) {   
  latent <- MASS::mvrnorm(
    n = n,
    mu = mu,
    Sigma = Sigma
  )
  
  composition <- exp(latent)
  composition <- composition / rowSums(composition)
  
  presence <- matrix(
    rbinom(n * s, size = 1, prob = presence_prob),
    nrow = n,
    ncol = s
  )
  composition[, seq_len(s)] <- composition[, seq_len(s), drop = FALSE] * presence
  
  row_total <- rowSums(composition)
  zero_rows <- row_total <= 0
  if (any(zero_rows)) stop("Unexpected all-zero composition.")
  
  composition / row_total
}

partial_linear_nuisance <- function(W, pseud) {
  log_W <- log(W + pseud)
  
  balance_1 <- log_W[, 2] - log_W[, 3]
  balance_2 <- log_W[, 4] - log_W[, 5]
  balance_3 <- log_W[, 6] - log_W[, 7]
  
  0.8 * sin(balance_1) + 0.5 * tanh(balance_2)^2 + 0.4 * tanh(balance_1) * tanh(balance_3)
}

get_true_regression <- function(X, target_dim, theta_plm, pseud) {
  reparam <- reparameterize(X = X, s = target_dim)
  L <- reparam$L
  W <- reparam$W
  
  h_value <- partial_linear_nuisance(W, pseud)
  f_original <- theta_plm * L + h_value
  f_perturbed <- theta_plm + h_value
  
  list(
    L = L,
    W = W,
    conditional_mean = f_original,
    perturbed_mean = f_perturbed,
    perturbation_contrast = f_perturbed - f_original
  )
}

## ------------------------------------------------------------------
## Numerical reference value
## ------------------------------------------------------------------
set.seed(20261021)
X_oracle <- generate_logistic_normal(
  n = n_large,
  p = p,
  s = s,
  presence_prob = presence_prob
)

truth_oracle <- get_true_regression(
  X = X_oracle,
  target_dim = target_dim,
  theta_plm = theta_plm,
  pseud = pseud
)

oracle_error <- rnorm(n_large, mean = 0, sd = sig)
Y_oracle <- truth_oracle$conditional_mean + oracle_error

oracle_numerator <- mean(truth_oracle$perturbation_contrast^2)
oracle_sigma_y_sq <- var(Y_oracle)
oracle_phi <- oracle_numerator / oracle_sigma_y_sq

## ------------------------------------------------------------------
## One finite-sample replicate
## ------------------------------------------------------------------
run_one_replication <- function(n_now, replication_id) {
  replication_seed <- 100L * n_now + replication_id
  set.seed(replication_seed)
  
  X <- generate_logistic_normal(
    n = n_now,
    p = p,
    s = s,
    presence_prob = presence_prob
  )
  
  truth <- get_true_regression(
    X = X,
    target_dim = target_dim,
    theta_plm = theta_plm,
    pseud = pseud
  )
  
  y <- truth$conditional_mean + rnorm(n_now, mean = 0, sd = sig)
    
  fit <- tryCatch(
    R2_partial_linear2(
      Y = y,
      L = truth$L,
      W = truth$W,  
      s = target_dim,  
      nfolds = n_folds,
      seed = replication_seed
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(data.frame(
      n = n_now,
      replication = replication_id,
      phi_hat = NA_real_
    ))
  }
  
  
  data.frame(
    n = n_now,
    replication = replication_id,
    phi_hat = fit$phi_hat
  )
}

## ------------------------------------------------------------------
## Monte Carlo experiment
## ------------------------------------------------------------------
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

clusterExport(cluster, "algorithm_file")
clusterEvalQ(cluster, {
  library(MASS)
  library(randomForest)
  source(algorithm_file)
  NULL
})

simulation_by_n <- vector("list", length(sample_sizes))
overall_start_time <- Sys.time()

for (n_index in seq_along(sample_sizes)) {
  n_now <- sample_sizes[n_index]
  n_start_time <- Sys.time()
  
  cat(
    sprintf(
      "\n[%d/%d] Starting n = %d at %s\n",
      n_index,
      length(sample_sizes),
      n_now,
      format(n_start_time, "%Y-%m-%d %H:%M:%S")
    )
  )
  flush.console()
  
  simulation_by_n[[n_index]] <- foreach(
    replication_id = seq_len(n_reps),
    .combine = rbind,
    .packages = c("MASS", "randomForest"),
    .export = c(
      "generate_logistic_normal",
      "partial_linear_nuisance",
      "get_true_regression",
      "run_one_replication",
      "p",
      "s",
      "presence_prob",
      "sig",
      "pseud",
      "target_dim",
      "n_folds",
      "theta_plm",
      "mu",
      "Sigma"
    )
  ) %dopar% {
    run_one_replication(
      n_now = n_now,
      replication_id = replication_id
    )
  }
  
  n_end_time <- Sys.time()
  elapsed_seconds <- as.numeric(
    difftime(n_end_time, n_start_time, units = "secs")
  )
  
  cat(
    sprintf(
      "[%d/%d] Finished n = %d; elapsed = %.2f minutes\n",
      n_index,
      length(sample_sizes),
      n_now,
      elapsed_seconds / 60
    )
  )
  flush.console()
}

simulation_results <- do.call(rbind, simulation_by_n)

stopCluster(cluster)
registerDoSEQ()

cat(
  sprintf(
    "\nAll sample sizes completed; total elapsed = %.2f minutes.\n",
    as.numeric(difftime(
      Sys.time(),
      overall_start_time,
      units = "mins"
    ))
  )
)

## ------------------------------------------------------------------
## Bias and variance summaries
## ------------------------------------------------------------------
summary_results <- do.call(
  rbind,
  lapply(sample_sizes, function(n_now) {
    results_now <- simulation_results[
      simulation_results$n == n_now &
        is.finite(simulation_results$phi_hat),
    ]
    
    estimates <- results_now$phi_hat
    mean_hat <- mean(estimates)
    bias_hat <- mean_hat - oracle_phi
    variance_hat <- var(estimates)
    bias_mcse <- sd(estimates) / sqrt(length(estimates))
    sqrt_n_bias_hat <- sqrt(n_now) * bias_hat
    sqrt_n_bias_mcse <- sqrt(n_now) * bias_mcse
    z_mc <- qnorm(0.975)
    
    data.frame(
      n = n_now,
      sqrt_n_bias = sqrt_n_bias_hat,
      sqrt_n_bias_lower = sqrt_n_bias_hat - z_mc * sqrt_n_bias_mcse,
      sqrt_n_bias_upper = sqrt_n_bias_hat + z_mc * sqrt_n_bias_mcse,
      n_variance = n_now * variance_hat
    )
  })
)

## ------------------------------------------------------------------
## Plotting
## ------------------------------------------------------------------
draw_convergence_figure <- function() {
  temp <- summary_results[order(summary_results$n), ]
  
  par(
    mfrow = c(1, 2),
    mar = c(4.5, 4.8, 3.0, 1.0),
    oma = c(0, 0, 0, 0)
  )
  
  plot(
    temp$n,
    temp$sqrt_n_bias,
    type = "b",
    pch = 16,
    lwd = 2,
    col = "#0072B2",
    ylim = c(-0.5, 0.5),
    xlab = "Sample size",
    ylab = expression(sqrt(n) %.% Bias(hat(phi)[1])),
    main = "Bias"
  )
  arrows(
    x0 = temp$n,
    y0 = temp$sqrt_n_bias_lower,
    x1 = temp$n,
    y1 = temp$sqrt_n_bias_upper,
    angle = 90,
    code = 3,
    length = 0.05,
    lwd = 1.2,
    col = "#0072B2"
  )
  abline(h = 0, lty = 3, lwd = 1.5, col = "#0072B2")
  grid()
  
  plot(
    temp$n,
    temp$n_variance,
    type = "b",
    pch = 16,
    lwd = 2,
    col = "#0072B2",
    xlab = "Sample size",
    ylab = expression(n %.% Var(hat(phi)[1])),
    main = "Variance"
  )
  grid()
}

png(
  filename = "VIPER_partial_linear_root_n.png",
  width = 1500,
  height = 750,
  res = 150
)
draw_convergence_figure()
dev.off()

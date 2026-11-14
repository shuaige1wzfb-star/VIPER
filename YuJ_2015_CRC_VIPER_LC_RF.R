# YuJ_2015 CRC: all-target VIPER / LC / RF inference and one top-10 figure.
# Data: the dated curatedMetagenomicData resource is fetched via ExperimentHub.
# ExperimentHub stores the download in its local cache.
# Sample metadata come from the installed curatedMetagenomicData package.
# VIPER uses the existing R2.R learners unchanged.
# Inference is pointwise at alpha = 0.05 against c = 0.01.
# LC/RF fit on original Z and evaluate that SAME fitted model at W_j.
# Displayed intervals are pointwise, not post-selection intervals.

rm(list = ls())
set.seed(2026)

suppressPackageStartupMessages({
  library(randomForest)
  library(doParallel)
  library(foreach)
  library(ggplot2)
})

# ------------------------- Settings -------------------------
# Run from the directory containing R2.R.
algorithm_file <- normalizePath("R2.R", mustWork = TRUE)
source(algorithm_file)

nfolds <- 3
base_seed <- 2026
top_k <- 10
alpha <- 0.05
c_threshold <- 0.01
q_lower <- 0.10
q_upper <- 0.90
ncores <- 8

plugin_pseud <- 1e-10
plugin_ridge_grid <- 10^seq(-4, 8, by = 1)
plugin_rf_ntree <- 300
plugin_n_boot <- 200L

# ============================================================================
# 1. Retrieve the abundance resource through ExperimentHub
# ============================================================================

resource_pattern <- "^2021-03-31[.]YuJ_2015[.]relative_abundance$"
cat("Retrieving YuJ_2015 relative abundance through ExperimentHub...\n")
abundance_resource <- curatedMetagenomicData::curatedMetagenomicData(
  pattern = resource_pattern,
  dryrun = FALSE,
  counts = FALSE,
  rownames = "long"
)

raw_abundance <- as.matrix(SummarizedExperiment::assay(abundance_resource[[1L]]))

# The downloaded matrix contains several taxonomic ranks.
# Keep species-level entries only; exclude higher ranks and strain-level rows.
species_rows <- grepl(
  "(^|\\|)s__[^|]+$",
  rownames(raw_abundance)
)


# When available, apply the same reference-tree compatibility restriction
# used by curatedMetagenomicData.
reference_tree <- get0(
  "phylogeneticTree",
  envir = asNamespace("curatedMetagenomicData"),
  inherits = FALSE
)

if (!is.null(reference_tree) &&
    !is.null(reference_tree$tip.label)) {
  unmatched_species <- sum(
    species_rows &
      !rownames(raw_abundance) %in% reference_tree$tip.label
  )
  
  if (unmatched_species > 0L) {
    cat(sprintf(
      "Species without package reference-tree matches excluded: %d\n",
      unmatched_species
    ))
  }
  
  species_rows <- species_rows &
    rownames(raw_abundance) %in% reference_tree$tip.label
}


abundance_species <- raw_abundance[
  species_rows,
  ,
  drop = FALSE
]
storage.mode(abundance_species) <- "numeric"

# ============================================================================
# 1B. Match the locally available sample metadata
# ============================================================================

crc_metadata <- subset(
  curatedMetagenomicData::sampleMetadata,
  study_name == "YuJ_2015" &
    study_condition %in% c("CRC", "control") &
    body_site == "stool"
)

if (nrow(crc_metadata) == 0L) {
  stop("No YuJ_2015 CRC/control stool samples were found in sampleMetadata.")
}

required_metadata <- c(
  "sample_id", "subject_id", "study_condition"
)

if (!all(required_metadata %in% names(crc_metadata))) {
  stop("Required sample metadata fields are missing.")
}

crc_metadata$sample_id <- as.character(crc_metadata$sample_id)
crc_metadata$subject_id <- trimws(
  as.character(crc_metadata$subject_id)
)

if (anyNA(crc_metadata$sample_id) ||
    any(!nzchar(crc_metadata$sample_id)) ||
    anyDuplicated(crc_metadata$sample_id)) {
  stop("Missing or duplicate sample identifiers in metadata.")
}

if (anyNA(crc_metadata$subject_id) ||
    any(!nzchar(crc_metadata$subject_id)) ||
    anyDuplicated(crc_metadata$subject_id)) {
  stop("Missing or repeated subject IDs; inspect participants before inference.")
}

missing_samples <- setdiff(
  crc_metadata$sample_id,
  colnames(abundance_species)
)

if (length(missing_samples) > 0L) {
  print(missing_samples)
  stop("Some requested metadata samples are missing from the local abundance file.")
}

# Explicitly align abundance columns and metadata rows.
abundance_species <- abundance_species[
  ,
  crc_metadata$sample_id,
  drop = FALSE
]

metadata <- as.data.frame(
  crc_metadata,
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$sample_id

stopifnot(
  identical(colnames(abundance_species), rownames(metadata))
)

cat(sprintf(
  "Locally loaded matrix: %d species x %d samples\n",
  nrow(abundance_species),
  ncol(abundance_species)
))
cat("Requested YuJ_2015 samples by outcome:\n")
print(table(metadata$study_condition, useNA = "ifany"))


# ============================================================================
# 1C. Preserve original species labels
# ============================================================================

original_feature_id <- rownames(abundance_species)

extract_original_species <- function(feature_id) {
  taxonomy_parts <- strsplit(
    feature_id,
    split = "|",
    fixed = TRUE
  )[[1]]
  
  species_part <- taxonomy_parts[
    grepl("^s__", taxonomy_parts)
  ]
  
  if (length(species_part) > 0L) {
    species_name <- species_part[length(species_part)]
  } else {
    species_name <- taxonomy_parts[length(taxonomy_parts)]
  }
  
  species_name <- sub("^s__", "", species_name)
  species_name <- gsub("_", " ", species_name, fixed = TRUE)
  
  species_name <- gsub(
    "\\bsp CAG ([0-9]+)\\b",
    "sp. CAG:\\1",
    species_name
  )
  
  species_name <- trimws(species_name)
  
  if (is.na(species_name) || !nzchar(species_name)) {
    species_name <- feature_id
  }
  
  species_name
}

species_label <- vapply(
  original_feature_id,
  extract_original_species,
  character(1)
)

duplicated_species <- duplicated(species_label) |
  duplicated(species_label, fromLast = TRUE)

if (any(duplicated_species)) {
  duplicate_table <- data.frame(
    Original_feature_id = original_feature_id[duplicated_species],
    Extracted_species = species_label[duplicated_species],
    stringsAsFactors = FALSE
  )
  
  print(duplicate_table, row.names = FALSE)
  stop(
    "Duplicated species labels were extracted. ",
    "Inspect duplicate_table before continuing."
  )
}

rownames(abundance_species) <- species_label

# ============================================================================
# 1D. Define the CRC/control outcome
# ============================================================================

condition_column_index <- match(
  "study_condition",
  tolower(colnames(metadata))
)

if (is.na(condition_column_index)) {
  stop("A study_condition column was not found in metadata.")
}

study_condition <- trimws(
  as.character(metadata[[condition_column_index]])
)

eligible_outcome <- study_condition %in% c("CRC", "control")

abundance_species <- abundance_species[
  ,
  eligible_outcome,
  drop = FALSE
]

metadata <- metadata[eligible_outcome, , drop = FALSE]
study_condition <- study_condition[eligible_outcome]

if (!all(c("CRC", "control") %in% study_condition)) {
  stop("Both CRC and control samples are required.")
}

Y <- as.numeric(study_condition == "CRC")

cat(sprintf("CRC samples: %d\n", sum(Y == 1)))
cat(sprintf("Control samples: %d\n", sum(Y == 0)))

# ============================================================================
# 4. Prepare the full composition and define eligible target species
# ============================================================================

X <- t(abundance_species)
if (any(!is.finite(X)) || any(X < 0)) {
  stop("Abundance values must be finite and nonnegative.")
}

valid_sample <- (
  is.finite(rowSums(X)) &
    rowSums(X) > 0
)

X <- X[
  valid_sample,
  ,
  drop = FALSE
]

Y <- Y[valid_sample]
metadata <- metadata[valid_sample, , drop = FALSE]

X <- X / rowSums(X)

# Apply the same numerical zero rule used by reparameterize().
X[X < 1e-9] <- 0

valid_sample_after_threshold <- rowSums(X) > 0

X <- X[
  valid_sample_after_threshold,
  ,
  drop = FALSE
]

Y <- Y[valid_sample_after_threshold]
metadata <- metadata[valid_sample_after_threshold, , drop = FALSE]
X <- X / rowSums(X)

if (length(unique(Y)) != 2L) {
  stop("The retained response does not contain both CRC and control samples.")
}

# Keep the original full composition, including any all-zero coordinates;
# the prevalence restriction below applies only to TARGETS.
stopifnot(identical(rownames(X), rownames(metadata)), nrow(X) == length(Y))
if (ncol(X) < 2L || any(rowSums(X >= 1e-9) < 2L)) {
  stop("At least two present components per sample are required for removal.")
}
if (min(table(Y)) < nfolds) stop("Too few cases or controls for the chosen folds.")

# Verify final analytic subjects before treating rows as independent observations.
subject_column <- match("subject_id", tolower(names(metadata)))
if (is.na(subject_column)) {
  stop("No subject_id field: verify participant identities before independent-row inference.")
}
subject_id <- trimws(as.character(metadata[[subject_column]]))
if (anyNA(subject_id) || any(!nzchar(subject_id))) {
  stop("Missing subject_id values: independence cannot be checked.")
}
if (anyDuplicated(subject_id)) {
  repeated_subjects <- unique(subject_id[duplicated(subject_id)])
  print(metadata[subject_id %in% repeated_subjects, , drop = FALSE])
  stop("Repeated subjects detected. Choose one observation per person before fitting.")
}
study_column <- match("study_name", tolower(names(metadata)))
if (!is.na(study_column) &&
    !all(as.character(metadata[[study_column]]) %in% "YuJ_2015")) {
  stop("The samples are not exclusively YuJ_2015.")
}
# YuJ_2015 identifies the Chinese discovery cohort, not the Danish validation set.
country_column <- match("country", tolower(names(metadata)))
if (!is.na(country_column)) {
  country <- trimws(as.character(metadata[[country_column]]))
  cat("Recorded countries:\n")
  print(table(country, useNA = "ifany"))
  if (anyNA(country) || any(!country %in% c("CHN", "China", "china"))) {
    stop("Unexpected country labels: verify that only the Chinese Yu cohort was loaded.")
  }
}
cat(sprintf("Unique subjects in final analysis: %d\n", length(unique(subject_id))))
if (nrow(X) != 128L || sum(Y == 1) != 74L || sum(Y == 0) != 54L) {
  warning("Counts differ from the published 74 CRC / 54 control cohort. Inspect metadata and preprocessing before interpreting results.")
}

species_taxonomy <- data.frame(
  Species = species_label, Original_feature_id = original_feature_id,
  stringsAsFactors = FALSE
)
n <- nrow(X)
p_total <- ncol(X)

presence_count <- colSums(X >= 1e-9)
absence_count <- n - presence_count
q_empirical <- presence_count / n

target_eligible <- (
  q_empirical >= q_lower &
    q_empirical <= q_upper
)

target_indices <- which(target_eligible)
p_target <- length(target_indices)

if (p_target == 0L) {
  stop("No species satisfy the target prevalence restriction.")
}

cat("\n============================================================\n")
cat("Analysis dimensions\n")
cat("============================================================\n")
cat(sprintf("Samples: %d\n", n))
cat(sprintf("Species retained in the full composition: %d\n", p_total))
cat(
  sprintf(
    "Eligible target species with %.2f <= q_hat <= %.2f: %d\n",
    q_lower,
    q_upper,
    p_target
  )
)
cat(
  sprintf(
    "Smallest present group among eligible targets: %d\n",
    min(presence_count[target_indices])
  )
)
cat(
  sprintf(
    "Smallest absent group among eligible targets: %d\n",
    min(absence_count[target_indices])
  )
)
cat(sprintf("CRC cases: %d\n", sum(Y == 1)))
cat(sprintf("Controls: %d\n", sum(Y == 0)))
cat(sprintf("CRC proportion: %.4f\n", mean(Y)))
cat(sprintf("Cross-fitting folds: %d\n", nfolds))
cat(sprintf("Parallel workers: %d\n", ncores))
cat(sprintf("All-zero coordinates retained in the full composition: %d\n",
            sum(colSums(X) == 0)))
cat("Output folder:", normalizePath(getwd(), winslash = "/"), "\n")

# 1C. q distribution over ALL coordinates, not only the eligible targets.
q_distribution <- data.frame(
  Dataset = "YuJ_2015_CRC_vs_control",
  Samples = n, Total_species = p_total,
  Group = 1:5,
  Interval = c("[0,0.1) or [0.9,1]", "[0.1,0.3)", "[0.3,0.5)",
               "[0.5,0.7)", "[0.7,0.9)"),
  Species_count = c(
    sum(q_empirical < 0.1 | q_empirical >= 0.9),
    sum(q_empirical >= 0.1 & q_empirical < 0.3),
    sum(q_empirical >= 0.3 & q_empirical < 0.5),
    sum(q_empirical >= 0.5 & q_empirical < 0.7),
    sum(q_empirical >= 0.7 & q_empirical < 0.9)
  )
)
q_distribution$Percent <- round(100 * q_distribution$Species_count / p_total, 2)
stopifnot(sum(q_distribution$Species_count) == p_total)
cat("\nPrevalence distribution over the full composition:\n")
print(q_distribution, row.names = FALSE)
flush.console()

check_rep <- reparameterize(X = X, s = target_indices[1])
if (!isTRUE(all.equal(
  as.numeric(check_rep$L), as.numeric(X[, target_indices[1]] < 1e-9)
))) {
  stop("reparameterize() must use L = 1 for absence, so q denotes presence.")
}


# 2. VIPER: fit ALL eligible targets using the unchanged external algorithm.
required_fit_columns <- c(
  "s_idx", "Species", "phi_hat", "sigma_phi_sq", "theta_hat",
  "q_hat", "present_n", "absent_n", "failed", "error_message"
)
cat("\nStarting parallel VIPER estimation for ALL eligible species...\n")
flush.console()
viper_start <- Sys.time()
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
results_list <- tryCatch({
  parallel::clusterExport(cl, "algorithm_file")
  parallel::clusterEvalQ(cl, {
    source(algorithm_file)
    NULL
  })
  foreach(
    s_idx = target_indices,
    .packages = "randomForest",
    .errorhandling = "stop"
  ) %dopar% {
    tryCatch({
      rep_X <- reparameterize(X = X, s = s_idx)
      fit_j <- R2_partial_linear(
        Y = Y, L = rep_X$L, W = rep_X$W, s = s_idx,
        nfolds = nfolds, seed = base_seed + s_idx
      )
      data.frame(
        s_idx = s_idx, Species = colnames(X)[s_idx],
        phi_hat = fit_j$phi_hat, sigma_phi_sq = fit_j$sigma_phi_sq,
        theta_hat = fit_j$theta_hat,
        q_hat = mean((1 - rep_X$L)^2),
        present_n = sum(rep_X$L == 0), absent_n = sum(rep_X$L == 1),
        failed = FALSE, error_message = "", stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        s_idx = s_idx, Species = colnames(X)[s_idx],
        phi_hat = NA_real_, sigma_phi_sq = NA_real_,
        theta_hat = NA_real_,
        q_hat = q_empirical[s_idx],
        present_n = presence_count[s_idx], absent_n = absence_count[s_idx],
        failed = TRUE, error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
  }
}, finally = {
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
})
fit_results <- do.call(rbind, results_list)
cat(sprintf("VIPER completed in %.2f minutes\n",
            as.numeric(difftime(Sys.time(), viper_start, units = "mins"))))

if (!all(required_fit_columns %in% names(fit_results)) ||
    nrow(fit_results) != p_target ||
    !identical(as.integer(fit_results$s_idx), as.integer(target_indices)) ||
    !identical(as.character(fit_results$Species), colnames(X)[target_indices])) {
  stop("VIPER results do not match the current eligible species; refit.")
}
if (any(fit_results$failed)) {
  cat("\nFailed VIPER fits (not counted as negative tests):\n")
  print(fit_results[fit_results$failed,
                    c("Species", "error_message"), drop = FALSE],
        row.names = FALSE)
}
# Crucially: DO NOT replace target_indices with the VIPER top 10.
# Every eligible target is evaluated by LC and RF below; top 10 is display-only.

# 4. LC/RF estimation and fully refitted paired bootstrap.
plugin_ridge_path <- function(y, x, lambda, basis) {
  x_mean <- colMeans(x)
  y_mean <- mean(y)
  z <- sweep(x, 2, x_mean, "-") %*% basis
  decomp <- svd(z)
  uy <- as.numeric(crossprod(decomp$u, y - y_mean))
  coef_path <- vapply(lambda, function(lam) {
    as.numeric(decomp$v %*% (decomp$d * uy /
                               (decomp$d^2 + lam)))
  }, numeric(ncol(basis)))
  coef_path <- matrix(coef_path, nrow = ncol(basis))
  beta <- basis %*% coef_path
  list(beta = beta, intercept = as.numeric(y_mean - x_mean %*% beta))
}

plugin_fit_lc <- function(y, x, lambda_grid, basis, original_ids) {
  # All bootstrap copies of one individual stay in the same INNER CV fold.
  ids <- unique(original_ids)
  inner_k <- min(3L, length(ids))
  if (inner_k < 2L) stop("Too few distinct training observations for LC CV.")
  id_fold <- sample(rep(seq_len(inner_k), length.out = length(ids)))
  cv_fold <- id_fold[match(original_ids, ids)]
  cv_error <- numeric(length(lambda_grid))
  for (v in seq_len(inner_k)) {
    train <- cv_fold != v
    path <- plugin_ridge_path(y[train], x[train, , drop = FALSE],
                              lambda_grid, basis)
    pred <- sweep(x[!train, , drop = FALSE] %*% path$beta,
                  2, path$intercept, "+")
    cv_error <- cv_error + colSums((pred - y[!train])^2)
  }
  lambda <- lambda_grid[which.min(cv_error)]
  fit <- plugin_ridge_path(y, x, lambda, basis)
  list(beta = as.numeric(fit$beta), intercept = fit$intercept[1], lambda = lambda)
}

plugin_compute <- function(X, y, targets, fold_id, original_ids, pseud, lambda_grid, ntree) {
  p <- ncol(X)
  n <- nrow(X)
  sy <- var(y)
  if (!is.finite(sy) || sy <= 0) stop("Nonpositive response variance.")
  # Orthonormal basis of {beta: sum(beta)=0}.
  basis <- stats::contr.helmert(p)
  basis <- sweep(basis, 2, sqrt(colSums(basis^2)), "/")
  log_X <- log(X + pseud)
  sq_lc <- sq_rf <- matrix(NA_real_, nrow = n, ncol = length(targets))
  for (k in sort(unique(fold_id))) {
    test <- which(fold_id == k)
    train <- which(fold_id != k)
    if (length(train) < 3L) stop("Too few training observations.")
    fit_lc <- plugin_fit_lc(y[train], log_X[train, , drop = FALSE],
                            lambda_grid, basis, original_ids[train])
    constant_y <- length(unique(y[train])) == 1L
    fit_rf <- if (!constant_y) randomForest::randomForest(
      x = X[train, , drop = FALSE], y = as.numeric(y[train]), ntree = ntree
    ) else NULL
    z_test <- X[test, , drop = FALSE]
    full_lc <- as.numeric(fit_lc$intercept +
                            log_X[test, , drop = FALSE] %*% fit_lc$beta)
    full_rf <- if (constant_y) rep(mean(y[train]), length(test)) else
      as.numeric(predict(fit_rf, z_test))
    for (j in seq_along(targets)) {
      W <- z_test
      W[, targets[j]] <- 0
      remaining <- rowSums(W)
      if (any(remaining <= 0)) {
        stop("A sample has only the target species present: its removal is undefined.")
      }
      W <- W / remaining
      perturbed_lc <- as.numeric(fit_lc$intercept + log(W + pseud) %*% fit_lc$beta)
      perturbed_rf <- if (constant_y) rep(mean(y[train]), length(test)) else
        as.numeric(predict(fit_rf, W))
      sq_lc[test, j] <- (full_lc - perturbed_lc)^2
      sq_rf[test, j] <- (full_rf - perturbed_rf)^2
    }
  }
  if (any(!is.finite(sq_lc)) || any(!is.finite(sq_rf))) {
    stop("Missing or nonfinite cross-fitted plug-in predictions.")
  }
  list(LC = colMeans(sq_lc) / sy, RF = colMeans(sq_rf) / sy)
}

if (plugin_n_boot < 2L) stop("plugin_n_boot must be at least 2.")
set.seed(base_seed + 300000L)
plugin_fold <- integer(n)
for (value in unique(Y)) {
  rows <- which(Y == value)
  plugin_fold[rows] <- sample(rep(seq_len(nfolds), length.out = length(rows)))
}
cat("\nLC/RF: fitting ALL eligible target species (not only VIPER top 10)...\n")
flush.console()
plugin_point <- plugin_compute(
  X, Y, target_indices, plugin_fold, seq_len(n), plugin_pseud,
  plugin_ridge_grid, plugin_rf_ntree
)
boot_lc <- boot_rf <- matrix(NA_real_, plugin_n_boot, p_target)
boot_messages <- rep("", plugin_n_boot)
cl_plugin <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl_plugin)
plugin_start <- Sys.time()
tryCatch({
  for (first in seq.int(1L, plugin_n_boot, by = ncores)) {
    batch <- seq.int(first, min(first + ncores - 1L, plugin_n_boot))
    cat(sprintf("LC/RF bootstrap: starting %d-%d of %d\n",
                min(batch), max(batch), plugin_n_boot))
    flush.console()
    boot_batch <- foreach(
      b = batch, .packages = "randomForest",
      .export = c("plugin_compute", "plugin_fit_lc", "plugin_ridge_path")
    ) %dopar% {
      tryCatch({
        set.seed(base_seed + 400000L + b)
        rows <- sample.int(n, n, replace = TRUE)
        # Copies inherit their original outer-fold IDs to prevent leakage
        # of an individual between training and validation sets.
        answer <- plugin_compute(
          X[rows, , drop = FALSE], Y[rows], target_indices,
          plugin_fold[rows], rows, plugin_pseud, plugin_ridge_grid,
          plugin_rf_ntree
        )
        list(LC = answer$LC, RF = answer$RF, error = "")
      }, error = function(e) {
        list(LC = rep(NA_real_, p_target), RF = rep(NA_real_, p_target),
             error = conditionMessage(e))
      })
    }
    for (i in seq_along(batch)) {
      boot_lc[batch[i], ] <- boot_batch[[i]]$LC
      boot_rf[batch[i], ] <- boot_batch[[i]]$RF
      boot_messages[batch[i]] <- boot_batch[[i]]$error
    }
    cat(sprintf("LC/RF bootstrap: completed %d/%d; %.2f minutes elapsed\n",
                max(batch), plugin_n_boot,
                as.numeric(difftime(Sys.time(), plugin_start, units = "mins"))))
    flush.console()
  }
}, finally = {
  parallel::stopCluster(cl_plugin)
  foreach::registerDoSEQ()
})


# 5. Common inference and result tables over ALL eligible targets.
# Point estimates are nonnegative. Negative Wald CI endpoints are not truncated.
# Two-sided CI uses z_(1-alpha/2); the one-sided threshold test uses z_(1-alpha).
make_result_table <- function(method, point, se,
                              valid_bootstrap = rep(NA_integer_, p_target)) {
  point <- as.numeric(point)
  se <- as.numeric(se)
  if (length(point) != p_target || length(se) != p_target) {
    stop("Result vector lengths do not match the eligible targets.")
  }
  valid_estimate <- is.finite(point) & point >= 0
  valid_inference <- valid_estimate & is.finite(se) & se > 0
  se[!valid_inference] <- NA_real_
  threshold_lower <- point - qnorm(1 - alpha) * se
  p_threshold <- pnorm((point - c_threshold) / se, lower.tail = FALSE)
  data.frame(
    s_idx = target_indices, Species = colnames(X)[target_indices],
    q_hat = unname(q_empirical[target_indices]),
    Method = method, estimate = point, se = se,
    ci_lower = point - qnorm(1 - alpha / 2) * se,
    ci_upper = point + qnorm(1 - alpha / 2) * se,
    threshold_lower = threshold_lower, p_threshold = p_threshold,
    selected_threshold = ifelse(valid_inference, p_threshold < alpha, NA),
    point_above_threshold = ifelse(valid_estimate, point > c_threshold, NA),
    valid_estimate = valid_estimate, valid_inference = valid_inference,
    valid_bootstrap = valid_bootstrap, stringsAsFactors = FALSE
  )
}

viper_point <- fit_results$phi_hat
viper_point[fit_results$failed] <- NA_real_
viper_se <- sqrt(pmax(fit_results$sigma_phi_sq, 0) / n)
viper_se[fit_results$failed | !is.finite(fit_results$sigma_phi_sq) |
           fit_results$sigma_phi_sq < 0] <- NA_real_
viper_results <- make_result_table("VIPER", viper_point, viper_se)

plugin_summarize <- function(point, boot, method) {
  if (ncol(boot) != p_target) stop(method, " bootstrap column count mismatch.")
  valid <- colSums(is.finite(boot))
  # Bootstrap SD is already an estimator SE, not an influence-function SD.
  se <- apply(boot, 2, function(x) sd(x[is.finite(x)]))
  enough <- valid >= max(2L, ceiling(0.8 * nrow(boot)))
  se[!enough | !is.finite(se) | se <= 0] <- NA_real_
  if (any(!enough)) {
    warning(method, ": too few valid bootstrap draws for ", sum(!enough),
            " species; their inference is unavailable.")
  }
  make_result_table(method, point, se, valid_bootstrap = valid)
}
comparison_results <- rbind(
  viper_results,
  plugin_summarize(plugin_point$LC, boot_lc, "LC"),
  plugin_summarize(plugin_point$RF, boot_rf, "RF")
)
rownames(comparison_results) <- NULL
stopifnot(nrow(comparison_results) == 3L * p_target)

if (any(nzchar(boot_messages))) {
  cat("\nLC/RF bootstrap failures:\n")
  print(table(boot_messages[nzchar(boot_messages)]))
}

# Retain useful VIPER-specific aliases for later inspection.
ranking_results <- fit_results
ranking_results$se_phi <- viper_results$se
ranking_results$phi_ci_lower <- viper_results$ci_lower
ranking_results$phi_ci_upper <- viper_results$ci_upper
ranking_results$threshold_ci_lower <- viper_results$threshold_lower
ranking_results$p_value_phi_threshold <- viper_results$p_threshold
ranking_results$selected_threshold <- viper_results$selected_threshold
ranking_results <- ranking_results[
  order(ranking_results$phi_hat, decreasing = TRUE, na.last = TRUE), , drop = FALSE
]
ranking_results$VIPER_rank <- NA_integer_
finite_rank <- which(is.finite(ranking_results$phi_hat) & !ranking_results$failed)
ranking_results$VIPER_rank[finite_rank] <- seq_along(finite_rank)

# 6. Print BOTH requested selections for ALL eligible targets.
methods <- c("VIPER", "LC", "RF")
method_summary <- do.call(rbind, lapply(methods, function(method_now) {
  temp <- comparison_results[comparison_results$Method == method_now, , drop = FALSE]
  data.frame(
    Method = method_now, Eligible_targets = nrow(temp),
    Valid_estimates = sum(temp$valid_estimate),
    Valid_inference = sum(temp$valid_inference),
    Significant_above_c = sum(temp$selected_threshold %in% TRUE),
    Point_estimates_above_c = sum(temp$point_above_threshold %in% TRUE)
  )
}))
cat(sprintf("\nALL eligible species: c = %.3f, pointwise alpha = %.3f\n",
            c_threshold, alpha))
print(method_summary, row.names = FALSE)
cat("Unavailable inference is not counted as a negative test.\n")

selected_signals <- comparison_results[
  comparison_results$selected_threshold %in% TRUE, , drop = FALSE
]
above_threshold <- comparison_results[
  comparison_results$point_above_threshold %in% TRUE, , drop = FALSE
]
for (method_now in methods) {
  selected_now <- selected_signals[
    selected_signals$Method == method_now, , drop = FALSE
  ]
  selected_now <- selected_now[
    order(selected_now$estimate, decreasing = TRUE), , drop = FALSE
  ]
  above_now <- above_threshold[
    above_threshold$Method == method_now, , drop = FALSE
  ]
  above_now <- above_now[
    order(above_now$estimate, decreasing = TRUE), , drop = FALSE
  ]
  cat(sprintf("\n%s: %d species SIGNIFICANTLY above c = %.3f\n",
              method_now, nrow(selected_now), c_threshold))
  print(selected_now[, c("Species", "q_hat", "estimate", "se",
                         "ci_lower", "ci_upper", "threshold_lower",
                         "p_threshold"), drop = FALSE],
        row.names = FALSE, digits = 5)
  cat(sprintf("\n%s: %d species with POINT ESTIMATES above c = %.3f\n",
              method_now, nrow(above_now), c_threshold))
  print(above_now[, c("Species", "q_hat", "estimate"), drop = FALSE],
        row.names = FALSE, digits = 5)
}

# Additional descriptive comparison: VIPER-selected species not selected by LC.
# "Not selected" is not evidence of absence; LC inference must be available.
lc_valid_species <- comparison_results$Species[
  comparison_results$Method == "LC" & comparison_results$valid_inference
]
lc_selected_species <- selected_signals$Species[selected_signals$Method == "LC"]
VIPER_only_vs_LC <- selected_signals[
  selected_signals$Method == "VIPER" &
    selected_signals$Species %in% lc_valid_species &
    !selected_signals$Species %in% lc_selected_species,
  , drop = FALSE
]
cat(sprintf("\nVIPER significant but LC not significant: %d species\n",
            nrow(VIPER_only_vs_LC)))
print(VIPER_only_vs_LC[, c("Species", "q_hat", "estimate", "se",
                           "threshold_lower", "p_threshold"), drop = FALSE],
      row.names = FALSE, digits = 5)

# 7. Same three-panel figure; the VIPER top 10 are used ONLY for display.
vip_top <- viper_results[viper_results$valid_estimate, , drop = FALSE]
vip_top <- head(
  vip_top[order(vip_top$estimate, decreasing = TRUE), , drop = FALSE], top_k
)
if (nrow(vip_top) == 0L) stop("No valid VIPER point estimates for plotting.")
plot_results <- comparison_results[
  comparison_results$Species %in% vip_top$Species, , drop = FALSE
]
plot_results$Species <- factor(plot_results$Species, levels = rev(vip_top$Species))
plot_results$Method <- factor(plot_results$Method, levels = methods)
p_same_target <- ggplot(plot_results, aes(x = estimate, y = Species)) +
  geom_vline(xintercept = c_threshold, linetype = "dashed", color = "gray45") +
  geom_segment(
    aes(x = ci_lower, xend = ci_upper, yend = Species),
    color = "gray45", na.rm = TRUE
  ) +
  geom_point(aes(color = selected_threshold), size = 2.4, na.rm = TRUE) +
  facet_wrap(~Method, nrow = 1, scales = "free_x") +
  scale_color_manual(
    values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"),
    na.value = "gray60",
    guide = "none"
  ) +
  labs(x = "Perturbation importance estimate", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black")
  )
if (interactive()) print(p_same_target)
ggsave("YuJ_2015_CRC_species_VIPER_LC_RF_comparison.png",
       p_same_target, width = 15, height = 7, dpi = 300, bg = "white")
cat("\nSaved YuJ_2015_CRC_species_VIPER_LC_RF_comparison.png\n")
# Inference is pointwise, under model assumptions; no multiple-testing adjustment.
# Unique subject IDs prevent repeated-record analysis but do not by themselves prove i.i.d.
# The top-10 figure is exploratory, not post-selection adjusted.


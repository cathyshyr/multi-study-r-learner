# Scenario B

# Load packages and data
source("Utils.R")
load("edat_orig.RData")

# Remove small studies with < 150 observations
edat_orig <- edat_orig[which(sapply(edat_orig, function(x) nrow(x)) >= 150)]

# Number of studies
K <- 4

# Total number of observations
n <- sum(sapply(edat_orig, nrow))

# Total number of simulation runs
nreps <- 1000

# Number of covariates
p <- 5

covariate_cols <- sapply(1:p, function(x) paste("X.", x, sep = ""))

# Studies
set.seed(1)
edat_source <- init_data(edat_orig[1:(K - 1)], ndat = K - 1, nvar = p)
for (i in 1:length(edat_source)) {
  dataset = edat_source[[i]]
  # dataset = dataset[1:n_k,]
  edat_source[[i]] = cbind(dataset[, 1:p])
}
edat_train <- do.call("rbind", edat_source[1:6])
edat_all <- do.call("rbind", edat_source[1:6])

X <- edat_all[, covariate_cols]
mu <- colMeans(X)
S_inv <- solve(cov(X))

# Compute Mahalanobis distance from global center
mahal_dist <- apply(X, 1, function(x) {
  sqrt(t(x - mu) %*% S_inv %*% (x - mu))
})

# Choose percentiles of distance
train_cutoff <- quantile(mahal_dist, 0.85)
test_cutoff  <- quantile(mahal_dist, 0.90)  

# Indices
train_idx <- which(mahal_dist <= train_cutoff)
test_idx  <- which(mahal_dist >= test_cutoff)

# Subset data
edat_train <- edat_all[train_idx, ]
edat_test  <- edat_all[test_idx, ]

# Create a grid of high-support test set points to evaluate on
qs_test <- apply(edat_test[, covariate_cols], 2, quantile, probs = c(0.05, 0.95))
grid_list <- lapply(1:p, function(j) seq(qs_test[1, j], qs_test[2, j], length.out = 5))
n_sample_grid <- 500
set.seed(3)
sampled_idx_list <- lapply(grid_list, function(vals) sample(vals, n_sample_grid, replace = TRUE))
edat_test <- data.frame(matrix(unlist(sampled_idx_list), ncol = p))
names(edat_test) <- names(edat_train)
nrow(edat_train)
nrow(edat_test)

######################################################################################
# Simulation function for calculating the MSE of conditional average treatment effect
# 
# Input: 
# n - total number of observations
# p - number of predictors
# K - number of source studies
# data - all studies combined
# covariate_cols - names of predictors
# cols_re_tau - column indices of the beta_tau vector with random effects
# cols_re_mu0 - column indices of the beta_mu0 vector with random effects
# cols_re_a - column indices of the beta_a vector with random effects
# beta_tau - coefficients for the heterogeneous treatment effect
# beta_mu0 - coefficients for the mean outcome model among untreated individuals
# beta_a - coefficients for the propensity score model
# beta_mu0_k - study-specific coefficients for the mean outcome model among untreated individuals
# beta_a_k - study-specific coefficients for the propensity score model
# beta_tau_k - study-specific coefficients for the heterogeneous treatment effect
# p_mat - n by K matrix of membership probabilities
# perturb_tau - variance of the random effects for beta_tau
# perturb_mu0 - variance of the random effects for beta_mu0
# perturb_a - variance of the random effects for beta_a
# k-folds - number of folds to perform cross-fitting 
# alpha - glmnet parameter (alpha = 1 => ridge, alpha = 0 => LASSO)
# iteration - iteration number
#
# Output:
# Data frame of the performance results (MSE, bias, variance) for the multi-study and study-specific R-learners
# across the oracle and plug-in settings
######################################################################################
sim_each <- function(n, p, K, data, covariate_cols, cols_re_tau, cols_re_mu0, cols_re_a, beta_tau, beta_mu0, beta_a, p_mat, perturb_tau, perturb_mu0, perturb_a, perturb_offset_tau, k_folds, alpha, iteration){
  check <- rep(0, K - 1)
  while(any(check == 0)){
    # Test set indices 
    ind_test <- which(data$S == K)
    
    # Initialize vectors to store coefficients
    X_int <- as.matrix(cbind(1, data[, covariate_cols]))
    
    beta_a_k <- beta_mu0_k <- vector("list", length = K - 1)
    # Compute study-specific tau_k, mu0_k and e_k
    for (k in 1:(K - 1)) {
      # mu0 
      sigma_mu0 <- rep(0, p + 1)
      sigma_mu0[cols_re_mu0[[k]]] <- perturb_mu0
      beta_mu0_k[[k]] <- beta_mu0
      beta_mu0_k[[k]] <- beta_mu0_k[[k]] + rnorm(p + 1, 0, sigma_mu0)
      
      # e
      sigma_a <- rep(0, p + 1)
      sigma_a[cols_re_a[[k]]] <- perturb_a
      beta_a_k[[k]] <- beta_a
      beta_a_k[[k]] <- beta_a_k[[k]] + rnorm(p + 1, 0, sigma_a)
    }
    
    # A|X, S=k is 0/1 Bernoulli distributed with probability e_k(X)
    true_e <- vector("list", length = K - 1)
    for(k in 1:(K - 1)){
      X_int_k <- as.matrix(cbind(1, data[which(data$S == k), covariate_cols]))
      true_e[[k]] <- expit(X_int_k %*% beta_a_k[[k]])
      data$A[which(data$S == k)] <- rbinom(nrow(data[which(data$S == k), ]), 1, true_e[[k]]) 
      check[k] <- ifelse(length(which(data$A[which(data$S == k)] == 1)) < 10 | length(which(data$A[which(data$S == k)] == 0)) < 10, 0, 1)
    }
  }
  
  data_train <- data[-ind_test, ]
  tau_k_list <- tau_list <- vector("list", length = K - 1)
  for (k in 1:(K - 1)) {
    offset_k <- perturb_offset_tau[[k]]
    tau_k_list[[k]] <- apply(data_train[, covariate_cols], 1, function(row) generate_tau_k(row, k, offset_k))
    tau_list[[k]] <- apply(data[, covariate_cols], 1, function(row) generate_tau_k(row, k, offset_k))
  }
  
  # Generate Y
  X_int_train <- X_int[-ind_test, ]
  data_train <- data[-ind_test, ]
  Y_mat <- sapply(1:(K - 1), function(k){
    (X_int_train %*% beta_mu0_k[[k]] + data_train$A * tau_k_list[[k]]) * p_mat[-ind_test, k]
  })
  data_train$Y <- apply(Y_mat, 1, sum) + rnorm(nrow(data_train), 0, 1)
  
  # Oracle m
  m_mat <- sapply(1:(K - 1), function(k){
    (X_int_train %*% beta_mu0_k[[k]] + expit(X_int_train %*% beta_a_k[[k]]) * tau_k_list[[k]]) * p_mat[-ind_test, k]
  })
  true_m <- apply(m_mat, 1, sum)
  
  # Estimate m (merge)
  m_estimate_merge <- rlearner_m_hat_merge(data = data_train, K = K - 1, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # Estimate a 
  a_estimate_merge <- rlearner_a_hat_merge(data = data_train, K = K - 1, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # Oracle tau
  true_tau <- rowSums(sapply(1:(K - 1), function(k) {
    tau_list[[k]] * p_mat[, k]
  }))
  
  true_tau <- true_tau - mean(true_tau)
  true_tau <- true_tau[ind_test]
  
  X_test <- data[ind_test, covariate_cols]
  X_test_spline <- generate_spline_basis(data[ind_test, covariate_cols])
  X_test_spline_int <- cbind(1, X_test_spline)
  covariate_cols_spline <- colnames(X_test_spline)
  p_spline <- ncol(X_test_spline)
  
  # Calculate additional metrics
  B <- 1000
  tau_boot_mat_rlearner_oracle <- tau_boot_mat_rlearner_merge <- tau_boot_mat_rlearner_oracle_SS <- tau_boot_mat_ms_rlearner_oracle <- tau_boot_mat_ms_rlearner_mMerge_aMerge <- tau_boot_mat_ms_rlearner_mSS_aSS <- tau_boot_mat_ms_rlearner_mMerge_aSS <- tau_boot_mat_ms_rlearner_mAvg_aSS <- tau_boot_mat_ms_rlearner_mStack_aSS <- tau_boot_mat_rlearner_estimate_SS <- matrix(NA, nrow = B, ncol = length(ind_test))
  for (b in 1:B) {
    X_k_b <- X_k_b_int <- X_k_b_spline <- X_k_b_int_spline <- Y_k_b <- W_k_b <- data_train_k_b <- data_train_k_b_spline <- p_mat_k_b <- true_e_k_b <- true_m_k_b <- vector("list", length = K - 1)
    # Sample with replacement from the training data
    for(k in 1:(K - 1)){
      train_idx_k <- which(data_train$S == k)  # training rows for study k
      idx <- sample(train_idx_k, replace = TRUE)
      X_k_b[[k]] <- data_train[idx, covariate_cols]
      X_k_b_spline[[k]] <- generate_spline_basis(X_k_b[[k]])
      X_k_b_int[[k]] <- as.matrix(cbind(1, X_k_b[[k]]))
      X_k_b_int_spline[[k]] <- as.matrix(cbind(1, X_k_b_spline[[k]]))
      Y_k_b[[k]] <- data_train[idx, ]$Y
      W_k_b[[k]] <- data_train[idx, ]$A
      p_mat_k_b[[k]] <- p_mat[idx, ]
      true_e_k_b[[k]] <- as.numeric(expit(X_k_b_int[[k]] %*% beta_a_k[[k]]))
      data_train_k_b[[k]] <- cbind(X_k_b[[k]], S = k, A = W_k_b[[k]], Y = Y_k_b[[k]])
      data_train_k_b_spline[[k]] <- cbind(X_k_b_spline[[k]], S = k, A = W_k_b[[k]], Y = Y_k_b[[k]])
    }
    # Combine to form training dataset
    X_b <- do.call("rbind", X_k_b)
    X_b_spline <- do.call("rbind", X_k_b_spline)
    Y_b <- unlist(Y_k_b)
    W_b <- unlist(W_k_b)
    data_train_b <- as.data.frame(do.call("rbind", data_train_k_b))
    data_train_b_spline <- as.data.frame(do.call("rbind", data_train_k_b_spline))
    X_int_b <- do.call("rbind", X_k_b_int)
    X_int_b_spline <- do.call("rbind", X_k_b_int_spline)
    tau_k_list_b <- vector("list")
    for(k in 1:(K - 1)){
      offset_k <- perturb_offset_tau[[k]]
      tau_k_list_b[[k]] <- apply(data_train_b[, covariate_cols], 1, function(row) generate_tau_k(row, k, offset_k))
    }
    p_mat_b <- do.call("rbind", p_mat_k_b)
    m_mat_b <- sapply(1:(K - 1), function(k){
      (X_int_b %*% beta_mu0_k[[k]] + expit(X_int_b %*% beta_a_k[[k]]) * tau_k_list_b[[k]]) * p_mat_b[, k]
    })
    true_m_b <- apply(m_mat_b, 1, sum)
    
    # Estimate m (merge)
    m_estimate_merge_b <- rlearner_m_hat_merge(data = data_train_b_spline, K = K - 1, p = p_spline, covariate_cols = covariate_cols_spline, k_folds = k_folds, alpha = alpha)
    
    # Estimate a 
    a_estimate_merge_b <- rlearner_a_hat_merge(data = data_train_b_spline, K = K - 1, p = p_spline, covariate_cols = covariate_cols_spline, k_folds = k_folds, alpha = alpha)
    
    ################################################
    #########   MS R-learner Merge Both    #########
    ################################################
    # Estimate ascertainment probability on training data
    rf_pmat_model_b <- randomForest(
      x = X_b,
      y = as.factor(data_train_b$S),
      ntree = 500,
      mtry = floor(sqrt(ncol(X_b))),
      probability = TRUE
    )
    
    p_pred_mat_train_b <- predict(rf_pmat_model_b, newdata = X_b, type = "prob")
    p_pred_mat_test_b <- predict(rf_pmat_model_b, newdata = X_test, type = "prob")
    
    X_res_mat_merge_b <- lapply(1:(K - 1), function(k){
      sweep(X_int_b_spline, 1, (data_train_b$A - a_estimate_merge_b) * p_pred_mat_train_b[, k], "*")
    })
    
    X_res_mat_merge_all_b <- as.data.frame(do.call("cbind", X_res_mat_merge_b))
    names(X_res_mat_merge_all_b) <- sapply(1:ncol(X_res_mat_merge_all_b), function(x) paste("X.", x, sep = ""))
    
    Y_res_merge_estimate_b <- data_train_b$Y - m_estimate_merge_b
    tau_boot_mat_ms_rlearner_mMerge_aMerge[b, ] <- ms_rlearner(X = X_res_mat_merge_all_b, 
                                                               Y = Y_res_merge_estimate_b, 
                                                               alpha = alpha, K = K - 1, p = p_spline, 
                                                               p_mat = p_pred_mat_test_b,
                                                               newX = X_test_spline_int)
    
    ################################################
    #########  MS R-learner Study-Specific #########
    ################################################
    m_estimate_study_specific_b <- rlearner_m_hat_stack(data = data_train_b_spline, K = K - 1, p = p_spline, covariate_cols = covariate_cols_spline, k_folds = k_folds, alpha = alpha, p_mat = p_pred_mat_train_b)
    m_estimate_avg_b <- m_estimate_study_specific_b$avg
    m_estimate_stack_b <- m_estimate_study_specific_b$stack_int
    m_estimate_ss_b <- m_estimate_study_specific_b$ss
    
    # Propensity score model e_k(x)
    a_hat_k_b <- m_hat_k_b <- vector("list", length = K - 1)
    a_hat_k_lambda_min_b <- m_hat_k_lambda_min_b <- vector()
    for(k in 1:(K - 1)){
      a_hat_k_b[[k]] <- glmnet::cv.glmnet(x = as.matrix(X_k_b_spline[[k]]),
                                          y = as.matrix(data_train_b$A[which(data_train_b$S == k)]),
                                          alpha = alpha,
                                          lambda = NULL,
                                          standardize = TRUE,
                                          keep = TRUE,
                                          family = "binomial",
                                          type.measure = "deviance")
      a_hat_k_lambda_min_b[k] = a_hat_k_b[[k]]$lambda[which.min(a_hat_k_b[[k]]$cvm[!is.na(colSums(a_hat_k_b[[k]]$fit.preval))])]
      
      m_hat_k_b[[k]] <- glmnet::cv.glmnet(x = as.matrix(X_k_b_spline[[k]]),
                                          y = as.matrix(data_train_b$Y[which(data_train_b$S == k)]),
                                          alpha = alpha,
                                          lambda = NULL,
                                          keep = TRUE)
      m_hat_k_lambda_min_b[k] = m_hat_k_b[[k]]$lambda[which.min(m_hat_k_b[[k]]$cvm[!is.na(colSums(m_hat_k_b[[k]]$fit.preval))])]
    }
    
    X_res_mat_avg_b <- lapply(1:(K - 1), function(k){
      sweep(X_int_b_spline, 1, (data_train_b$A - predict(a_hat_k_b[[k]], newx = as.matrix(X_b_spline), s = a_hat_k_lambda_min_b[k], type = "response")) * p_pred_mat_train_b[, k], "*")
    })
    
    X_res_mat_avg_all_b <- as.data.frame(do.call("cbind", X_res_mat_avg_b))
    names(X_res_mat_avg_all_b) <- sapply(1:ncol(X_res_mat_avg_all_b), function(x) paste("X.", x, sep = ""))
    
    # Merge m, Study-specific A
    tau_boot_mat_ms_rlearner_mMerge_aSS[b, ] <- ms_rlearner(X = X_res_mat_avg_all_b, Y = Y_res_merge_estimate_b, alpha = alpha, K = K - 1, p = p_spline, p_mat = p_pred_mat_test_b,
                                                            newX = X_test_spline_int)
    
    # Avg m, Study-specific A
    Y_res_avg_estimate_b <- data_train_b$Y - m_estimate_avg_b
    tau_boot_mat_ms_rlearner_mAvg_aSS[b, ] <- ms_rlearner(X = X_res_mat_avg_all_b, Y = Y_res_avg_estimate_b, alpha = alpha, K = K - 1, p = p_spline, p_mat = p_pred_mat_test_b,
                                                          newX = X_test_spline_int)
    
    # Stack m, Study-specific A
    Y_res_stack_estimate_b <- data_train_b$Y - m_estimate_stack_b
    tau_boot_mat_ms_rlearner_mStack_aSS[b, ] <- ms_rlearner(X = X_res_mat_avg_all_b, Y = Y_res_stack_estimate_b, alpha = alpha, K = K - 1, p = p_spline, p_mat = p_pred_mat_test_b,
                                                            newX = X_test_spline_int)
    
    # Study-specific m, Study-specific A
    Y_res_ss_estimate_b <- data_train_b$Y - m_estimate_ss_b
    tau_boot_mat_ms_rlearner_mSS_aSS[b, ] <- ms_rlearner(X = X_res_mat_avg_all_b, Y = Y_res_ss_estimate_b, alpha = alpha, K = K - 1, p = p_spline, p_mat = p_pred_mat_test_b,
                                                         newX = X_test_spline_int)
    
    
    ########################################
    ######### R-learner Oracle SS ##########
    ########################################
    
    ##########################################
    ######### R-learner Estimate SS ##########
    ##########################################
    rlearner_oracle_SS_b <- rlearner_oracle_SS_pred_b <- rlearner_estimate_SS_b <- rlearner_estimate_SS_pred_b <- vector("list", length = K - 1)
    for(k in 1:(K - 1)){
      rlearner_oracle_SS_b[[k]] <- rlasso(x = as.matrix(X_k_b_spline[[k]]),
                                          w = as.matrix(data_train_b$A[which(data_train_b$S == k)]),
                                          y = as.matrix(data_train_b$Y[which(data_train_b$S == k)]),
                                          lambda_choice = "lambda.min",
                                          p_hat = true_e_k_b[[k]],
                                          m_hat = true_m_k_b[[k]],
                                          alpha = alpha)
      rlearner_oracle_SS_pred_b[[k]] <- predict(rlearner_oracle_SS_b[[k]], newx = as.matrix(X_test_spline))
      
      rlearner_estimate_SS_b[[k]] <- rlasso(x = as.matrix(X_k_b_spline[[k]]),
                                            w = as.matrix(data_train_b$A[which(data_train_b$S == k)]),
                                            y = as.matrix(data_train_b$Y[which(data_train_b$S == k)]),
                                            lambda_choice = "lambda.min",
                                            p_hat = predict(a_hat_k_b[[k]], newx = as.matrix(X_k_b_spline[[k]]), s = a_hat_k_lambda_min_b[k], type = "response"),
                                            m_hat = predict(m_hat_k_b[[k]], newx = as.matrix(X_k_b_spline[[k]]), s = m_hat_k_lambda_min_b[k], type = "response"),
                                            alpha = alpha)
      rlearner_estimate_SS_pred_b[[k]] <- predict(rlearner_estimate_SS_b[[k]], newx = as.matrix(X_test_spline))
    }
    
    tau_boot_mat_rlearner_oracle_SS[b, ] <- apply(do.call("cbind", rlearner_oracle_SS_pred_b) * p_mat[ind_test, ], 1, sum)
    tau_boot_mat_rlearner_estimate_SS[b, ] <- apply(do.call("cbind", rlearner_estimate_SS_pred_b) * p_pred_mat_test_b, 1, sum)
    
    # ################################################
    # #########   MS R-learner Oracle Both  #########
    # ################################################
    # Training set predictors
    X_res_mat_with_int_b <- lapply(1:(K - 1), function(k){
      sweep(X_int_b_spline, 1, (data_train_b$A - expit(X_int_b %*% beta_a_k[[k]])) * p_mat_b[, k], "*")
    })
    
    X_res_all_oracle_b <- as.data.frame(do.call("cbind", X_res_mat_with_int_b))
    names(X_res_all_oracle_b) <- sapply(1:ncol(X_res_all_oracle_b), function(x) paste("X.", x, sep = ""))
    
    # Training set outcomes
    Y_res_merge_oracle_b <- data_train_b$Y - true_m_b
    
    # Fit multi-study R-learner
    tau_boot_mat_ms_rlearner_oracle[b, ] <- ms_rlearner(X = X_res_all_oracle_b, Y = Y_res_merge_oracle_b, alpha = alpha, K = K - 1, p = p_spline, p_mat = p_mat[ind_test, ],
                                                        newX = X_test_spline_int)
    
    
  }
  
  # #######################################
  # ######### R-learner Oracle SS #########
  # #######################################
  se_tau_rlearner_oracle_SS <- apply(tau_boot_mat_rlearner_oracle_SS, 2, sd)
  ci_lower_rlearner_oracle_SS <- apply(tau_boot_mat_rlearner_oracle_SS, 2, quantile, probs = 0.025)
  ci_upper_rlearner_oracle_SS <- apply(tau_boot_mat_rlearner_oracle_SS, 2, quantile, probs = 0.975)
  coverage_rlearner_oracle_SS <- mean((true_tau >= ci_lower_rlearner_oracle_SS) & (true_tau <= ci_upper_rlearner_oracle_SS))
  type1_error_rlearner_oracle_SS <- mean((true_tau < ci_lower_rlearner_oracle_SS | true_tau > ci_upper_rlearner_oracle_SS))
  ci_length_rlearner_oracle_SS <- mean(ci_upper_rlearner_oracle_SS - ci_lower_rlearner_oracle_SS)
  tau_hat_mean_rlearner_oracle_SS <- colMeans(tau_boot_mat_rlearner_oracle_SS)  # length n_test
  bias_squared_pointwise_rlearner_oracle_SS <- (tau_hat_mean_rlearner_oracle_SS - true_tau)^2
  variance_pointwise_rlearner_oracle_SS <- apply(tau_boot_mat_rlearner_oracle_SS, 2, var)
  mse_pointwise_rlearner_oracle_SS <- bias_squared_pointwise_rlearner_oracle_SS + variance_pointwise_rlearner_oracle_SS
  bias_squared_rlearner_oracle_SS <- mean(bias_squared_pointwise_rlearner_oracle_SS)
  variance_hat_tau_rlearner_oracle_SS <- mean(variance_pointwise_rlearner_oracle_SS)
  mse_rlearner_oracle_SS <- mean(mse_pointwise_rlearner_oracle_SS)
  
  #######################################
  ######### R-learner Estimate SS #######
  #######################################
  se_tau_rlearner_estimate_SS <- apply(tau_boot_mat_rlearner_estimate_SS, 2, sd)
  ci_lower_rlearner_estimate_SS <- apply(tau_boot_mat_rlearner_estimate_SS, 2, quantile, probs = 0.025)
  ci_upper_rlearner_estimate_SS <- apply(tau_boot_mat_rlearner_estimate_SS, 2, quantile, probs = 0.975)
  coverage_rlearner_estimate_SS <- mean((true_tau >= ci_lower_rlearner_estimate_SS) & (true_tau <= ci_upper_rlearner_estimate_SS))
  type1_error_rlearner_estimate_SS <- mean((true_tau < ci_lower_rlearner_estimate_SS | true_tau > ci_upper_rlearner_estimate_SS))
  ci_length_rlearner_estimate_SS <- mean(ci_upper_rlearner_estimate_SS - ci_lower_rlearner_estimate_SS)
  tau_hat_mean_rlearner_estimate_SS <- colMeans(tau_boot_mat_rlearner_estimate_SS)  # length n_test
  bias_squared_pointwise_rlearner_estimate_SS <- (tau_hat_mean_rlearner_estimate_SS - true_tau)^2
  variance_pointwise_rlearner_estimate_SS <- apply(tau_boot_mat_rlearner_estimate_SS, 2, var)
  mse_pointwise_rlearner_estimate_SS <- bias_squared_pointwise_rlearner_estimate_SS + variance_pointwise_rlearner_estimate_SS
  bias_squared_rlearner_estimate_SS <- mean(bias_squared_pointwise_rlearner_estimate_SS)
  variance_hat_tau_rlearner_estimate_SS <- mean(variance_pointwise_rlearner_estimate_SS)
  mse_rlearner_estimate_SS <- mean(mse_pointwise_rlearner_estimate_SS)
  
  # ############################################
  # ######### MS R-learner Oracle Both #########
  # ############################################
  se_tau_ms_rlearner_oracle <- apply(tau_boot_mat_ms_rlearner_oracle, 2, sd)
  ci_lower_ms_rlearner_oracle <- apply(tau_boot_mat_ms_rlearner_oracle, 2, quantile, probs = 0.025)
  ci_upper_ms_rlearner_oracle <- apply(tau_boot_mat_ms_rlearner_oracle, 2, quantile, probs = 0.975)
  coverage_ms_rlearner_oracle <- mean((true_tau >= ci_lower_ms_rlearner_oracle) & (true_tau <= ci_upper_ms_rlearner_oracle))
  type1_error_ms_rlearner_oracle <- mean((true_tau < ci_lower_ms_rlearner_oracle | true_tau > ci_upper_ms_rlearner_oracle))
  ci_length_ms_rlearner_oracle <- mean(ci_upper_ms_rlearner_oracle - ci_lower_ms_rlearner_oracle)
  tau_hat_mean_ms_rlearner_oracle <- colMeans(tau_boot_mat_ms_rlearner_oracle)  # length n_test
  bias_squared_pointwise_ms_rlearner_oracle <- (tau_hat_mean_ms_rlearner_oracle - true_tau)^2
  variance_pointwise_ms_rlearner_oracle <- apply(tau_boot_mat_ms_rlearner_oracle, 2, var)
  mse_pointwise_ms_rlearner_oracle <- bias_squared_pointwise_ms_rlearner_oracle + variance_pointwise_ms_rlearner_oracle
  bias_squared_ms_rlearner_oracle <- mean(bias_squared_pointwise_ms_rlearner_oracle)
  variance_hat_tau_ms_rlearner_oracle <- mean(variance_pointwise_ms_rlearner_oracle)
  mse_ms_rlearner_oracle <- mean(mse_pointwise_ms_rlearner_oracle)
  
  ############################################
  ######### MS R-learner Merge Both  #########
  ############################################
  se_tau_ms_rlearner_mMerge_aMerge <- apply(tau_boot_mat_ms_rlearner_mMerge_aMerge, 2, sd)
  ci_lower_ms_rlearner_mMerge_aMerge <- apply(tau_boot_mat_ms_rlearner_mMerge_aMerge, 2, quantile, probs = 0.025)
  ci_upper_ms_rlearner_mMerge_aMerge <- apply(tau_boot_mat_ms_rlearner_mMerge_aMerge, 2, quantile, probs = 0.975)
  coverage_ms_rlearner_mMerge_aMerge <- mean((true_tau >= ci_lower_ms_rlearner_mMerge_aMerge) & (true_tau <= ci_upper_ms_rlearner_mMerge_aMerge))
  type1_error_ms_rlearner_mMerge_aMerge <- mean((true_tau < ci_lower_ms_rlearner_mMerge_aMerge | true_tau > ci_upper_ms_rlearner_mMerge_aMerge))
  ci_length_ms_rlearner_mMerge_aMerge <- mean(ci_upper_ms_rlearner_mMerge_aMerge - ci_lower_ms_rlearner_mMerge_aMerge)
  tau_hat_mean_ms_rlearner_mMerge_aMerge <- colMeans(tau_boot_mat_ms_rlearner_mMerge_aMerge)  # length n_test
  bias_squared_pointwise_ms_rlearner_mMerge_aMerge <- (tau_hat_mean_ms_rlearner_mMerge_aMerge - true_tau)^2
  variance_pointwise_ms_rlearner_mMerge_aMerge <- apply(tau_boot_mat_ms_rlearner_mMerge_aMerge, 2, var)
  mse_pointwise_ms_rlearner_mMerge_aMerge <- bias_squared_pointwise_ms_rlearner_mMerge_aMerge + variance_pointwise_ms_rlearner_mMerge_aMerge
  bias_squared_ms_rlearner_mMerge_aMerge <- mean(bias_squared_pointwise_ms_rlearner_mMerge_aMerge)
  variance_hat_tau_ms_rlearner_mMerge_aMerge <- mean(variance_pointwise_ms_rlearner_mMerge_aMerge)
  mse_ms_rlearner_mMerge_aMerge <- mean(mse_pointwise_ms_rlearner_mMerge_aMerge)
  
  ############################################
  ######### MS R-learner mSS aSS  #########
  ############################################
  se_tau_ms_rlearner_mSS_aSS <- apply(tau_boot_mat_ms_rlearner_mSS_aSS, 2, sd)
  ci_lower_ms_rlearner_mSS_aSS <- apply(tau_boot_mat_ms_rlearner_mSS_aSS, 2, quantile, probs = 0.025)
  ci_upper_ms_rlearner_mSS_aSS <- apply(tau_boot_mat_ms_rlearner_mSS_aSS, 2, quantile, probs = 0.975)
  coverage_ms_rlearner_mSS_aSS <- mean((true_tau >= ci_lower_ms_rlearner_mSS_aSS) & (true_tau <= ci_upper_ms_rlearner_mSS_aSS))
  type1_error_ms_rlearner_mSS_aSS <- mean((true_tau < ci_lower_ms_rlearner_mSS_aSS | true_tau > ci_upper_ms_rlearner_mSS_aSS))
  ci_length_ms_rlearner_mSS_aSS <- mean(ci_upper_ms_rlearner_mSS_aSS - ci_lower_ms_rlearner_mSS_aSS)
  tau_hat_mean_ms_rlearner_mSS_aSS <- colMeans(tau_boot_mat_ms_rlearner_mSS_aSS)  # length n_test
  bias_squared_pointwise_ms_rlearner_mSS_aSS <- (tau_hat_mean_ms_rlearner_mSS_aSS - true_tau)^2
  variance_pointwise_ms_rlearner_mSS_aSS <- apply(tau_boot_mat_ms_rlearner_mSS_aSS, 2, var)
  mse_pointwise_ms_rlearner_mSS_aSS <- bias_squared_pointwise_ms_rlearner_mSS_aSS + variance_pointwise_ms_rlearner_mSS_aSS
  bias_squared_ms_rlearner_mSS_aSS <- mean(bias_squared_pointwise_ms_rlearner_mSS_aSS)
  variance_hat_tau_ms_rlearner_mSS_aSS <- mean(variance_pointwise_ms_rlearner_mSS_aSS)
  mse_ms_rlearner_mSS_aSS <- mean(mse_pointwise_ms_rlearner_mSS_aSS)
  
  ############################################
  ######### MS R-learner mMerge aSS     ######
  ############################################
  se_tau_ms_rlearner_mMerge_aSS <- apply(tau_boot_mat_ms_rlearner_mMerge_aSS, 2, sd)
  ci_lower_ms_rlearner_mMerge_aSS <- apply(tau_boot_mat_ms_rlearner_mMerge_aSS, 2, quantile, probs = 0.025)
  ci_upper_ms_rlearner_mMerge_aSS <- apply(tau_boot_mat_ms_rlearner_mMerge_aSS, 2, quantile, probs = 0.975)
  coverage_ms_rlearner_mMerge_aSS <- mean((true_tau >= ci_lower_ms_rlearner_mMerge_aSS) & (true_tau <= ci_upper_ms_rlearner_mMerge_aSS))
  type1_error_ms_rlearner_mMerge_aSS <- mean((true_tau < ci_lower_ms_rlearner_mMerge_aSS | true_tau > ci_upper_ms_rlearner_mMerge_aSS))
  ci_length_ms_rlearner_mMerge_aSS <- mean(ci_upper_ms_rlearner_mMerge_aSS - ci_lower_ms_rlearner_mMerge_aSS)
  tau_hat_mean_ms_rlearner_mMerge_aSS <- colMeans(tau_boot_mat_ms_rlearner_mMerge_aSS)  # length n_test
  bias_squared_pointwise_ms_rlearner_mMerge_aSS <- (tau_hat_mean_ms_rlearner_mMerge_aSS - true_tau)^2
  variance_pointwise_ms_rlearner_mMerge_aSS <- apply(tau_boot_mat_ms_rlearner_mMerge_aSS, 2, var)
  mse_pointwise_ms_rlearner_mMerge_aSS <- bias_squared_pointwise_ms_rlearner_mMerge_aSS + variance_pointwise_ms_rlearner_mMerge_aSS
  bias_squared_ms_rlearner_mMerge_aSS <- mean(bias_squared_pointwise_ms_rlearner_mMerge_aSS)
  variance_hat_tau_ms_rlearner_mMerge_aSS <- mean(variance_pointwise_ms_rlearner_mMerge_aSS)
  mse_ms_rlearner_mMerge_aSS <- mean(mse_pointwise_ms_rlearner_mMerge_aSS)
  
  ############################################
  ######### MS R-learner mAvg aSS     ######
  ############################################
  se_tau_ms_rlearner_mAvg_aSS <- apply(tau_boot_mat_ms_rlearner_mAvg_aSS, 2, sd)
  ci_lower_ms_rlearner_mAvg_aSS <- apply(tau_boot_mat_ms_rlearner_mAvg_aSS, 2, quantile, probs = 0.025)
  ci_upper_ms_rlearner_mAvg_aSS <- apply(tau_boot_mat_ms_rlearner_mAvg_aSS, 2, quantile, probs = 0.975)
  coverage_ms_rlearner_mAvg_aSS <- mean((true_tau >= ci_lower_ms_rlearner_mAvg_aSS) & (true_tau <= ci_upper_ms_rlearner_mAvg_aSS))
  type1_error_ms_rlearner_mAvg_aSS <- mean((true_tau < ci_lower_ms_rlearner_mAvg_aSS | true_tau > ci_upper_ms_rlearner_mAvg_aSS))
  ci_length_ms_rlearner_mAvg_aSS <- mean(ci_upper_ms_rlearner_mAvg_aSS - ci_lower_ms_rlearner_mAvg_aSS)
  tau_hat_mean_ms_rlearner_mAvg_aSS <- colMeans(tau_boot_mat_ms_rlearner_mAvg_aSS)  # length n_test
  bias_squared_pointwise_ms_rlearner_mAvg_aSS <- (tau_hat_mean_ms_rlearner_mAvg_aSS - true_tau)^2
  variance_pointwise_ms_rlearner_mAvg_aSS <- apply(tau_boot_mat_ms_rlearner_mAvg_aSS, 2, var)
  mse_pointwise_ms_rlearner_mAvg_aSS <- bias_squared_pointwise_ms_rlearner_mAvg_aSS + variance_pointwise_ms_rlearner_mAvg_aSS
  bias_squared_ms_rlearner_mAvg_aSS <- mean(bias_squared_pointwise_ms_rlearner_mAvg_aSS)
  variance_hat_tau_ms_rlearner_mAvg_aSS <- mean(variance_pointwise_ms_rlearner_mAvg_aSS)
  mse_ms_rlearner_mAvg_aSS <- mean(mse_pointwise_ms_rlearner_mAvg_aSS)
  
  ############################################
  ######### MS R-learner mStack aSS     ######
  ############################################
  se_tau_ms_rlearner_mStack_aSS <- apply(tau_boot_mat_ms_rlearner_mStack_aSS, 2, sd)
  ci_lower_ms_rlearner_mStack_aSS <- apply(tau_boot_mat_ms_rlearner_mStack_aSS, 2, quantile, probs = 0.025)
  ci_upper_ms_rlearner_mStack_aSS <- apply(tau_boot_mat_ms_rlearner_mStack_aSS, 2, quantile, probs = 0.975)
  coverage_ms_rlearner_mStack_aSS <- mean((true_tau >= ci_lower_ms_rlearner_mStack_aSS) & (true_tau <= ci_upper_ms_rlearner_mStack_aSS))
  type1_error_ms_rlearner_mStack_aSS <- mean((true_tau < ci_lower_ms_rlearner_mStack_aSS | true_tau > ci_upper_ms_rlearner_mStack_aSS))
  ci_length_ms_rlearner_mStack_aSS <- mean(ci_upper_ms_rlearner_mStack_aSS - ci_lower_ms_rlearner_mStack_aSS)
  tau_hat_mean_ms_rlearner_mStack_aSS <- colMeans(tau_boot_mat_ms_rlearner_mStack_aSS)  # length n_test
  bias_squared_pointwise_ms_rlearner_mStack_aSS <- (tau_hat_mean_ms_rlearner_mStack_aSS - true_tau)^2
  variance_pointwise_ms_rlearner_mStack_aSS <- apply(tau_boot_mat_ms_rlearner_mStack_aSS, 2, var)
  mse_pointwise_ms_rlearner_mStack_aSS <- bias_squared_pointwise_ms_rlearner_mStack_aSS + variance_pointwise_ms_rlearner_mStack_aSS
  bias_squared_ms_rlearner_mStack_aSS <- mean(bias_squared_pointwise_ms_rlearner_mStack_aSS)
  variance_hat_tau_ms_rlearner_mStack_aSS <- mean(variance_pointwise_ms_rlearner_mStack_aSS)
  mse_ms_rlearner_mStack_aSS <- mean(mse_pointwise_ms_rlearner_mStack_aSS)
  
  out <- data.frame(
    # R-learner oracle SS
    rlearner_oracle_SS_coverage = coverage_rlearner_oracle_SS,
    rlearner_oracle_SS_type_I_error = type1_error_rlearner_oracle_SS,
    rlearner_oracle_SS_CI_length = ci_length_rlearner_oracle_SS,
    rlearner_oracle_SS_MSE = mse_rlearner_oracle_SS,
    rlearner_oracle_SS_bias_sq = bias_squared_rlearner_oracle_SS,
    rlearner_oracle_SS_variance_hat = variance_hat_tau_rlearner_oracle_SS,
    
    # R-learner SS
    rlearner_estimate_SS_coverage = coverage_rlearner_estimate_SS,
    rlearner_estimate_SS_type_I_error = type1_error_rlearner_estimate_SS,
    rlearner_estimate_SS_CI_length = ci_length_rlearner_estimate_SS,
    rlearner_estimate_SS_MSE = mse_rlearner_estimate_SS,
    rlearner_estimate_SS_bias_sq = bias_squared_rlearner_estimate_SS,
    rlearner_estimate_SS_variance_hat = variance_hat_tau_rlearner_estimate_SS,
    
    # MS R-learner Oracle
    ms_rlearner_oracle_coverage = coverage_ms_rlearner_oracle,
    ms_rlearner_oracle_type_I_error = type1_error_ms_rlearner_oracle,
    ms_rlearner_oracle_CI_length = ci_length_ms_rlearner_oracle,
    ms_rlearner_oracle_MSE = mse_ms_rlearner_oracle,
    ms_rlearner_oracle_bias_sq = bias_squared_ms_rlearner_oracle,
    ms_rlearner_oracle_variance_hat = variance_hat_tau_ms_rlearner_oracle,
    
    # MS R-learner mMerge aMerge
    ms_rlearner_mMerge_aMerge_coverage = coverage_ms_rlearner_mMerge_aMerge,
    ms_rlearner_mMerge_aMerge_type_I_error = type1_error_ms_rlearner_mMerge_aMerge,
    ms_rlearner_mMerge_aMerge_CI_length = ci_length_ms_rlearner_mMerge_aMerge,
    ms_rlearner_mMerge_aMerge_MSE = mse_ms_rlearner_mMerge_aMerge,
    ms_rlearner_mMerge_aMerge_bias_sq = bias_squared_ms_rlearner_mMerge_aMerge,
    ms_rlearner_mMerge_aMerge_variance_hat = variance_hat_tau_ms_rlearner_mMerge_aMerge,
    
    # MS R-learner mSS aSS
    ms_rlearner_mSS_aSS_coverage = coverage_ms_rlearner_mSS_aSS,
    ms_rlearner_mSS_aSS_type_I_error = type1_error_ms_rlearner_mSS_aSS,
    ms_rlearner_mSS_aSS_CI_length = ci_length_ms_rlearner_mSS_aSS,
    ms_rlearner_mSS_aSS_MSE = mse_ms_rlearner_mSS_aSS,
    ms_rlearner_mSS_aSS_bias_sq = bias_squared_ms_rlearner_mSS_aSS,
    ms_rlearner_mSS_aSS_variance_hat = variance_hat_tau_ms_rlearner_mSS_aSS,
    
    # MS R-learner mMerge aSS
    ms_rlearner_mMerge_aSS_coverage = coverage_ms_rlearner_mMerge_aSS,
    ms_rlearner_mMerge_aSS_type_I_error = type1_error_ms_rlearner_mMerge_aSS,
    ms_rlearner_mMerge_aSS_CI_length = ci_length_ms_rlearner_mMerge_aSS,
    ms_rlearner_mMerge_aSS_MSE = mse_ms_rlearner_mMerge_aSS,
    ms_rlearner_mMerge_aSS_bias_sq = bias_squared_ms_rlearner_mMerge_aSS,
    ms_rlearner_mMerge_aSS_variance_hat = variance_hat_tau_ms_rlearner_mMerge_aSS,
    
    # MS R-learner mAvg aSS
    ms_rlearner_mAvg_aSS_coverage = coverage_ms_rlearner_mAvg_aSS,
    ms_rlearner_mAvg_aSS_type_I_error = type1_error_ms_rlearner_mAvg_aSS,
    ms_rlearner_mAvg_aSS_CI_length = ci_length_ms_rlearner_mAvg_aSS,
    ms_rlearner_mAvg_aSS_MSE = mse_ms_rlearner_mAvg_aSS,
    ms_rlearner_mAvg_aSS_bias_sq = bias_squared_ms_rlearner_mAvg_aSS,
    ms_rlearner_mAvg_aSS_variance_hat = variance_hat_tau_ms_rlearner_mAvg_aSS,
    
    # MS R-learner mStack aSS
    ms_rlearner_mStack_aSS_coverage = coverage_ms_rlearner_mStack_aSS,
    ms_rlearner_mStack_aSS_type_I_error = type1_error_ms_rlearner_mStack_aSS,
    ms_rlearner_mStack_aSS_CI_length = ci_length_ms_rlearner_mStack_aSS,
    ms_rlearner_mStack_aSS_MSE = mse_ms_rlearner_mStack_aSS,
    ms_rlearner_mStack_aSS_bias_sq = bias_squared_ms_rlearner_mStack_aSS,
    ms_rlearner_mStack_aSS_variance_hat = variance_hat_tau_ms_rlearner_mStack_aSS,
    
    perturb_tau = perturb_tau,
    perturb_a = perturb_a,
    perturb_mu0 = perturb_mu0,
    p = p)
  
  metrics_to_display <- c("coverage", "MSE", "bias_sq", "variance_hat")
  methods <- c("rlearner_oracle_SS",
               "rlearner_estimate_SS", 
               "ms_rlearner_oracle",
               "ms_rlearner_mMerge_aSS", 
               "ms_rlearner_mSS_aSS", 
               "ms_rlearner_mAvg_aSS", 
               "ms_rlearner_mStack_aSS")
  
  summary_table <- do.call(rbind, lapply(methods, function(meth) {
    sapply(metrics_to_display, function(metric) {
      out[[paste0(meth, "_", metric)]]
    })
  }))
  summary_df <- as.data.frame(summary_table)
  summary_df$Method <- methods
  summary_df <- summary_df[, c("Method", metrics_to_display)]
  print(summary_df, row.names = FALSE, digits = 4)
  return(out)
}
# Scenario B
# true tau_k is nonlinear (cubic spline with knots)

# Load packages and data
source("Utils.R")
load("edat_orig.RData")

# Remove small studies with < 150 observations
edat_orig <- edat_orig[which(sapply(edat_orig, function(x) nrow(x)) >= 150)]

# Number of studies
K <- 4

# Total number of observations
n_k <- 150
n <- n_k * (K)

# Total number of simulation runs
nreps <- 500

# Number of covariates
p <- 40

covariate_cols <- sapply(1:p, function(x) paste("X.", x, sep = ""))

# Studies
set.seed(1)
edat_source <- init_data(edat_orig[1:5], ndat = K, nvar = p)
for (i in 1:length(edat_source)) {
  dataset = edat_source[[i]]
  dataset = dataset[1:n_k,]
  edat_source[[i]] = cbind(dataset[, 1:p])
}
edat_all <- do.call("rbind", edat_source)

######################################################################################
# Simulation function for calculating the MSE of conditional average treatment effect
# 
# Input: 
# n - total number of observations
# p - number of predictors
# K - number of source studies
# data - all studies combined
# covariate_cols - names of predictors
# X_int_basis - basis-expanded predictor matrix with intercept
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
#
# Output:
# Data frame of the performance results (MSE, bias, variance) for the multi-study and study-specific R-learners
# across the oracle and plug-in settings
######################################################################################
sim_each <- function(n, p, K, data, covariate_cols, X_int_basis, cols_re_tau, cols_re_mu0, cols_re_a, beta_tau, beta_mu0, beta_a, beta_mu0_k, beta_a_k, beta_tau_k, p_mat, perturb_tau, perturb_mu0, perturb_a, k_folds, alpha){
  check <- rep(0, K)
  while(any(check == 0)){
    # Test set indices
    ind_test <- sample(1:n, 0.3 * n, replace = FALSE)
    
    # Initialize vectors to store coefficients
    tau_k <- mu0_k <- p_k <- vector("list", length = K)
    X_int <- as.matrix(cbind(1, data[, covariate_cols]))
    for(k in 1:(K)){
      # Tau
      sigma_tau <- rep(0, ncol(X_int_basis))
      sigma_tau[cols_re_tau[[k]]] <- perturb_tau
      beta_tau_k[[k]] <- beta_tau_k[[k]] + rnorm(ncol(X_int_basis), 0, sigma_tau)
      tau_k[[k]] <- as.matrix(X_int_basis[data$S == k, ]) %*% beta_tau_k[[k]]
      
      # mu0 := E(Y|X)
      sigma_mu0 <- rep(0, p + 1)
      sigma_mu0[cols_re_mu0[[k]]] <- perturb_mu0
      beta_mu0_k[[k]] <- beta_mu0_k[[k]] + rnorm(p + 1, 0, sigma_mu0)
      mu0_k[[k]] <- X_int[data$S == k, ] %*% beta_mu0_k[[k]] 
      
      # e := P(A|X)
      sigma_a <- rep(0, p + 1)
      sigma_a[cols_re_a[[k]]] <- perturb_a
      beta_a_k[[k]] <- beta_a_k[[k]] + rnorm(p + 1, 0, sigma_a)
    }
    
    # A|X, S=k is 0/1 Bernoulli distributed with probability e_k(X)
    true_e <- vector("list", length = K)
    for(k in 1:K){
      X_int_k <- as.matrix(cbind(1, data[which(data$S == k), covariate_cols]))
      true_e[[k]] <- expit(X_int_k %*% beta_a_k[[k]])
      data$A[which(data$S == k)] <- rbinom(nrow(data[which(data$S == k), ]), 1, true_e[[k]]) 
      check[k] <- ifelse(length(which(data$A[which(data$S == k)] == 1)) < 10 | length(which(data$A[which(data$S == k)] == 0)) < 10, 0, 1)
    }
  }
  
  
  # Generate Y
  Y_mat <- sapply(1:K, function(k){
    (X_int %*% beta_mu0_k[[k]] + data$A * as.matrix(X_int_basis) %*% beta_tau_k[[k]]) * p_mat[, k]
  })
  data$Y <- apply(Y_mat, 1, sum) + rnorm(n, 0, 1)
  
  # Oracle m
  m_mat <- sapply(1:K, function(k){
    (X_int %*% beta_mu0_k[[k]] + expit(X_int %*% beta_a_k[[k]]) * as.matrix(X_int_basis) %*% beta_tau_k[[k]]) * p_mat[, k]
  })
  true_m <- apply(m_mat, 1, sum)
  
  # Estimate m (merge)
  m_estimate_merge <- rlearner_m_hat_merge(data = data[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # Estimate a 
  a_estimate_merge <- rlearner_a_hat_merge(data = data[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # Oracle tau
  tau_mat <- sapply(1:K, function(k){
    (as.matrix(X_int_basis) %*% beta_tau_k[[k]]) * p_mat[, k]
  })
  true_tau <- apply(tau_mat, 1, sum)[ind_test]
  
  
  ##############################################################
  #                      R-learner Oracle                      #
  ##############################################################
  X_basis <- X_int_basis[, 2:ncol(X_int_basis)]
  
  rlearner_oracle <- rlasso(x = as.matrix(X_basis[-ind_test, ]), 
                            w = data[-ind_test, ]$A, 
                            y = data[-ind_test, ]$Y,
                            lambda_choice = "lambda.min",
                            p_hat = unlist(true_e)[-ind_test],
                            m_hat = true_m[-ind_test],
                            alpha = alpha)
  
  rlearner_oracle_pred <- predict(rlearner_oracle, newx = as.matrix(X_basis[ind_test, ]))
  
  ################################################################
  #                      R-learner Estimate                      #
  ################################################################
  rlearner_merge <- rlasso(x = as.matrix(X_basis[-ind_test, ]), 
                           w = data[-ind_test, ]$A, 
                           y = data[-ind_test, ]$Y,
                           lambda_choice = "lambda.min",
                           p_hat = a_estimate_merge,
                           m_hat = m_estimate_merge,
                           alpha = alpha)
  
  rlearner_merge_pred <- predict(rlearner_merge, newx = as.matrix(X_basis[ind_test, ]))
  
  #############################################################################
  #                      R-learner Oracle Study-Specific                      #
  #############################################################################
  rlearner_oracle_SS <- rlearner_oracle_SS_pred <- vector("list", length = K)
  for(k in 1:K){
    data_k_ind <- setdiff(which(data$S == k), ind_test)
    data_k <- data[-ind_test, ]
    data_k <- data_k[which(data_k$S == k), ]
    rlearner_oracle_SS[[k]] <- rlasso(x = as.matrix(data_k[, covariate_cols]), 
                                      w = data_k$A, 
                                      y = data_k$Y,
                                      lambda_choice = "lambda.min",
                                      p_hat = true_e[[k]][which(as.numeric(rownames(true_e[[k]])) %in% data_k_ind)],
                                      m_hat = true_m[data_k_ind],
                                      alpha = alpha)
    rlearner_oracle_SS_pred[[k]] <- predict(rlearner_oracle_SS[[k]], newx = as.matrix(data[ind_test, covariate_cols]))
  }
  
  rlearner_oracle_SS_pred <- apply(do.call("cbind", rlearner_oracle_SS_pred) * p_mat[ind_test, ], 1, sum)
  
  #################################################################################
  #                      Multi-Study R-Learner Oracle (both)                      # 
  #################################################################################
  
  X_res_mat_with_int <- lapply(1:K, function(k){
    sweep(X_int_basis[-ind_test, ], 1, (data[-ind_test, ]$A - expit(X_int[-ind_test, ] %*% beta_a_k[[k]])) * p_mat[-ind_test, k], "*")
  })
  
  
  X_res_all_oracle <- as.data.frame(do.call("cbind", X_res_mat_with_int))
  names(X_res_all_oracle) <- sapply(1:ncol(X_res_all_oracle), function(x) paste("X.", x, sep = ""))
  
  Y_res_merge_oracle <- data[-ind_test, ]$Y - true_m[-ind_test]
  ms_rlearner_oracle_pred <- ms_rlearner(X = X_res_all_oracle, Y = Y_res_merge_oracle, alpha = alpha, K = K, p = ncol(X_int_basis) - 1, p_mat = p_mat[ind_test, ],
                                         newX = as.matrix(X_int_basis[ind_test, ]))
  
  ####################################################################################
  #                      Multi-Study R-Learner Merge (Both)                          # 
  ####################################################################################
  
  # Estimate ascertainment probability on training data
  p_mod <- glmnet::cv.glmnet(x = as.matrix(data[-ind_test, covariate_cols]),
                             y = as.matrix(data$S[-ind_test]),
                             alpha = alpha,
                             lambda = NULL,
                             standardize = TRUE,
                             keep = TRUE,
                             family = "multinomial")
  p_hat_lambda_min = p_mod$lambda[which.min(p_mod$cvm[!is.na(colSums(p_mod$fit.preval))])]
  p_pred_mat_train <- as.data.frame(predict(p_mod, newx = as.matrix(data[-ind_test, covariate_cols]), s = p_hat_lambda_min, type = "response"))
  p_pred_mat_test <- as.data.frame(predict(p_mod, newx = as.matrix(data[ind_test, covariate_cols]), s = p_hat_lambda_min, type = "response"))
  
  X_res_mat_merge <- lapply(1:K, function(k){
    sweep(X_int_basis[-ind_test, ], 1, (data[-ind_test, ]$A - a_estimate_merge) * p_pred_mat_train[, k], "*")
  })
  
  X_res_mat_merge_all <- as.data.frame(do.call("cbind", X_res_mat_merge))
  names(X_res_mat_merge_all) <- sapply(1:ncol(X_res_mat_merge_all), function(x) paste("X.", x, sep = ""))
  
  Y_res_merge_estimate <- data[-ind_test, ]$Y - m_estimate_merge
  ms_rlearner_mMerge_aMerge_pred <- ms_rlearner(X = X_res_mat_merge_all, 
                                                Y = Y_res_merge_estimate, 
                                                alpha = alpha, K = K, p = ncol(X_int_basis) - 1, 
                                                p_mat = p_pred_mat_test,
                                                newX = as.matrix(X_int_basis[ind_test, ]))
  
  ####################################################################################
  #               Multi-Study R-Learner Nuisance Study-Specific                      # 
  ####################################################################################
  m_estimate_study_specific <- rlearner_m_hat_stack(data = data[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha, p_mat = p_pred_mat_train)
  m_estimate_avg <- m_estimate_study_specific$avg
  m_estimate_stack <- m_estimate_study_specific$stack_int
  m_estimate_ss <- m_estimate_study_specific$ss
  
  # Propensity score model e_k(x)
  a_hat_k <- vector("list", length = K)
  a_hat_k_lambda_min <- vector()
  data_train <- data[-ind_test, ]
  for(k in 1:K){
    a_hat_k[[k]] <- glmnet::cv.glmnet(x = as.matrix(data_train[which(data_train$S == k), covariate_cols]),
                                      y = as.matrix(data_train$A[which(data_train$S == k)]),
                                      alpha = alpha,
                                      lambda = NULL,
                                      standardize = TRUE,
                                      keep = TRUE,
                                      family = "binomial",
                                      type.measure = "deviance")
    a_hat_k_lambda_min[k] = a_hat_k[[k]]$lambda[which.min(a_hat_k[[k]]$cvm[!is.na(colSums(a_hat_k[[k]]$fit.preval))])]
  }
  
  X_res_mat_avg <- lapply(1:K, function(k){
    sweep(X_int_basis[-ind_test, ], 1, (data[-ind_test, ]$A - predict(a_hat_k[[k]], newx = as.matrix(data[-ind_test, covariate_cols]), s = a_hat_k_lambda_min[k], type = "response")) * p_pred_mat_train[, k], "*")
  })
  
  X_res_mat_avg_all <- as.data.frame(do.call("cbind", X_res_mat_avg))
  names(X_res_mat_avg_all) <- sapply(1:ncol(X_res_mat_avg_all), function(x) paste("X.", x, sep = ""))
  
  # Merge m, Study-specific A
  ms_rlearner_mMerge_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_merge_estimate, alpha = alpha, K = K, p = ncol(X_int_basis) - 1, p_mat = p_pred_mat_test,
                                             newX = as.matrix(X_int_basis[ind_test, ]))
  
  # Avg m, Study-specific A
  Y_res_avg_estimate <- data[-ind_test, ]$Y - m_estimate_avg
  ms_rlearner_mAvg_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_avg_estimate, alpha = alpha, K = K, p = ncol(X_int_basis) - 1, p_mat = p_pred_mat_test,
                                           newX = as.matrix(X_int_basis[ind_test, ]))
  
  # Stack m, Study-specific A
  Y_res_stack_estimate <- data[-ind_test, ]$Y - m_estimate_stack
  ms_rlearner_mStack_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_stack_estimate, alpha = alpha, K = K, p = ncol(X_int_basis) - 1, p_mat = p_pred_mat_test,
                                             newX = as.matrix(X_int_basis[ind_test, ]))
  
  # Study-specific m, Study-specific A
  Y_res_ss_estimate <- data[-ind_test, ]$Y - m_estimate_ss
  ms_rlearner_mSS_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_ss_estimate, alpha = alpha, K = K, p = ncol(X_int_basis) - 1, p_mat = p_pred_mat_test,
                                          newX = as.matrix(X_int_basis[ind_test, ]))
  
  
  #############################################################################
  #                      R-learner Estimate Study-Specific                    #
  #############################################################################
  rlearner_estimate_SS <- rlearner_estimate_SS_pred <- vector("list", length = K)
  for(k in 1:K){
    data_k_ind <- setdiff(which(data$S == k), ind_test)
    data_k <- data[-ind_test, ]
    data_k <- data_k[which(data_k$S == k), ]
    rlearner_estimate_SS[[k]] <- rlasso(x = as.matrix(data_k[, covariate_cols]), 
                                        w = data_k$A, 
                                        y = data_k$Y,
                                        lambda_choice = "lambda.min",
                                        p_hat = predict(a_hat_k[[k]], newx = as.matrix(data_k[, covariate_cols]), s = a_hat_k_lambda_min[k], type = "response"),
                                        m_hat = m_estimate_study_specific$ss_by_study[[k]],
                                        alpha = alpha)
    rlearner_estimate_SS_pred[[k]] <- predict(rlearner_estimate_SS[[k]], newx = as.matrix(data[ind_test, covariate_cols]))
  }
  
  rlearner_estimate_SS_pred <- apply(do.call("cbind", rlearner_estimate_SS_pred) * p_pred_mat_test, 1, sum)
  
  
  out <- data.frame(# R-learner oracle merge
    rlearner_oracle_merge_MSE = mean((true_tau - rlearner_oracle_pred)^2), 
    rlearner_oracle_merge_bias2 = mean(rlearner_oracle_pred - true_tau)^2,
    rlearner_oracle_merge_var = var(rlearner_oracle_pred),
    # R-learner oracle SS
    rlearner_oracle_SS_MSE = mean((true_tau - rlearner_oracle_SS_pred)^2), 
    rlearner_oracle_SS_bias2 = mean(rlearner_oracle_SS_pred - true_tau)^2,
    rlearner_oracle_SS_var = var(rlearner_oracle_SS_pred),
    # R-learner estimate merge
    rlearner_estimate_merge_MSE = mean((true_tau - rlearner_merge_pred)^2), 
    rlearner_estimate_merge_bias2 = mean(rlearner_merge_pred - true_tau)^2,
    rlearner_estimate_merge_var = var(rlearner_merge_pred),
    # R-learner estimate SS
    rlearner_estimate_SS_MSE = mean((true_tau - rlearner_estimate_SS_pred)^2), 
    rlearner_estimate_SS_bias2 = mean(rlearner_estimate_SS_pred - true_tau)^2,
    rlearner_estimate_SS_var = var(rlearner_estimate_SS_pred),
    # MS-rlearner oracle (merge both)
    ms_rlearner_oracle_MSE = mean((true_tau - ms_rlearner_oracle_pred)^2),
    ms_rlearner_oracle_bias2 = mean(ms_rlearner_oracle_pred - true_tau)^2,
    ms_rlearner_oracle_var = var(ms_rlearner_oracle_pred),
    # MS-rlearner mMerge aMerge
    ms_rlearner_mMerge_aMerge_MSE = mean((true_tau - ms_rlearner_mMerge_aMerge_pred)^2),
    ms_rlearner_mMerge_aMerge_bias2 = mean(ms_rlearner_mMerge_aMerge_pred - true_tau)^2,
    ms_rlearner_mMerge_aMerge_var = var(ms_rlearner_mMerge_aMerge_pred),
    # MS-rlearner mMerge aSS
    ms_rlearner_mMerge_aSS_MSE = mean((true_tau - ms_rlearner_mMerge_aSS_pred)^2),
    ms_rlearner_mMerge_aSS_bias2 = mean(ms_rlearner_mMerge_aSS_pred - true_tau)^2,
    ms_rlearner_mMerge_aSS_var = var(ms_rlearner_mMerge_aSS_pred),
    # MS-rlearner mAvg aSS
    ms_rlearner_mAvg_aSS_MSE = mean((true_tau - ms_rlearner_mAvg_aSS_pred)^2),
    ms_rlearner_mAvg_aSS_bias2 = mean(ms_rlearner_mAvg_aSS_pred - true_tau)^2,
    ms_rlearner_mAvg_aSS_var = var(ms_rlearner_mAvg_aSS_pred),
    # MS-rlearner mStack aSS
    ms_rlearner_mStack_aSS_MSE = mean((true_tau - ms_rlearner_mStack_aSS_pred)^2),
    ms_rlearner_mStack_aSS_bias2 = mean(ms_rlearner_mStack_aSS_pred - true_tau)^2,
    ms_rlearner_mStack_aSS_var = var(ms_rlearner_mStack_aSS_pred),
    # MS-rlearner mSS aSS
    ms_rlearner_mSS_aSS_MSE = mean((true_tau - ms_rlearner_mSS_aSS_pred)^2),
    ms_rlearner_mSS_aSS_bias2 = mean(ms_rlearner_mSS_aSS_pred - true_tau)^2,
    ms_rlearner_mSS_aSS_var = var(ms_rlearner_mSS_aSS_pred),
    perturb_tau = perturb_tau,
    perturb_a = perturb_a,
    perturb_mu0 = perturb_mu0,
    p = p)
  return(out)
}

sim_multiple <- function(edat_all, n, p, K, covariate_cols, perturb_tau, perturb_mu0, perturb_a, k_folds, alpha, numRE, nreps){
  # Coefficients for generating mu0(x) = E[Y(0)|X=x]
  beta_mu0 <- rnorm(p + 1, mean = 0, sd = 1)
  
  # Coefficients for generating e(x) = E[A|X=x]
  beta_a <- rnorm(p + 1, mean = 0, sd = 1)
  
  # Basis-expanded predictor matrix
  spline.design = bs(edat_all[, 1])
  colnames(spline.design) = paste0("bs1_", 1:ncol(spline.design))
  edat_all_basis = cbind(edat_all[, 2:ncol(edat_all)], spline.design)
  X_int_basis <- cbind(1, edat_all_basis)
  
  # Coefficients for generating tau(x) = E[Y(1)-Y(0)|X=x]
  beta_tau <- rnorm(p + ncol(spline.design), mean = 0, sd = 1)

  # Coefficients for generating p(k|x) = P(S=k|X=x)
  beta_p_k <- vector("list", length = K - 1)
  for(k in 1:(K - 1)){
    beta_p_k[[k]] <- rnorm(p + 1, mean = 0, sd = 1)   
    beta_p_k[[k]][sample(2:(p + 1), round(p * (4/5)))] <- 0
  }
  
  # p := P(S|X = x) ascertainment probability
  X_int <- cbind(1, edat_all)
  
  denom_mat <- sapply(1:length(beta_p_k), function(x){
    exp(as.matrix(X_int) %*% beta_p_k[[x]])
  })
  denom <- apply(denom_mat, 1, sum)
  
  p_vec_dat <- data.frame(sapply(1:length(beta_p_k), function(x){
    exp(as.matrix(X_int) %*% beta_p_k[[x]])/(1 + denom)
  }))
  names(p_vec_dat) <- sapply(2:K, function(x) paste("p", x, sep = ""))
  p1 <- 1 - apply(p_vec_dat, 1, sum)
  data <- data.frame(cbind(edat_all, p1 = p1, p_vec_dat))
  data <- data %>% 
    rowwise %>%
    dplyr::mutate(S = sample(1:K, size = 1, replace = TRUE, 
                             prob = c_across(matches('^p\\d+$')))) %>% 
    ungroup
  data <- as.data.frame(data[, c(covariate_cols, "S")])
  # data <- data[order(data$S), ]
  p_mat <- cbind(p1, p_vec_dat)
  
  beta_mu0_k <- beta_a_k <- beta_tau_k <- coef_ind <- vector("list", length = K)
  # For selecting study-specific numRE
  coef_ind[[1]] <- c()
  
  # Beta_mu0_k
  for(k in 1:(K)){
    beta_mu0_k[[k]] <- beta_mu0
  }
  
  # Beta_a_k
  for(k in 1:(K)){
    beta_a_k[[k]] <- beta_a
  }
  
  # Beta_tau_k
  for(k in 1:(K)){
    beta_tau_k[[k]] <- beta_tau
  }
  
  # Indices of random effects
  cols_re_tau <- cols_re_mu0 <-cols_re_a <- vector("list", length = K)
  for(k in 1:(K)){
    cols_re_tau[[k]] <- sample(1:(p + 1), numRE)
    
    cols_re_mu0[[k]] <- sample(1:(p + 1), numRE)
    
    cols_re_a[[k]] <- sample(1:(p + 1), numRE)
  }
  
  data$A <- NA
  
  registerDoMC(cores = 48)
  results = foreach (j = 1:nreps, .combine = rbind, .errorhandling = "remove") %dopar% {
    print(paste("Perturb_a =", perturb_a, "Perturb_mu0 =", perturb_mu0, "Perturb_tau =", perturb_tau, "Iteration =", j))
    sim_each(n = n, p = p, K = K, data = data, covariate_cols = covariate_cols, X_int_basis = X_int_basis,
             cols_re_tau = cols_re_tau, cols_re_mu0 = cols_re_mu0, cols_re_a = cols_re_a,
             beta_tau = beta_tau, beta_mu0 = beta_mu0, beta_a = beta_a, beta_mu0_k = beta_mu0_k, beta_a_k = beta_a_k, beta_tau_k = beta_tau_k, p_mat = p_mat,
             perturb_tau = perturb_tau, perturb_mu0 = perturb_mu0, perturb_a = perturb_a,
             k_folds = k_folds, alpha = alpha)
  }
}

perturb_a <- c(0, 1, 2)
perturb_tau <- seq(0, 2, 0.5)
perturb_mu0 <- c(0, 1, 2)
k_folds <- 3
alpha <- 0.5


# Number of random effects
numRE <- 8

results <- vector("list", length = length(perturb_a))
ind_a <- 1:length(perturb_a)
ind_mu0 <- 1:length(perturb_mu0)
ind_t <- 1:length(perturb_tau)
for(a in ind_a){
  results[[a]] <- vector("list", length = length(perturb_mu0)) 
  for(m in ind_mu0){
    results[[a]][[m]] <- vector("list", length = length(perturb_tau))  
    for(t in ind_t){
      results[[a]][[m]][[t]] <- as.data.frame(sim_multiple(edat_all = edat_all, n = n, p = p, K = K, covariate_cols = covariate_cols, 
                                                           perturb_tau = perturb_tau[t], perturb_mu0 = perturb_mu0[m], perturb_a = perturb_a[a], 
                                                           k_folds = k_folds, alpha = alpha,
                                                           numRE = numRE,
                                                           nreps = nreps))
    }
  }
}


save(results, perturb_tau, perturb_mu0, perturb_a, numRE, K, n, p, file = paste0("Scenario_B", ".RData"))


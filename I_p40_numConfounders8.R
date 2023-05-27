# Load packages and data
source("../Utils.R")
load("../edat_orig.RData")

# Remove small studies with < 100 observations
edat_orig <- edat_orig[which(sapply(edat_orig, function(x) nrow(x)) >= 100)]

# Number of studies
K <- 4

# Total number of observations
n_k <- 100
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
# cols_re_tau - column indices of the beta_tau vector with random effect
# cols_re_m - column indices of the beta_m vector with random effect
# cols_re_a - column indices of the beta_a vector with random effect
# beta_tau - coefficients for the heterogeneous treatment effect
# beta_m - coefficients for the mean outcome regression model
# beta_a - coefficients for the propensity score model
# beta_m_k - study-specific coefficients for the mean outcome regression model
# beta_a_k - study-specific coefficients for the propensity score model
# beta_p - list of coefficients for the study label model
# perturb_tau - variance of the random effects for beta_tau
# perturb_mu0 - variance of the random effects for beta_mu0
# perturb_a - variance of the random effects for beta_a
# k-folds - number of folds to perform cross-fitting 
# alpha - glmnet parameter (alpha = 1 => ridge, alpha = 0 => LASSO)
# ps - propensity score
# prop_overlap - proportion of overlap in confounders across studies
#
# Output:
# list of two:
# 1) mse 2) bias, variance
######################################################################################
sim_each <- function(n, p, K, data, covariate_cols, cols_re_tau, cols_re_mu0, cols_re_a, beta_tau, beta_mu0, beta_a, beta_mu0_k, beta_a_k, beta_tau_k, p_mat, perturb_tau, perturb_mu0, perturb_a, k_folds, alpha, prop_overlap){
  check <- rep(0, K)
  while(any(check == 0)){
    # Test set indices
    ind_test <- sample(1:n, 0.3 * n, replace = FALSE)
    
    # Initialize vectors to store coefficients
    tau_k <- mu0_k <- p_k <- vector("list", length = K)
    X_int <- as.matrix(cbind(1, data[, covariate_cols]))
    for(k in 1:(K)){
      # Tau
      sigma_tau <- rep(0, p + 1)
      sigma_tau[cols_re_tau[[k]]] <- perturb_tau
      beta_tau_k[[k]] <- beta_tau_k[[k]] + rnorm(p + 1, 0, sigma_tau)
      tau_k[[k]] <- X_int[data$S == k, ] %*% beta_tau_k[[k]]
      
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
      X_int_k <- as.matrix(cbind(1, data[data$S == k, covariate_cols]))
      true_e[[k]] <- rep(0.5, nrow(X_int_k))
      data$A[data$S == k] <- rbinom(nrow(X_int_k), 1, true_e[[k]]) 
      check[k] <- ifelse(length(which(data$A[which(data$S == k)] == 1)) < 10 | length(which(data$A[which(data$S == k)] == 0)) < 10, 0, 1)
    }
  }
  
  
  
  # Generate Y
  Y_mat <- sapply(1:K, function(k){
    (X_int %*% beta_mu0_k[[k]] + data$A * X_int %*% beta_tau_k[[k]]) * p_mat[, k]
  })
  data$Y <- apply(Y_mat, 1, sum) + rnorm(n, 0, 1)
  
  # Oracle m
  m_mat <- sapply(1:K, function(k){
    (X_int %*% beta_mu0_k[[k]] + 0.5 * X_int %*% beta_tau_k[[k]]) * p_mat[, k]
  })
  true_m <- apply(m_mat, 1, sum)
  
  # Estimate m (merge)
  m_estimate_merge <- rlearner_m_hat_merge(data = data[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # # Estimate a 
  # a_estimate_merge <- rlearner_a_hat_merge(data = data, K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # Oracle tau
  tau_mat <- sapply(1:K, function(k){
    (X_int %*% beta_tau_k[[k]]) * p_mat[, k]
  })
  true_tau <- apply(tau_mat, 1, sum)[ind_test]
  
  ##############################################################
  #                      R-learner Oracle                      #
  ##############################################################

  rlearner_oracle <- rlasso(x = as.matrix(data[-ind_test, covariate_cols]), 
                            w = data[-ind_test, ]$A, 
                            y = data[-ind_test, ]$Y,
                            lambda_choice = "lambda.min",
                            p_hat = unlist(true_e)[-ind_test],
                            m_hat = true_m[-ind_test],
                            alpha = alpha)
  
  rlearner_oracle_pred <- predict(rlearner_oracle, newx = as.matrix(data[ind_test, covariate_cols]))
  
  ################################################################
  #                      R-learner Estimate                      #
  ################################################################
  rlearner_merge <- rlasso(x = as.matrix(data[-ind_test, covariate_cols]), 
                           w = data[-ind_test, ]$A, 
                           y = data[-ind_test, ]$Y,
                           lambda_choice = "lambda.min",
                           p_hat = unlist(true_e)[-ind_test],
                           m_hat = m_estimate_merge,
                           alpha = alpha)
  
  rlearner_merge_pred <- predict(rlearner_merge, newx = as.matrix(data[ind_test, covariate_cols]))
  
  #################################################################################
  #                      Multi-Study R-Learner Oracle (both)                      # 
  #################################################################################
  
  # Training set predictors
  X_res_mat_with_int <- lapply(1:K, function(k){
    sweep(X_int[-ind_test, ], 1, (data[-ind_test, ]$A - 0.5) * p_mat[-ind_test, k], "*")
  })
  
  X_res_all_oracle <- as.data.frame(do.call("cbind", X_res_mat_with_int))
  names(X_res_all_oracle) <- sapply(1:ncol(X_res_all_oracle), function(x) paste("X.", x, sep = ""))
  
  # Training set outcomes
  Y_res_merge_oracle <- data[-ind_test, ]$Y - true_m[-ind_test]
  
  # Fit multi-study R-learner
  ms_rlearner_oracle_pred <- ms_rlearner(X = X_res_all_oracle, Y = Y_res_merge_oracle, alpha = alpha, K = K, p = p, p_mat = p_mat[ind_test, ],
                                         newX = X_int[ind_test, ])
  
  
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
    sweep(X_int[-ind_test, ], 1, (data[-ind_test, ]$A - 0.5) * p_pred_mat_train[, k], "*")
  })
  
  X_res_mat_merge_all <- as.data.frame(do.call("cbind", X_res_mat_merge))
  names(X_res_mat_merge_all) <- sapply(1:ncol(X_res_mat_merge_all), function(x) paste("X.", x, sep = ""))
  
  Y_res_merge_estimate <- data[-ind_test, ]$Y - m_estimate_merge
  ms_rlearner_mMerge_aMerge_pred <- ms_rlearner(X = X_res_mat_merge_all, 
                                                Y = Y_res_merge_estimate, 
                                                alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                                newX = X_int[ind_test, ])
  
  ####################################################################################
  #               Multi-Study R-Learner Nuisance Study-Specific                      # 
  ####################################################################################
  m_estimate_ss <- rlearner_m_hat_stack(data = data[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha, p_mat = p_pred_mat_train)
  m_estimate_avg <- m_estimate_ss$avg
  m_estimate_stack <- m_estimate_ss$stack_int
  m_estimate_ss <- m_estimate_ss$ss
  
  
  X_res_mat_avg <- lapply(1:K, function(k){
    sweep(X_int[-ind_test, ], 1, (data[-ind_test, ]$A - 0.5) * p_pred_mat_train[, k], "*")
  })
  
  X_res_mat_avg_all <- as.data.frame(do.call("cbind", X_res_mat_avg))
  names(X_res_mat_avg_all) <- sapply(1:ncol(X_res_mat_avg_all), function(x) paste("X.", x, sep = ""))
  
  # Merge m, Study-specific A
  ms_rlearner_mMerge_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_merge_estimate, alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                             newX = X_int[ind_test, ])
  
  # Avg m, Study-specific A
  Y_res_avg_estimate <- data[-ind_test, ]$Y - m_estimate_avg
  ms_rlearner_mAvg_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_avg_estimate, alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                           newX = X_int[ind_test, ])
  
  # Stack m, Study-specific A
  Y_res_stack_estimate <- data[-ind_test, ]$Y - m_estimate_stack
  ms_rlearner_mStack_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_stack_estimate, alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                             newX = X_int[ind_test, ])
  
  # Study-specific m, Study-specific A
  Y_res_ss_estimate <- data[-ind_test, ]$Y - m_estimate_ss
  ms_rlearner_mSS_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_ss_estimate, alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                          newX = X_int[ind_test, ])
  
  
  out <- data.frame(rlearner_oracle_pred = mean((true_tau - rlearner_oracle_pred)^2), 
                    rlearner_merge_pred = mean((true_tau - rlearner_merge_pred)^2), 
                    ms_rlearner_oracle_pred = mean((true_tau - ms_rlearner_oracle_pred)^2),
                    ms_rlearner_mMerge_aMerge_pred = mean((true_tau - ms_rlearner_mMerge_aMerge_pred)^2),
                    ms_rlearner_mMerge_aSS_pred = mean((true_tau - ms_rlearner_mMerge_aSS_pred)^2),
                    ms_rlearner_mAvg_aSS_pred = mean((true_tau - ms_rlearner_mAvg_aSS_pred)^2),
                    ms_rlearner_mStack_aSS_pred = mean((true_tau - ms_rlearner_mStack_aSS_pred)^2),
                    ms_rlearner_mSS_aSS_pred = mean((true_tau - ms_rlearner_mSS_aSS_pred)^2),
                    perturb_tau = perturb_tau,
                    perturb_a = perturb_a,
                    perturb_mu0 = perturb_mu0,
                    prop_overlap = prop_overlap,
                    p = p)
  return(out)
}

sim_multiple <- function(n, p, K, covariate_cols, perturb_tau, perturb_mu0, perturb_a, k_folds, alpha, prop_overlap, num_confounders, nreps){
  # Coefficients for generating mu0(x) = E[Y(0)|X=x]
  beta_mu0 <- rnorm(p + 1, mean = 0, sd = 1)
  # Select the variables to 0 out their effects (rest are the confounders common across studies)
  cols_re_mu0 <- sample(2:(p + 1), p - num_confounders + round((1 - prop_overlap) * num_confounders))
  beta_mu0[cols_re_mu0] <- 0
  # Coefficients for generating e(x) = E[A|X=x]
  beta_a <- rnorm(p + 1, mean = 0, sd = 1)
  beta_a[cols_re_mu0] <- 0
  # Coefficients for generating tau(x) = E[Y(1)-Y(0)|X=x]
  beta_tau <- rnorm(p + 1, mean = 2, sd = 1)
  beta_tau[cols_re_mu0] <- 0
  # Coefficients for generating p(k|x) = P(S=k|X=x)
  set.seed(5)
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
  # For selecting study-specific confounders
  coef_ind[[1]] <- c()
  
  # Beta_mu0_k
  for(k in 1:(K)){
    new_coefs <- rep(0, p + 1)
    coef_ind[[k]] <- sample(setdiff(cols_re_mu0, unlist(coef_ind)), round((1 - prop_overlap) * num_confounders))
    new_coefs[coef_ind[[k]]] <- rnorm(length(coef_ind[[k]]), 0, 1)
    beta_mu0_k[[k]] <- beta_mu0 + new_coefs
  }
  
  # Beta_a_k
  for(k in 1:(K)){
    new_coefs <- rep(0, p + 1)
    new_coefs[coef_ind[[k]]] <- rnorm(length(coef_ind[[k]]), 0, 1)
    beta_a_k[[k]] <- beta_a + new_coefs
  }
  
  # Beta_tau_k
  for(k in 1:(K)){
    new_coefs <- rep(0, p + 1)
    new_coefs[coef_ind[[k]]] <- rnorm(length(coef_ind[[k]]), 0, 1)
    beta_tau_k[[k]] <- beta_tau + new_coefs
  }
  
  # Indices of random effects
  cols_re_tau <- cols_re_mu0 <-cols_re_a <- vector("list", length = K)
  for(k in 1:(K)){
    cols_re_tau[[k]] <- which(beta_tau_k[[k]] != 0)
    cols_re_tau[[k]] <- cols_re_tau[[k]][-which(cols_re_tau[[k]] == 1)]
    
    cols_re_mu0[[k]] <- which(beta_mu0_k[[k]] != 0)
    cols_re_mu0[[k]] <- cols_re_mu0[[k]][-which(cols_re_mu0[[k]] == 1)]
    
    cols_re_a[[k]] <- which(beta_a_k[[k]] != 0)
    cols_re_a[[k]] <- cols_re_a[[k]][-which(cols_re_a[[k]] == 1)]
  }
  
  
  registerDoMC(cores = 48)
  results = foreach (j = 1:nreps, .combine = rbind) %dopar% {
    print(paste("Prop_overlap =", prop_overlap, "Perturb_tau =", perturb_tau, "Iteration =", j))
    sim_each(n = n, p = p, K = K, data = data, covariate_cols = covariate_cols, 
             cols_re_tau = cols_re_tau, cols_re_mu0 = cols_re_mu0, cols_re_a = cols_re_a,
             beta_tau = beta_tau, beta_mu0 = beta_mu0, beta_a = beta_a, beta_mu0_k = beta_mu0_k, beta_a_k = beta_a_k, beta_tau_k = beta_tau_k, p_mat = p_mat,
             perturb_tau = perturb_tau, perturb_mu0 = perturb_mu0, perturb_a = perturb_a,
             k_folds = k_folds, alpha = alpha, prop_overlap = prop_overlap)
  }
}

perturb_a <- 0
perturb_tau <- seq(0, 3, 0.75)
perturb_mu0 <- 0
k_folds <- 3
alpha <- 0.5

# Proportion of overlap among confounders
prop_overlap <- c(0, 0.5, 1)
# Number of confounders
num_confounders <- 8

results <- vector("list", length = length(prop_overlap))
set.seed(123)
ind_o <- 1:length(prop_overlap)
ind_t <- 1:length(perturb_tau)
for(o in ind_o){
 for(t in ind_t){
  results[[o]][[t]] = as.data.frame(sim_multiple(n = n, p = p, K = K, covariate_cols = covariate_cols, 
                                            perturb_tau = perturb_tau[t], perturb_mu0 = perturb_mu0, perturb_a = perturb_a, k_folds = k_folds, alpha = alpha,
                                            prop_overlap = prop_overlap[o], num_confounders = num_confounders,
                                            nreps = nreps))
  }
}

save(results, perturb_tau, perturb_mu0, perturb_a, prop_overlap, num_confounders, K, n, p, file = paste0("Simulation_I_p40_numConfounders8", ".RData"))


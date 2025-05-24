library(data.table)
library(nnls)
library(purrr)
library(MASS)
library(tidyverse)
library(extraDistr)
library(Rcpp)
library(doMC)
library(glmnet)
library(zetadiv)
library(zeallot)
library(grf)
library(rlearner)
library(curatedBreastData)
library(reshape2)
library(splines)
library(doParallel)
library(randomForest)

### Sample datasets from curatedOvarianData
# 
# Input:
# edat_orig: list of datasets
# ndat: number of datasets to sample
# nvar: number of predictors to sample
#
# Output:
# edat: list of data sets
#
init_data <- function(edat_orig, ndat, nvar){
  
  edat <- edat_orig
  edat <- edat[sample(1:length(edat), ndat)] # Randomize dataset order
  
  idx <- sample(1:ncol(edat[[1]]), nvar)
  for(i in 1:ndat){
    edat[[i]] <- edat[[i]][,idx]
    edat[[i]] <- as.data.frame(edat[[i]])
    colnames(edat[[i]]) <- paste0("X.", 1:nvar)
  }
  return(edat)
}

# Expit function
expit <- function(x){
  exp(x)/(1 + exp(x))
}

# Checks whether each fold contains at least 8 of each class
#
# Input
# foldid - indicies for each fold
# dataA - data set's vector of treatment assignments
#
# Output
# boolean value indicating whether there are at least an acceptable number of each class per fold
check_obs <- function(foldid, dataA, num){
  dat <- data.frame(foldid = foldid, dataA = dataA)
  uniq_fold <- unique(foldid)
  check <- vector()
  for(i in uniq_fold){
    tab <- table(dat$dataA[dat$foldid == i])
    if(length(tab[tab > num]) < 2){
      check[i] <- FALSE
    }else if(length(tab[tab > num]) == 2){
      check[i] <- TRUE
    }
  }
  if(any(check == FALSE)){
    return(FALSE)
  }else return(TRUE)
}

### Bootstrap MSPE ratio for two models
#
# Input:
# nboot: number of bootstrap iterations
# results.df: dataframe formatted like the output of sim_multi()
# col1: column of results.df corresponding to first model
# col2: column of results.df corresponding to second model
#
# Output:
# avg.mspe.merged: mean mspe for merged model
# avg.mspe.ens: mean mspe for ensemble model
# avg.mspe.ens/avg.mspe.merged: MSPE ratio comparing ensemble to merged
#
boot_ci <- function(nboot, results.df, col1, col2, seed = 1) {
  set.seed(seed)
  results.df = results.df[which(results.df[,col1] > 0 & results.df[, col2]>0),]
  avg.mspe.merged = rep(NA, nboot)
  avg.mspe.ens = rep(NA, nboot)
  for (i in 1:nboot) {
    boot.sample.ind = sample(1:nrow(results.df), nrow(results.df), replace = T)
    avg.mspe.merged[i] = mean(results.df[boot.sample.ind, col1])
    avg.mspe.ens[i] = mean(results.df[boot.sample.ind, col2])
  }
  return(list(avg.mspe.merged, avg.mspe.ens, avg.mspe.ens/avg.mspe.merged))
}


#####################################################################################################
# Calculates E[Y|X] predictions using stacking
# Input
# data - source data
# K - number of source data sets
# p - number of covariates
# covariate_cols - covariate columns
# k_folds - number of folds for cross-validation
# alpha - tuning parameter for glmnet
# p_mat - ascertainment probability matrix
#
# Output:
# list of three: 1) m_hat based on stacking with intercept, 2) m_hat based on stacking w/o intercept and 3) m_hat based on avging
#
#####################################################################################################
rlearner_m_hat_stack <- function(data, K, p, covariate_cols, k_folds, alpha, p_mat){
  for(i in 1:K){
    data[data$S == i, covariate_cols] <- as.data.frame(scale(data[data$S == i, covariate_cols], center = TRUE, scale = TRUE))
  }
  
  # Stacking
  m_hat_k <- m_hat_k_pred_on_data <- m_hat_k_pred_on_data_k <- m_hat_pred_stack_mat <- vector("list", length = K)
  lambda_min <- vector()
  for(k in 1:K){
    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(k_folds), length = length(data$A[data$S ==k])))
    
    m_hat_k[[k]] = glmnet::cv.glmnet(x = as.matrix(data[data$S == k, covariate_cols]), y = as.matrix(data$Y[data$S == k]),
                                     foldid = foldid,
                                     keep = TRUE,
                                     lambda = NULL,
                                     alpha = alpha)
    
    if(length(covariate_cols) < nrow(data)){
      # Fit cross-fitting OLS with lambda = 0
      # if(any(m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == 0) == FALSE){
      #   # If lambda = 0 did not converge, then use the smallest lambda that's not 0
      #   m_hat_k_pred_on_data[[k]] = m_hat_k[[k]]$fit.preval[,!is.na(colSums(m_hat_k[[k]]$fit.preval))][, m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == min(m_hat_k[[k]]$lambda)]
      # }else if(any(m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == 0) == TRUE){
      #   m_hat_k_pred_on_data[[k]] = m_hat_k[[k]]$fit.preval[,!is.na(colSums(m_hat_k[[k]]$fit.preval))][, m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == 0]
      # }
      m_hat_k_pred_on_data[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S <= K, covariate_cols]), s = 0)
      # lambda_min[k] = m_hat_k[[k]]$lambda[which.min(m_hat_k[[k]]$cvm[!is.na(colSums(m_hat_k[[k]]$fit.preval))])]
      # m_hat_k_pred_on_data[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S <= K, covariate_cols]), s = lambda_min[k])
    }else if(length(covariate_cols) >= nrow(data)){
      # Fit cross-fitting LASSO 
      lambda_min[k] = m_hat_k[[k]]$lambda[which.min(m_hat_k[[k]]$cvm[!is.na(colSums(m_hat_k[[k]]$fit.preval))])]
      # # Obtain cross-fitted values
      # # fit.preval contains the pre-validated fits, i.e. "this means these fits are computed with this observation and the rest of its fold omitted."
      # m_hat_k_pred_on_data[[k]] = m_hat_k[[k]]$fit.preval[,!is.na(colSums(m_hat_k[[k]]$fit.preval))][, m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == lambda_min[k]]
      m_hat_k_pred_on_data[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S <= K, covariate_cols]), s = lambda_min[k])
    }
  }
  for(j in 1:K){
    m_hat_k_pred <- vector("list", length = K)
    for(k in 1:K){
      data_pred <- as.matrix(data[which(data$S == k), c(covariate_cols)])
      if(length(covariate_cols) < nrow(data)){
        m_hat_k_pred[[k]] <- predict(m_hat_k[[j]], newx = data_pred, s = 0)
      }else if(length(covariate_cols) >= nrow(data)){
        m_hat_k_pred[[k]] <- predict(m_hat_k[[j]], newx = data_pred, s = lambda_min[j])
      }
    }
    m_hat_pred_stack_mat[[j]] <- do.call("c", m_hat_k_pred)
  }
  m_hat_pred_stack_mat <- do.call("cbind", m_hat_pred_stack_mat)
  
  # Stacked m with intercept
  m_stack_int <- nnls::nnls(cbind(rep(1, nrow(m_hat_pred_stack_mat)), m_hat_pred_stack_mat), data[data$S <= K, ]$Y)$x
  
  # Stacked m without intercept
  m_stack_noint <- nnls::nnls(m_hat_pred_stack_mat, data[data$S <= K, ]$Y)$x
  
  # Normalize so that they sum up to 1
  if(sum(m_stack_int[-1]) == 0){
    m_stack_int = m_stack_int
  }else m_stack_int[-1] <- m_stack_int[-1]/sum(m_stack_int[-1])
  
  if(sum(m_stack_noint) == 0){
    m_stack_noint = m_stack_noint
  }else m_stack_noint <- m_stack_noint/sum(m_stack_noint)
  
  # Obtain predictions on source data data and cbind them together, then use stacking weights to weight them
  m_hat_pred_on_data <- do.call("cbind", m_hat_k_pred_on_data)
  m_pred_on_data_int <- apply(m_hat_pred_on_data, 1, function(x){m_stack_int[1] + sum(m_stack_int[-1] * x)})
  m_pred_on_data_noint <- apply(m_hat_pred_on_data, 1, function(x){sum(m_stack_noint * x)})
  m_pred_on_data_avg <- apply(m_hat_pred_on_data, 1, function(x){sum(1/K * x)})
  m_pred_on_data_ss <- apply(m_hat_pred_on_data * p_mat, 1, sum)
  return(list(stack_int = m_pred_on_data_int, stack_noint = m_pred_on_data_noint, avg = m_pred_on_data_avg, ss = m_pred_on_data_ss))
}

rlearner_m_hat_stack_categorical <- function(data, K, p, covariate_cols, covariate_cols_continuous, k_folds, alpha, p_mat){
  for(i in 1:K){
    data[data$S == i, covariate_cols_continuous] <- as.data.frame(scale(data[data$S == i, covariate_cols_continuous], center = TRUE, scale = TRUE))
  }
  
  # Stacking
  m_hat_k <- m_hat_k_pred_on_data <- m_hat_pred_stack_mat <- m_hat_k_pred_on_data_k <- vector("list", length = K)
  lambda_min <- vector()
  for(k in 1:K){
    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(k_folds), length = length(data$A[data$S ==k])))
    
    m_hat_k[[k]] = glmnet::cv.glmnet(x = as.matrix(data[data$S == k, covariate_cols]), y = as.matrix(data$Y[data$S == k]),
                                     foldid = foldid,
                                     keep = TRUE,
                                     lambda = NULL,
                                     alpha = alpha)
    
    if(length(covariate_cols) < nrow(data)){
      # Fit cross-fitting OLS with lambda = 0
      # if(any(m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == 0) == FALSE){
      #   # If lambda = 0 did not converge, then use the smallest lambda that's not 0
      #   m_hat_k_pred_on_data[[k]] = m_hat_k[[k]]$fit.preval[,!is.na(colSums(m_hat_k[[k]]$fit.preval))][, m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == min(m_hat_k[[k]]$lambda)]
      # }else if(any(m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == 0) == TRUE){
      #   m_hat_k_pred_on_data[[k]] = m_hat_k[[k]]$fit.preval[,!is.na(colSums(m_hat_k[[k]]$fit.preval))][, m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == 0]
      # }
      m_hat_k_pred_on_data[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S <= K, covariate_cols]), s = 0)
      m_hat_k_pred_on_data_k[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S == k, covariate_cols]), s = 0)
    }else if(length(covariate_cols) >= nrow(data)){
      # Fit cross-fitting LASSO 
      lambda_min[k] = m_hat_k[[k]]$lambda[which.min(m_hat_k[[k]]$cvm[!is.na(colSums(m_hat_k[[k]]$fit.preval))])]
      # # Obtain cross-fitted values
      # # fit.preval contains the pre-validated fits, i.e. "this means these fits are computed with this observation and the rest of its fold omitted."
      # m_hat_k_pred_on_data[[k]] = m_hat_k[[k]]$fit.preval[,!is.na(colSums(m_hat_k[[k]]$fit.preval))][, m_hat_k[[k]]$lambda[!is.na(colSums(m_hat_k[[k]]$fit.preval))] == lambda_min[k]]
      m_hat_k_pred_on_data[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S <= K, covariate_cols]), s = lambda_min[k])
      m_hat_k_pred_on_data_k[[k]] = predict(m_hat_k[[k]], newx = as.matrix(data[data$S == k, covariate_cols]), s = lambda_min[k])
    }
  }
  for(j in 1:K){
    m_hat_k_pred <- vector("list", length = K)
    for(k in 1:K){
      data_pred <- as.matrix(data[which(data$S == k), c(covariate_cols)])
      if(length(covariate_cols) < nrow(data)){
        m_hat_k_pred[[k]] <- predict(m_hat_k[[j]], newx = data_pred, s = 0)
      }else if(length(covariate_cols) >= nrow(data)){
        m_hat_k_pred[[k]] <- predict(m_hat_k[[j]], newx = data_pred, s = lambda_min[j])
      }
    }
    m_hat_pred_stack_mat[[j]] <- do.call("c", m_hat_k_pred)
  }
  m_hat_pred_stack_mat <- do.call("cbind", m_hat_pred_stack_mat)
  
  # Stacked m with intercept
  m_stack_int <- nnls::nnls(cbind(rep(1, nrow(m_hat_pred_stack_mat)), m_hat_pred_stack_mat), data[data$S <= K, ]$Y)$x
  
  # Stacked m without intercept
  m_stack_noint <- nnls::nnls(m_hat_pred_stack_mat, data[data$S <= K, ]$Y)$x
  
  # Normalize so that they sum up to 1
  if(sum(m_stack_int[-1]) == 0){
    m_stack_int = m_stack_int
  }else m_stack_int[-1] <- m_stack_int[-1]/sum(m_stack_int[-1])
  
  if(sum(m_stack_noint) == 0){
    m_stack_noint = m_stack_noint
  }else m_stack_noint <- m_stack_noint/sum(m_stack_noint)
  
  # Obtain predictions on source data data and cbind them together, then use stacking weights to weight them
  m_hat_pred_on_data <- do.call("cbind", m_hat_k_pred_on_data)
  m_pred_on_data_int <- apply(m_hat_pred_on_data, 1, function(x){m_stack_int[1] + sum(m_stack_int[-1] * x)})
  m_pred_on_data_noint <- apply(m_hat_pred_on_data, 1, function(x){sum(m_stack_noint * x)})
  m_pred_on_data_avg <- apply(m_hat_pred_on_data, 1, function(x){sum(1/K * x)})
  m_pred_on_data_ss <- apply(m_hat_pred_on_data * p_mat, 1, sum)
  return(list(stack_int = m_pred_on_data_int, stack_noint = m_pred_on_data_noint, avg = m_pred_on_data_avg, ss = m_pred_on_data_ss, ss_by_study = m_hat_k_pred_on_data_k))
}

#####################################################################################################
# Calculates E[Y|X] predictions using merging
# Input
# data - source data
# K - number of source data sets
# p - number of covariates
# covariate_cols - covariate columns
#
# Output:
# m_hat predictions
#
#####################################################################################################
rlearner_m_hat_merge <- function(data, K, p, covariate_cols, k_folds, alpha){
  # Center and scale merged data
  data[data$S <= K, covariate_cols] <- as.data.frame(scale(data[data$S <= K, covariate_cols], center = TRUE, scale = TRUE))
  
  # Divide the data into folds for cross fitting
  foldid = sample(rep(seq(k_folds), length = length(data$A)))
  
  m_hat_mod = glmnet::cv.glmnet(x = as.matrix(data[data$S <= K, covariate_cols]), y = as.matrix(data$Y[data$S <= K]),
                                foldid = foldid,
                                keep = TRUE,
                                alpha = alpha)
  
  # Fit cross-fitting LASSO
  # Obtain the lambda based on lowest cvm (cross-validated error)
  lambda_min = m_hat_mod$lambda[which.min(m_hat_mod$cvm[!is.na(colSums(m_hat_mod$fit.preval))])]
    
  # Obtain cross-fitted values
  # fit.preval contains the pre-validated fits, i.e. "this means these fits are computed with this observation and the rest of its fold omitted."
  m_hat_mod_pred = m_hat_mod$fit.preval[,!is.na(colSums(m_hat_mod$fit.preval))][, m_hat_mod$lambda[!is.na(colSums(m_hat_mod$fit.preval))] == lambda_min]
  
  return(m_hat_mod_pred)
}

#####################################################################################################
# Calculates E[A|X] predictions using merging
# Input
# data - source data
# K - number of source data sets
# p - number of covariates
# covariate_cols - covariate columns
#
# Output:
# m_hat predictions
#
#####################################################################################################
rlearner_a_hat_merge <- function(data, K, p, covariate_cols, k_folds, alpha){
  # Center and scale merged data
  data[data$S <= K, covariate_cols] <- as.data.frame(scale(data[data$S <= K, covariate_cols], center = TRUE, scale = TRUE))
  
  # Divide the data into folds for cross fitting
  foldid = sample(rep(seq(k_folds), length = length(data$A)))
  
  # Fit cross-fitting LASSO with lambda = 0
  a_hat_mod = glmnet::cv.glmnet(x = as.matrix(data[data$S <= K, covariate_cols]), y = as.matrix(data$A[data$S <= K]),
                                  foldid = foldid,
                                  family="binomial",
                                  type.measure="deviance",
                                  keep = TRUE,
                                  lambda = NULL,
                                  alpha = alpha)
    
  # Obtain the lambda based on lowest cvm (cross-validated error)
  a_lambda_min = a_hat_mod$lambda[which.min(a_hat_mod$cvm[!is.na(colSums(a_hat_mod$fit.preval))])]
    
  # theta_hat are the cross-fitted values on the linear scale
  theta_hat = a_hat_mod$fit.preval[,!is.na(colSums(a_hat_mod$fit.preval))][, a_hat_mod$lambda[!is.na(colSums(a_hat_mod$fit.preval))] == a_lambda_min]
    
  # expit(theta_hat) to get everything on the logistic scale
  a_hat_mod_pred = 1/(1 + exp(-theta_hat))
  
  return(a_hat_mod_pred)
}

rlearner_a_hat_merge_categorical <- function(data, K, p, covariate_cols, covariate_cols_continuous, k_folds, alpha){
  # Center and scale merged data
  data[data$S <= K, covariate_cols_continuous] <- as.data.frame(scale(data[data$S <= K, covariate_cols_continuous], center = TRUE, scale = TRUE))
  
  # Divide the data into folds for cross fitting
  foldid = sample(rep(seq(k_folds), length = length(data$A)))
  
  # Fit cross-fitting LASSO with lambda = 0
  a_hat_mod = glmnet::cv.glmnet(x = as.matrix(data[data$S <= K, covariate_cols]), y = as.matrix(data$A[data$S <= K]),
                                foldid = foldid,
                                family="binomial",
                                type.measure="deviance",
                                keep = TRUE,
                                lambda = NULL,
                                alpha = alpha)
  
  # Obtain the lambda based on lowest cvm (cross-validated error)
  a_lambda_min = a_hat_mod$lambda[which.min(a_hat_mod$cvm[!is.na(colSums(a_hat_mod$fit.preval))])]
  
  # theta_hat are the cross-fitted values on the linear scale
  theta_hat = a_hat_mod$fit.preval[,!is.na(colSums(a_hat_mod$fit.preval))][, a_hat_mod$lambda[!is.na(colSums(a_hat_mod$fit.preval))] == a_lambda_min]
  
  # expit(theta_hat) to get everything on the logistic scale
  a_hat_mod_pred = 1/(1 + exp(-theta_hat))
  
  return(a_hat_mod_pred)
}



#####################################################################################################
# Calculates E(A|X) predictions using stacking
# Input
# data - source data
# K - number of source data sets
# p - number of covariates
# covariate_cols - covariate columns
#
# Output:
# list of three: 1) a_hat based on stacking with intercept, 2) a_hat based on stacking w/o intercept and 3) a_hat based on avging
#
#####################################################################################################
rlearner_a_hat_stack <- function(data, K, p, covariate_cols, k_folds, alpha){
  for(i in 1:K){
    data[data$S == i, covariate_cols] <- as.data.frame(scale(data[data$S == i, covariate_cols], center = TRUE, scale = TRUE))
  }
  
  # Stacking
  stack_p_mod <- NA
  while(is.na(stack_p_mod[1])){
    a_hat_k <- a_hat_k_pred_on_data <- a_hat_pred_stack_mat <- foldid <- vector("list", length = K)
    lambda_min <- vector()
    for(k in 1:K){
      foldid[[k]] <- NA
      check <- FALSE
      repeat{
        foldid[[k]] = sample(rep(seq(k_folds), length = length(data$A[data$S == k])))
        check <- check_obs(foldid = foldid[[k]], dataA = data$A[data$S ==k], num = 3)
        if(check){
          break
        }
      }
      
        # Fit cross-fitting glm LASSO 
        a_hat_k[[k]] = glmnet::cv.glmnet(x = as.matrix(data[data$S == k, covariate_cols]), y = as.matrix(data$A[data$S == k]),
                                         foldid = foldid[[k]],
                                         family="binomial",
                                         type.measure="deviance",
                                         keep = TRUE,
                                         lambda = NULL,
                                         alpha = alpha)
        
        # Obtain the lambda based on lowest cvm (cross-validated error)
        lambda_min[k] = a_hat_k[[k]]$lambda[which.min(a_hat_k[[k]]$cvm[!is.na(colSums(a_hat_k[[k]]$fit.preval))])]
        
        # Obtain cross-fitted values on the linear scale
        # a_hat_k_pred_on_data[[k]] = a_hat_k[[k]]$fit.preval[,!is.na(colSums(a_hat_k[[k]]$fit.preval))][, a_hat_k[[k]]$lambda[!is.na(colSums(a_hat_k[[k]]$fit.preval))] == lambda_min[k]]
        a_hat_k_pred_on_data[[k]] = predict(a_hat_k[[k]], newx = as.matrix(data[data$S <= K, covariate_cols]), s = lambda_min[k])
      }
    
    for(j in 1:K){
      a_hat_k_pred <- vector("list", length = K)
      for(k in 1:K){
        data_pred <- as.matrix(data[which(data$S == k), c(covariate_cols)])
        if(length(covariate_cols) < nrow(data)){
          a_hat_k_pred[[k]] <- predict(a_hat_k[[j]], newx = data_pred, s = 0)
        }else if(length(covariate_cols) >= nrow(data)){
          a_hat_k_pred[[k]] <- predict(a_hat_k[[j]], newx = data_pred, s = lambda_min[j])
        }
        
      }
      a_hat_pred_stack_mat[[j]] <- do.call("c", a_hat_k_pred)
    }
    a_hat_pred_stack_mat <- do.call("cbind", a_hat_pred_stack_mat)
    
    # Fit stacking classifier (stacked on linear scale)
    stack_p_mod <- tryCatch({
      glmnet(x = a_hat_pred_stack_mat, y = data$A[data$S <= K], family = "binomial", type.measure="deviance", lower.limits = 0, alpha = 0, lambda = 0)
    },
    error = function(cond){
      return(NA)
    },
    warning = function(cond){
      return(NA)
    })
  }
  
  if(sum(coefficients(stack_p_mod)[-1]) > 0){
    stack_p_mod_coefs <- c(coefficients(stack_p_mod)[1], coefficients(stack_p_mod)[-1]/(sum(coefficients(stack_p_mod)[-1])))
  }else stack_p_mod_coefs <- coefficients(stack_p_mod)
  
  # Obtain predictions on test data and cbind them together, then use stacking weights to weight them (on the linear scale first) then expit the results
  a_hat_pred_on_data <- do.call("cbind", a_hat_k_pred_on_data)
  p_pred_on_data_int <- expit(apply(a_hat_pred_on_data, 1, function(x){stack_p_mod_coefs[1] + sum(stack_p_mod_coefs[-1] * x)}))
  p_pred_on_data_avg <- expit(apply(a_hat_pred_on_data, 1, function(x){sum(1/K * x)}))
  
  return(list(stack_int = p_pred_on_data_int, avg = p_pred_on_data_avg))
}


ms_rlearner <- function(X, Y, alpha, K, p, p_mat, newX){
  ms_rlearner <- glmnet::cv.glmnet(x = as.matrix(X),
                                   y = as.matrix(Y),
                                   alpha = alpha,
                                   lambda = NULL,
                                   standardize = TRUE,
                                   intercept = FALSE,
                                   standardize.response = TRUE)
  
  # Remove the first entry that corresponds to the intercept (which is empty)
  tau_beta <- as.vector(t(coef(ms_rlearner, s = "lambda.min")[-1]))
  # p + 1 coefficients per study (p covariates and 1 intercept)
  tau_beta_k <- split(tau_beta, ceiling(seq_along(tau_beta)/(p + 1)))
  
  ms_rlearner_pred_list <- vector("list", length = K)
  for(k in 1:K){
    ms_rlearner_pred_list[[k]] <- newX %*% as.matrix(tau_beta_k[[k]]) * p_mat[, k]
  }
  ms_rlearner_pred_mat <- do.call("cbind", ms_rlearner_pred_list)
  ms_rlearner_pred <- apply(ms_rlearner_pred_mat, 1, sum)
  return(ms_rlearner_pred)
}

make_basis <- function(k, p) replace(numeric(p), k, 1)

generate_tau_k <- function(X_row, k, offset_vec) {
  X_perturbed <- X_row + offset_vec
  (sin((2) * X_perturbed[1]) + log(abs(X_perturbed[2] * 1)) +
    cos(X_perturbed[3]) * X_perturbed[4]^2)/33
}

generate_spline_basis <- function(X_df, df_spline = 3) {
  spline_cols <- lapply(seq_along(X_df), function(j) {
    col_j <- X_df[[j]]
    if (length(unique(col_j)) < df_spline) {
      # Reduce df if not enough unique values
      df_effective <- max(1, length(unique(col_j)) - 1)
    } else {
      df_effective <- df_spline
    }
    bs(col_j, df = df_effective)
  })
  spline_mat <- do.call(cbind, spline_cols)
  colnames(spline_mat) <- paste0("spline_", seq_len(ncol(spline_mat)))
  return(as.matrix(spline_mat))
}

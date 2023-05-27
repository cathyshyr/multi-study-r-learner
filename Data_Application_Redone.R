source("Utils.R")
library(curatedBreastData)
data(clinicalData)
data(curatedBreastDataExprSetList)

# Study 21997
study_GSE21997 <- clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 21997, ]
PID1 <- study_GSE21997$patient_ID

# Study 25065
study_GSE25065 <- clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 25065, ]
PID2 <- study_GSE25065$patient_ID

# Outcome variable: pathological complete response (pCR)
studies <- list(study_GSE21997, study_GSE25065)

# K = number of studies
K <- length(studies)

# alpha for glmnet
alpha <- 0

load("proc_curatedBreastDataExprSetList.RData")

study_IDs <- sapply(1:length(curatedBreastDataExprSetList), function(x) unique(curatedBreastDataExprSetList[[x]]$study_ID.x))
study_IDs_selected1 <- which(study_IDs %in% c(21997))
study_IDs_selected2 <- which(study_IDs %in% c(25065))

esets1 <- do.call("cbind", lapply(study_IDs_selected1, function(x) exprs(proc_curatedBreastDataExprSetList[[x]])))
esets2 <- do.call("cbind", lapply(study_IDs_selected2, function(x) exprs(proc_curatedBreastDataExprSetList[[x]])))
esets <- list(esets1, esets2)

esets_PID1 <- as.numeric(colnames(esets1))
esets_PID2 <- as.numeric(colnames(esets2))

# Take the intersection of genes
cn <- lapply(esets, rownames)
cn_int <- Reduce(intersect, cn)
# Retain the oncotype DX genes
oncotypeDx <- c("Ki67", "STK15", "SURVIVIN", "CCNB1", "MYBL2", "MMP11", "CTSL2",
                "GRB2", "HER2", "ER", "PGR", "BCL2", "SCUBE2", "GSTM1",
                "CD68", "BAG1", "ACTB", "GAPDH", "RPLPO", "GUS", "TFRC")
oncotypeDx_ID <- which(cn_int %in% oncotypeDx)
oncotypeDx_genes <- cn_int[oncotypeDx_ID]
# Randomly select up to 100 genes including oncotype DX
set.seed(3)
cn_int_100 <- sample(setdiff(cn_int, oncotypeDx_genes), 100 - length(oncotypeDx_ID))
cn_int_final <- c(oncotypeDx_genes, cn_int_100)

# Covariates: age, ER_preTrt (1 if positive, 0 if not), PR_preTrt (1 if positive, 0 if not), HER2_preTrt (1 if positive, 0 if not)
#             hist_grade (only study_GSE21997 and study_GSE25065), path (only study_GSE21997 and study_GSE20194)
data <- vector('list', length = length(studies))
for(i in 1:length(studies)){
  tmp <- studies[[i]]
  data[[i]] <- tmp[ , c("age", "pCR", "ER_preTrt", "PR_preTrt", "HER2_preTrt", "doxorubicin", "hist_grade")]
  if(i == 1){
    data[[i]] <- cbind(data[[i]], t(esets1[cn_int, as.character(tmp$patient_ID)]))
  }else if(i == 2){
    data[[i]] <- cbind(data[[i]], t(esets2[cn_int, as.character(tmp$patient_ID)]))
  }
  data[[i]] <- data[[i]][complete.cases(data[[i]]), ]
  data[[i]]$S <- i
}
data_all <- do.call("rbind", data)

# Change to numeric variable
data_all$doxorubicin <- as.numeric(levels(data_all$doxorubicin))[data_all$doxorubicin ]
data_all$hist_grade <- as.numeric(data_all$hist_grade)

# Divide data into training and test
set.seed(3)
ind_test <- sample(nrow(data_all), 100)

# Names of covariates
covariate_cols <- c("age", "ER_preTrt", "PR_preTrt", "HER2_preTrt", "hist_grade", cn_int_final)
p <- length(covariate_cols)

set.seed(3)
beta_tau1 <- rnorm(p + 1, 0, 0.3)
beta_tau0 <- rnorm(p + 1, 0, 0.3)
ind0 <- 15:106
beta_tau1[ind0] <- 0
beta_tau0[ind0] <- 0

data_application <- function(data_all, K, p, alpha, beta_tau1, beta_tau0, perturb_tau, covariate_cols, ind_test, ind0){
  beta_tau_k_A1 <- beta_tau_k_A0 <- tau_k <- pA1 <- pA0 <- vector("list", length = K)
  X_int <- cbind(1, data_all[, covariate_cols])
  for(k in 1:K){
    # Tau
    sigma_tau <- rep(perturb_tau, p + 1)
    if(k == 1){
      sigma_tau[c(1, ind0)] <- 0
      beta_tau_k_A1[[k]] <- beta_tau1 + rnorm(p + 1, 0, sigma_tau)
      beta_tau_k_A0[[k]] <- beta_tau0 + rnorm(p + 1, 0, sigma_tau)
    }else if(k == 2){
      # Set HER2+'s random effect to 0 for study 2 because it only has HER2- patients
      sigma_tau[c(1, 5, ind0)] <- 0
      beta_tau1[5] <- 0
      beta_tau0[5] <- 0
      beta_tau_k_A1[[k]] <- beta_tau1 + rnorm(p + 1, 0, sigma_tau)
      beta_tau_k_A0[[k]] <- beta_tau0 + rnorm(p + 1, 0, sigma_tau)
    }
   
    pA1[[k]] <- expit(as.matrix(X_int[data_all$S == k, ]) %*% beta_tau_k_A1[[k]]/50)
    pA0[[k]] <- expit(as.matrix(X_int[data_all$S == k, ]) %*% beta_tau_k_A0[[k]]/50)
    tau_k[[k]] <-  pA1[[k]] - pA0[[k]]
    data_all$true_tau[data_all$S == k] <- tau_k[[k]]
    data_all$Y1[data_all$S == k] <- as.numeric(rbernoulli(n = length(pA1[[k]]), p = pA1[[k]]))
    data_all$Y0[data_all$S == k] <- as.numeric(rbernoulli(n = length(pA0[[k]]), p = pA0[[k]]))
    data_all$Y[data_all$S == k] <- ifelse(data_all$doxorubicin[data_all$S == k] == 1, 
                                          data_all$Y1[data_all$S == k], 
                                          data_all$Y0[data_all$S == k])
  }
  
  # Divide data into training and test
  ind_test <- sample(nrow(data_all), 30)
  check <- 0
  while(check == 0){
    data_test <- data_all[ind_test, ]
    data_train <- data_all[-ind_test, ]
    test_true_tau <- data_test$true_tau
    # Check if study 2 has all 1 or 0s
    if(all(data_train$Y[data_train$S == 2] == 1) | all(data_train$Y[data_train$S == 2] == 0)){
      check = 0
    }else{
      check = 1
    }
  }
  
  # Estimate m 
  m_estimate_merge <- glmnet::cv.glmnet(x = data.matrix(data_train[, covariate_cols]),
                                        y = as.matrix(data_train$Y),
                                        alpha = alpha,
                                        lambda = NULL,
                                        standardize = TRUE,
                                        keep = TRUE,
                                        family = "binomial",
                                        type.measure = "deviance")
  m_hat_lambda_min = m_estimate_merge$lambda[which.min(m_estimate_merge$cvm[!is.na(colSums(m_estimate_merge$fit.preval))])]
  m_estimate_merge_pred <- predict(m_estimate_merge, newx = data.matrix(data_train[, covariate_cols]), s = m_hat_lambda_min, type = "response")
  
  # Fit multi-study R-learner on training set
  # Estimate ascertainment probability
  p_mod <- glmnet::cv.glmnet(x = data.matrix(data_train[, covariate_cols]),
                             y = as.matrix(data_train$S),
                             alpha = alpha,
                             lambda = NULL,
                             standardize = TRUE,
                             keep = TRUE,
                             family = "multinomial")
  p_hat_lambda_min = p_mod$lambda[which.min(p_mod$cvm[!is.na(colSums(p_mod$fit.preval))])]
  p_pred_mat_train <- as.data.frame(predict(p_mod, newx = data.matrix(data_train[, covariate_cols]), s = p_hat_lambda_min, type = "response"))
  p_pred_mat_test <- as.data.frame(predict(p_mod, newx = data.matrix(data_test[, covariate_cols]), s = p_hat_lambda_min, type = "response"))
  
  # Propensity score model e_k(x) for study 2
  ind2 <- which(data_train$S == 2)
  a_hat_k <- glmnet::cv.glmnet(x = data.matrix(data_train[ind2, covariate_cols]),
                                      y = as.matrix(data_train$doxorubicin[ind2]),
                                      alpha = alpha,
                                      lambda = NULL,
                                      standardize = TRUE,
                                      keep = TRUE,
                                      family = "binomial",
                                      type.measure = "deviance")
  a_hat_k_lambda_min = a_hat_k$lambda[which.min(a_hat_k$cvm[!is.na(colSums(a_hat_k$fit.preval))])]

  
  X_int_train <- cbind(1, data_train[, covariate_cols])

  X_res_mat_avg1 <- sweep(X_int_train, 1, (data_train$doxorubicin - 0.5) * p_pred_mat_train[, 1], "*")
  X_res_mat_avg2 <- sweep(X_int_train, 1, (data_train$doxorubicin - predict(a_hat_k, newx = data.matrix(data_train[, covariate_cols]), s = a_hat_k_lambda_min, type = "response")) * p_pred_mat_train[, 2], "*")
  X_res_mat_avg <- list(X_res_mat_avg1, X_res_mat_avg2)
  
  X_res_mat_avg_all <- as.data.frame(do.call("cbind", X_res_mat_avg))
  names(X_res_mat_avg_all) <- sapply(1:ncol(X_res_mat_avg_all), function(x) paste("X.", x, sep = ""))
  
  Y_res_merge_estimate <- data_train$Y - m_estimate_merge_pred
  
  # Estimate tau
  ms_rlearner <- glmnet::cv.glmnet(x = as.matrix(X_res_mat_avg_all),
                                   y = as.matrix(Y_res_merge_estimate),
                                   alpha = alpha,
                                   lambda = NULL,
                                   standardize = TRUE,
                                   intercept = FALSE,
                                   standardize.response = TRUE)
  tau_beta <- as.vector(t(coef(ms_rlearner, s = "lambda.min")[-1]))
  tau_beta_k <- split(tau_beta, ceiling(seq_along(tau_beta)/(p + 1)))
  
  X_int_test <- cbind(1, as.matrix(data_test[, covariate_cols]))
  ms_rlearner_pred_list <- vector("list", length = K)
  for(k in 1:K){
    ms_rlearner_pred_list[[k]] <- as.matrix(X_int_test) %*% as.matrix(tau_beta_k[[k]]) * p_pred_mat_test[, k]
  }
  ms_rlearner_pred_mat <- do.call("cbind", ms_rlearner_pred_list)
  ms_rlearner_pred <- 2 * expit(apply(ms_rlearner_pred_mat, 1, sum)) - 1
  
  
  # Fit R-learner on training set and test on test set
  rlearner_mod <- rlasso(x = data.matrix(data_train[, covariate_cols]),
                         w = data_train$doxorubicin,
                         y = data_train$Y,
                         lambda_choice = "lambda.min",
                         alpha = alpha, 
                         m_hat = m_estimate_merge_pred, 
                         p_hat = c(rep(0.5, length(which(data_train$S == 1))),
                                  predict(a_hat_k, newx = data.matrix(data_train[which(data_train$S == 2), covariate_cols]), s = a_hat_k_lambda_min, type = "response")))
  
  rlearner_pred <-  2 * expit(predict(rlearner_mod, newx = as.matrix(data_test[, covariate_cols]))) - 1
  
  out <- data.frame(rlearner =  mean((test_true_tau - rlearner_pred)^2),
                    ms_rlearner = mean((test_true_tau - ms_rlearner_pred)^2),
                    test_true_tau = test_true_tau,
                    rlearner_pred = rlearner_pred, 
                    ms_rlearner_pred = ms_rlearner_pred)
  return(out)
}

data_application_multiple <- function(data_all, K, p, alpha, beta_tau1, beta_tau0, perturb_tau, covariate_cols, ind_test, ind0, nreps){
  registerDoMC(cores = 48)
  results = foreach (j = 1:nreps, .combine = rbind) %dopar% {
    print(paste("Perturb_tau =", perturb_tau, "Iteration =", j))
    data_application(data_all = data_all, K = K, p = p, alpha = alpha, beta_tau1 = beta_tau1, beta_tau0 = beta_tau0, perturb_tau = perturb_tau, 
                     covariate_cols = covariate_cols, ind_test = ind_test, ind0 = ind0)
  }
}


nreps <- 500
results <- vector("list", length = length(perturb_tau))
perturb_tau <- seq(0, 2, 0.5)
ind_t <- 1:length(perturb_tau)
for(t in ind_t){
  results[[t]] = as.data.frame(data_application_multiple(data_all = data_all, K = K, p = p, alpha = alpha, beta_tau1 = beta_tau1, 
                                                         beta_tau0 = beta_tau0, perturb_tau = perturb_tau[t], 
                                                         covariate_cols = covariate_cols, ind_test = ind_test, ind0 = ind0,
                                                         nreps = nreps))
}


save(results, perturb_tau, data_all, K, p, beta_tau1, beta_tau0, file = paste0("data_application_0806", ".RData"))

q(save = "no")

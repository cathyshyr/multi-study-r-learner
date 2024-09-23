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

# Study 20194
study_GSE20194 <- clinicalData$clinicalTable[clinicalData$clinicalTable$study_ID == 20194, ]
PID3 <- study_GSE20194$patient_ID

# Outcome variable: pathological complete response (pCR)
studies <- list(study_GSE21997, study_GSE25065, study_GSE20194)

# K = number of studies
K <- length(studies)

# alpha for glmnet
alpha <- 0

load("./proc_curatedBreastDataExprSetList.RData")

study_IDs <- sapply(1:length(curatedBreastDataExprSetList), function(x) unique(curatedBreastDataExprSetList[[x]]$study_ID.x))
study_IDs_selected1 <- which(study_IDs %in% c(21997))
study_IDs_selected2 <- which(study_IDs %in% c(25065))
study_IDs_selected3 <- which(study_IDs %in% c(20194))

esets1 <- do.call("cbind", lapply(study_IDs_selected1, function(x) exprs(proc_curatedBreastDataExprSetList[[x]])))
esets2 <- do.call("cbind", lapply(study_IDs_selected2, function(x) exprs(proc_curatedBreastDataExprSetList[[x]])))
esets3 <- do.call("cbind", lapply(study_IDs_selected3, function(x) exprs(proc_curatedBreastDataExprSetList[[x]])))
esets <- list(esets1, esets2, esets3)

esets_PID1 <- as.numeric(colnames(esets1))
esets_PID2 <- as.numeric(colnames(esets2))
esets_PID3 <- as.numeric(colnames(esets3))

# Take the intersection of genes
cn <- lapply(esets, rownames)
cn_int <- Reduce(intersect, cn)
# Retain the oncotype DX genes
oncotypeDx_andOthers <- c("Ki67", "STK15", "SURVIVIN", "CCNB1", "MYBL2", "MMP11", "CTSL2",
                          "GRB2", "HER2", "ER", "PGR", "BCL2", "SCUBE2", "GSTM1",
                          "CD68", "BAG1", "ACTB", "GAPDH", "RPLPO", "GUS", "TFRC",
                          ####### Other genes that have been implicated in modifying the treatment effect #####,
                          "TP53", "BRCA1", "BRCA2", "ERBB2", "ESR1", "PGR", "CYP2D6", 
                          "ABCB1", "TOP2A", "MKI67", "BCL2", "CDKN1A", "ATM", "PARP1", 
                          "EGFR", "PTEN", "MYC", "VEGFA", "CCND1", "FGFR1", "PIK3CA", 
                          "AKT1", "KRAS", "HRAS", "NRAS", "NOTCH1", "STAT3", "JAK2", 
                          "SRC", "EP300", "CREBBP", "MDM2", "FOXO3", "GATA3", "BRAF", 
                          "TGFBR2", "SMAD4", "RUNX1", "ETS1", "PDGFRB", "CDH1", "ERBB3", 
                          "ERRB4", "SHH", "GLI1", "SMO", "PTCH1", "RB1", "CDK4", "CDK6", 
                          "MCL1", "BAX", "BCLXL", "BIM", "NOXA", "PUMA", "GSK3B", "MTOR", 
                          "TSC1", "TSC2", "PIK3R1", "PTPN11", "SOS1", "MAPK1", "MAPK3", 
                          "JUN", "FOS", "SP1", "NF1", "NF2", "SUFU", "HHIP", "PTCH2", 
                          "WNT1", "CTNNB1", "AXIN1", "APC", "LRP6", "GSK3A", "TCF7L2", 
                          "CDH2", "CDH3", "TWIST1", "SNAI1", "ZEB1", "VIM", "MMP9", 
                          "MMP2", "ITGB1", "COL1A1", "FN1", "ABL1", "ABL2", "AKT2", 
                          "ALK", "AP1G1", "APC", "AR", "ARID1A", "ARID1B", "ATM", 
                          "ATR", "AURKA", "AURKB", "AXL", "BCL2A1", "BCL2L12", "BCL7A", 
                          "BCL9", "BCR", "BIRC2", "BIRC3", "BIRC5", "BLM", "BRAF", 
                          "BRIP1", "BTK", "C11orf30", "C16orf3", "C18orf1", "C19orf12", 
                          "CBL", "CCNE1", "CD22", "CD79A", "CD79B", "CDC42", "CDH11", 
                          "CDK2", "CDK3", "CDK5", "CDK7", "CDKN2A", "CDKN2B", "CEBPA", 
                          "CHEK1", "CHEK2", "CHUK", "CIC", "CLDN1", "CLDN3", "CLDN4", 
                          "CLDN7", "CLPTM1L", "CREB1", "CSF1R", "CSNK2A1", "CSNK2B", 
                          "CTBP1", "CTNNB1", "CUL3", "CXCL12", "CXCR4", "CYLD", "DAXX", 
                          "DDR1", "DDR2", "DUSP4", "DUSP6", "E2F1", "E2F3", "E2F4", 
                          "E2F5", "E2F6", "E2F7", "E2F8", "EED", "EIF4A2", "EIF4EBP1", 
                          "ELF3", "EPHA2", "EPHB4", "ERBB4", "ERCC1", "ERCC2", "ERCC3", 
                          "ERCC4", "ERCC5", "ERCC6", "EZH2", "FAT1", "FAT2", "FAT3", 
                          "FAT4", "FBLN5", "FEN1", "FGF1", "FGF2", "FGFR2", "FGFR3", 
                          "FGFR4", "FLT1", "FLT3", "FLT4", "FOXA1", "FOXA2", "FOXM1", 
                          "FOXP1", "FSTL3", "FZR1", "GATA1", "GATA2", "GATA4", "GATA6", 
                          "GJA1", "GNAS", "GNG4", "GRB2", "GRIN2A", "GRM1", "GSK3B", 
                          "HGF", "HIF1A", "HSP90AA1", "HSP90AB1", "ICAM1", "IGF1R", 
                          "IGF2", "IKZF1", "IL7R", "INPP4B", "IRAK1", "IRAK2", "IRS2", 
                          "JAK1", "JAK2", "JUNB", "KDM5A", "KDM5B", "KDM6A", "KDR", 
                          "KIT", "KMT2A", "KMT2C", "KMT2D", "KRAS", "LAMP3", "LATS1", 
                          "LATS2", "LCK", "LEF1", "LIMK1", "LMO2", "LOXL2", "LRIG1", 
                          "LRP1", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K11", 
                          "MAP3K4", "MAP3K5", "MAPK14", "MAPK8", "MAPK9", 
                          "MDM2", "MDM4", "MEF2B", "MEN1", "MET", "MIR21", "MITF", 
                          "MLH1", "MLH3", "MMP11", "MMP13", "MMP14", "MMP15", "MMP16", 
                          "MMP19", "MMP3", "MMP7", "MMP8", "MPL")
oncotypeDx_ID <- which(cn_int %in% oncotypeDx_andOthers)
oncotypeDx_genes <- cn_int[oncotypeDx_ID]
cn_int_final <- c(oncotypeDx_genes)

# Covariates: age, ER_preTrt (1 if positive, 0 if not), PR_preTrt (1 if positive, 0 if not), HER2_preTrt (1 if positive, 0 if not)
data <- vector('list', length = length(studies))
for(i in 1:length(studies)){
  tmp <- studies[[i]]
  data[[i]] <- tmp[ , c("age", "pCR", "ER_preTrt", "PR_preTrt", "HER2_preTrt", "doxorubicin")]
  if(i == 1){
    data[[i]] <- cbind(data[[i]], t(esets1[cn_int, as.character(tmp$patient_ID)]))
  }else if(i == 2){
    data[[i]] <- cbind(data[[i]], t(esets2[cn_int, as.character(tmp$patient_ID)]))
  }else if(i == 3){
    data[[i]] <- cbind(data[[i]], t(esets3[cn_int, as.character(tmp$patient_ID)]))
  }
  data[[i]] <- data[[i]][complete.cases(data[[i]]), ]
  data[[i]]$S <- i
}
data_all <- do.call("rbind", data)


# Change to numeric variable
data_all$doxorubicin <- as.numeric(levels(data_all$doxorubicin))[data_all$doxorubicin]
names(data_all)[names(data_all) == "doxorubicin"] <- "A"
names(data_all)[names(data_all) == "pCR"] <- "Y"

# Divide data into training and test
set.seed(3)
ind_test <- c()
for (i in unique(data_all$S)) {
  # Get the indices of the rows corresponding to the current study
  study_indices <- which(data_all$S == i)
  
  # Sample 30% of these indices
  study_test_indices <- sample(study_indices, round(0.3 * length(study_indices)))
  
  # Add the sampled indices to the ind_test vector
  ind_test <- c(ind_test, study_test_indices)
}

# Names of covariates
covariate_cols <- c("age", "ER_preTrt", "PR_preTrt", "HER2_preTrt", cn_int_final)
covariate_cols_continuous <- c("age", cn_int_final)
p <- length(covariate_cols)

set.seed(3)
beta_tau <- rnorm(p + 1, 0, 1)


data_application <- function(data_all, K, p, alpha, beta_tau, perturb_tau, covariate_cols, covariate_cols_continuous, numRE, k_folds, ind_test){
  beta_tau_k <- tau_k <- vector("list", length = K)
  X_int <- cbind(1, data_all[, covariate_cols])
  for(k in 1:K){
    # Tau
    sigma_tau <- rep(perturb_tau, numRE)
    beta_tau_k[[k]] <- beta_tau
    beta_tau_k[[k]][sample(p + 1, numRE)] <- rnorm(numRE, 0, sigma_tau)
    tau_k[[k]] <- as.matrix(X_int[data_all$S == k, ]) %*% as.matrix(beta_tau_k[[k]])
    tau_k[[k]] <- 2 * expit(tau_k[[k]]/50) - 1
  }
  # Sample R based on |tau(x)|, if R = 0, then Y(0) = Y(1) = Y. If R = 1, then {Y(0), Y(1)} = {0, 1} if tau > 0, and vice-versa. 
  # Then let Y = Y(A)
  data_all$true_tau <- do.call("rbind", tau_k)
  data_all$R <- rbern(nrow(data_all), abs(data_all$true_tau))
  data_all$Y1 <- data_all$Y0 <- data_all$Ystar <- NA
  data_all$Y1[data_all$R == 0] <- data_all$Y0[data_all$R == 0] <- data_all$Ystar[data_all$R == 0] <- data_all$Y[data_all$R == 0]
  data_all$Y1[data_all$R == 1 & data_all$true_tau > 0] <- 1
  data_all$Y0[data_all$R == 1 & data_all$true_tau > 0] <- 0
  data_all$Y1[data_all$R == 1 & data_all$true_tau < 0] <- 0
  data_all$Y0[data_all$R == 1 & data_all$true_tau < 0] <- 1
  data_all$Ystar[data_all$R == 1 & data_all$A == 1] <- data_all$Y1[data_all$R == 1 & data_all$A == 1]
  data_all$Ystar[data_all$R == 1 & data_all$A == 0] <- data_all$Y0[data_all$R == 1 & data_all$A == 0]
  data_all$Y <- data_all$Ystar
  
  # Define train and test sets
  data_test <- data_all[ind_test, ]
  data_train <- data_all[-ind_test, ]
  test_true_tau <- data_test$true_tau
  
  ####################################################################################
  #               Multi-Study R-Learner Nuisance Study-Specific                      # 
  ####################################################################################
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
  
  # Estimate m 
  m_estimate_study_specific <- rlearner_m_hat_stack_categorical(data = data_train, K = K, p = p, covariate_cols = covariate_cols, covariate_cols_continuous = covariate_cols_continuous, k_folds = k_folds, alpha = alpha, p_mat = p_pred_mat_train)
  m_estimate_merge <- rlearner_m_hat_merge(data = data_all[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, k_folds = k_folds, alpha = alpha)
  
  # Propensity score model e_k(x) for study 2 (GSE25065) and study 3 (GSE20194)
  ind2 <- which(data_train$S == 2)
  a_hat_k2 <- glmnet::cv.glmnet(x = data.matrix(data_train[ind2, covariate_cols]),
                                y = as.matrix(data_train$A[ind2]),
                                alpha = alpha,
                                lambda = NULL,
                                standardize = TRUE,
                                keep = TRUE,
                                family = "binomial",
                                type.measure = "deviance")
  a_hat_k2_lambda_min = a_hat_k2$lambda[which.min(a_hat_k2$cvm[!is.na(colSums(a_hat_k2$fit.preval))])]
  
  ind3 <- which(data_train$S == 3)
  a_hat_k3 <- glmnet::cv.glmnet(x = data.matrix(data_train[ind3, covariate_cols]),
                                y = as.matrix(data_train$A[ind3]),
                                alpha = alpha,
                                lambda = NULL,
                                standardize = TRUE,
                                keep = TRUE,
                                family = "binomial",
                                type.measure = "deviance")
  a_hat_k3_lambda_min = a_hat_k3$lambda[which.min(a_hat_k3$cvm[!is.na(colSums(a_hat_k3$fit.preval))])]
  
  
  X_int_train <- cbind(1, data_train[, covariate_cols])
  
  X_res_mat_avg1 <- sweep(X_int_train, 1, (data_train$A - 0.5) * p_pred_mat_train[, 1], "*")
  X_res_mat_avg2 <- sweep(X_int_train, 1, (data_train$A - predict(a_hat_k2, newx = data.matrix(data_train[, covariate_cols]), s = a_hat_k2_lambda_min, type = "response")) * p_pred_mat_train[, 2], "*")
  X_res_mat_avg3 <- sweep(X_int_train, 1, (data_train$A - predict(a_hat_k3, newx = data.matrix(data_train[, covariate_cols]), s = a_hat_k3_lambda_min, type = "response")) * p_pred_mat_train[, 3], "*")
  
  X_res_mat_avg <- list(X_res_mat_avg1, X_res_mat_avg2, X_res_mat_avg3)
  
  X_res_mat_avg_all <- as.data.frame(do.call("cbind", X_res_mat_avg))
  names(X_res_mat_avg_all) <- sapply(1:ncol(X_res_mat_avg_all), function(x) paste("X.", x, sep = ""))
  
  Y_res_ss_estimate <- data_all[-ind_test, ]$Y - m_estimate_merge
  X_int_test <- cbind(1, as.matrix(data_test[, covariate_cols]))
  ms_rlearner_mSS_aSS_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_ss_estimate, alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                          newX = X_int_test)
  # ms_rlearner_pred <- 2 * expit(apply(ms_rlearner_pred_mat, 1, sum)) - 1
  ms_rlearner_pred <- 2 * expit(ms_rlearner_mSS_aSS_pred) - 1
  
  
  
  #############################################################################
  #                      R-learner Estimate Study-Specific                    #
  #############################################################################
  rlearner_estimate_SS <- rlearner_estimate_SS_pred <- vector("list", length = K)
  for(k in 1:K){
    data_k <- data_train[which(data_train$S == k), ]
    if(k == 1){
      # First study is an RCT
      p_hat = 0.5
    }else if(k == 2){
      # Second study is an OS
      p_hat = predict(a_hat_k2, newx = as.matrix(data_k[, covariate_cols]), s = a_hat_k2_lambda_min, type = "response")
    }else if(k == 3){
      # Thirs study is an OS
      p_hat = predict(a_hat_k3, newx = as.matrix(data_k[, covariate_cols]), s = a_hat_k3_lambda_min, type = "response")
    }
    rlearner_estimate_SS[[k]] <- rlasso(x = as.matrix(data_k[, covariate_cols]), 
                                        w = data_k$A, 
                                        y = data_k$Y,
                                        lambda_choice = "lambda.min",
                                        p_hat = p_hat,
                                        m_hat = m_estimate_study_specific$ss_by_study[[k]][, 1],
                                        alpha = alpha)
    rlearner_estimate_SS_pred[[k]] <- predict(rlearner_estimate_SS[[k]], newx = as.matrix(data_test[, covariate_cols]))
  }
  
  rlearner_pred <- apply(do.call("cbind", rlearner_estimate_SS_pred) * p_pred_mat_test, 1, sum)
  rlearner_pred <-  2 * expit(rlearner_pred) - 1
  
  ####################################################################################
  #                                 R-Learner Merged                                 # 
  ####################################################################################
  # Estimate a 
  a_estimate_merge <- rlearner_a_hat_merge_categorical(data = data_all[-ind_test, ], K = K, p = p, covariate_cols = covariate_cols, covariate_cols_continuous = covariate_cols_continuous, k_folds = k_folds, alpha = alpha)
  
  rlearner_merge <- rlasso(x = as.matrix(data_all[-ind_test, covariate_cols]), 
                           w = data_all[-ind_test, ]$A, 
                           y = data_all[-ind_test, ]$Y,
                           lambda_choice = "lambda.min",
                           p_hat = a_estimate_merge,
                           m_hat = m_estimate_merge,
                           alpha = alpha)
  
  rlearner_merge_pred <- predict(rlearner_merge, newx = as.matrix(data_all[ind_test, covariate_cols]))
  rlearner_merge_pred <- 2 * expit(rlearner_merge_pred) - 1
  
  out <- data.frame(rlearner =  mean((test_true_tau - rlearner_pred)^2),
                    rlearner_merge = mean((test_true_tau - rlearner_merge_pred)^2), 
                    ms_rlearner = mean((test_true_tau - ms_rlearner_pred)^2),
                    test_true_tau = test_true_tau,
                    rlearner_pred = rlearner_pred, 
                    rlearner_merge_pred = rlearner_merge_pred, 
                    ms_rlearner_pred = ms_rlearner_pred,
                    perturb_tau = perturb_tau)
  return(out)
}

data_application_multiple <- function(data_all, K, p, alpha, beta_tau, perturb_tau, covariate_cols, covariate_cols_continuous, numRE, ind_test, nreps){
  registerDoMC(cores = 48)
  
  results = foreach(j = 1:nreps, .combine = rbind) %dopar% {
    print(paste("Perturb_tau =", perturb_tau, "Iteration =", j))
    
    # Using tryCatch to handle errors
    result <- tryCatch({
      data_application(data_all = data_all, K = K, p = p, alpha = alpha, beta_tau = beta_tau, 
                       perturb_tau = perturb_tau,
                       covariate_cols = covariate_cols, covariate_cols_continuous = covariate_cols_continuous, numRE = numRE, k_folds = k_folds, ind_test = ind_test)
    }, error = function(e) {
      cat(paste("Error in iteration", j, ":", e$message, "\n"))
      return(NULL)  # Return NULL if there is an error
    })
    return(result)
  }
  results <- results[!sapply(results, is.null),]
  return(results)
}

nreps <- 500
perturb_tau <- seq(0, 5, 1)
results <- vector("list", length = length(perturb_tau))
ind_t <- 1:length(perturb_tau)
numRE <- 20
k_folds <- 3
for(t in ind_t){
  results[[t]] = as.data.frame(data_application_multiple(data_all = data_all, K = K, p = p, alpha = alpha, beta_tau = beta_tau, 
                                                         perturb_tau = perturb_tau[t], 
                                                         covariate_cols = covariate_cols, covariate_cols_continuous = covariate_cols_continuous,
                                                         ind_test = ind_test,
                                                         numRE = numRE, nreps = nreps))
}


save(results, perturb_tau, data_all, K, p, beta_tau, numRE, file = paste0("breastDataApplication_p100_numRE20_Final_nrep500", ".RData"))

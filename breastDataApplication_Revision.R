source("Utils.R")
library(curatedBreastData)
data(clinicalData)
data(curatedBreastDataExprSetList)
load("proc_curatedBreastDataExprSetList.RData")


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

# alpha for glmnet
alpha <- 0
k_folds <- 3

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
    data[[i]] <- cbind(data[[i]], t(esets1[cn_int_final, as.character(tmp$patient_ID)]))
  }else if(i == 2){
    data[[i]] <- cbind(data[[i]], t(esets2[cn_int_final, as.character(tmp$patient_ID)]))
  }else if(i == 3){
    data[[i]] <- cbind(data[[i]], t(esets3[cn_int_final, as.character(tmp$patient_ID)]))
  }
  data[[i]] <- data[[i]][complete.cases(data[[i]]), ]
  data[[i]]$S <- i
}
data_all <- do.call("rbind", data)


# Change to numeric variable
data_all$doxorubicin <- as.numeric(levels(data_all$doxorubicin))[data_all$doxorubicin]
names(data_all)[names(data_all) == "doxorubicin"] <- "A"
names(data_all)[names(data_all) == "pCR"] <- "Y"

# Names of covariates
covariate_cols <- c("age", "ER_preTrt", "PR_preTrt", "HER2_preTrt", cn_int_final)
covariate_cols_continuous <- c("age", cn_int_final)
p <- length(covariate_cols)

test_study <- 3

if(test_study == 3){
  ind_train <- which(data_all$S < 3)
  ind_test <- which(data_all$S == 3)
}else if(test_study == 2){
  ind_train <- which(data_all$S != 2)
  ind_test <- which(data_all$S == 2)
}

train_study <- setdiff(1:3, test_study)

X_int <- cbind(1, data_all[, covariate_cols])

# Define train and test sets
data_test <- data_all[ind_test, ]
data_train <- data_all[-ind_test, ]

K <- length(unique(data_train$S))

####################################################################################
#               Multi-Study R-Learner Nuisance Study-Specific                      # 
####################################################################################
# Estimate ascertainment probability
set.seed(3)
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
if(test_study == 3){
  ind2 <- which(data_train$S == 2)
  set.seed(3)
  a_hat_k2 <- glmnet::cv.glmnet(x = data.matrix(data_train[ind2, covariate_cols]),
                              y = as.matrix(data_train$A[ind2]),
                              alpha = alpha,
                              lambda = NULL,
                              standardize = TRUE,
                              keep = TRUE,
                              family = "binomial",
                              type.measure = "deviance")
  a_hat_k2_lambda_min = a_hat_k2$lambda[which.min(a_hat_k2$cvm[!is.na(colSums(a_hat_k2$fit.preval))])]
}else if(test_study == 2){
  ind3 <- which(data_train$S == 3)
  set.seed(3)
  a_hat_k3 <- glmnet::cv.glmnet(x = data.matrix(data_train[ind3, covariate_cols]),
                              y = as.matrix(data_train$A[ind3]),
                              alpha = alpha,
                              lambda = NULL,
                              standardize = TRUE,
                              keep = TRUE,
                              family = "binomial",
                              type.measure = "deviance")
  a_hat_k3_lambda_min = a_hat_k3$lambda[which.min(a_hat_k3$cvm[!is.na(colSums(a_hat_k3$fit.preval))])]
}

X_int_train <- cbind(1, data_train[, covariate_cols])

X_res_mat_avg1 <- sweep(X_int_train, 1, (data_train$A - 0.5) * p_pred_mat_train[, 1], "*")
if(test_study == 3){
  X_res_mat_avg2 <- sweep(X_int_train, 1, (data_train$A - predict(a_hat_k2, newx = data.matrix(data_train[, covariate_cols]), s = a_hat_k2_lambda_min, type = "response")) * p_pred_mat_train[, 2], "*")
  X_res_mat_avg <- list(X_res_mat_avg1, X_res_mat_avg2)
}else if(test_study == 2){
  X_res_mat_avg3 <- sweep(X_int_train, 1, (data_train$A - predict(a_hat_k3, newx = data.matrix(data_train[, covariate_cols]), s = a_hat_k3_lambda_min, type = "response")) * p_pred_mat_train[, 3], "*")
  X_res_mat_avg <- list(X_res_mat_avg1, X_res_mat_avg3)
}

X_res_mat_avg_all <- as.data.frame(do.call("cbind", X_res_mat_avg))
names(X_res_mat_avg_all) <- sapply(1:ncol(X_res_mat_avg_all), function(x) paste("X.", x, sep = ""))

Y_res_ss_estimate <- data_all[-ind_test, ]$Y - m_estimate_merge
X_int_test <- cbind(1, as.matrix(data_test[, covariate_cols]))
ms_rlearner_pred <- ms_rlearner(X = X_res_mat_avg_all, Y = Y_res_ss_estimate, alpha = alpha, K = K, p = p, p_mat = p_pred_mat_test,
                                        newX = X_int_test)
ms_rlearner_pred <- 2 * expit(ms_rlearner_pred) - 1

#############################################################################
#                      R-learner Estimate Study-Specific                    #
#############################################################################
rlearner_estimate_SS <- rlearner_estimate_SS_pred <- vector("list", length = K)
for(k in train_study){
  data_k <- data_train[which(data_train$S == k), ]
  if(k == 1){
    # First study is an RCT
    p_hat = 0.5
  }else if(k == 2){
    # Second study is an OS
    p_hat = predict(a_hat_k2, newx = as.matrix(data_k[, covariate_cols]), s = a_hat_k2_lambda_min, type = "response")
  }else if(k == 3){
    # Third study is an OS
    p_hat = predict(a_hat_k3, newx = as.matrix(data_k[, covariate_cols]), s = a_hat_k3_lambda_min, type = "response")
  }
  rlearner_estimate_SS[[k]] <- rlasso(x = as.matrix(data_k[, covariate_cols]), 
                                      w = data_k$A, 
                                      y = data_k$Y,
                                      lambda_choice = "lambda.min",
                                      alpha = alpha)
  rlearner_estimate_SS_pred[[k]] <- predict(rlearner_estimate_SS[[k]], newx = as.matrix(data_test[, covariate_cols]))
}

rlearner_pred <- apply(do.call("cbind", rlearner_estimate_SS_pred) * p_pred_mat_test, 1, sum)
rlearner_pred <-  2 * expit(rlearner_pred) - 1

#############################################################################
#                      Visualization                                        #
#############################################################################
summary_stats <- data.frame(
  Method = c("Multi-study R-learner", "Study-specific R-learner"),
  Mean = c(mean(ms_rlearner_pred), mean(rlearner_pred)),
  SD = c(sd(ms_rlearner_pred), sd(rlearner_pred)),
  Min = c(min(ms_rlearner_pred), min(rlearner_pred)),
  Max = c(max(ms_rlearner_pred), max(rlearner_pred))
)

print(summary_stats)

df_preds <- data.frame(
  Study = as.factor(data_test$S),
  MS_Rlearner = ms_rlearner_pred,
  SS_Rlearner = rlearner_pred
)

# Default
df_preds$Subtype <- "Other"

# Triple Negative: ER-, PR-, HER2-
df_preds$Subtype[data_test$ER_preTrt == 0 & data_test$PR_preTrt == 0 & data_test$HER2_preTrt == 0] <- "Triple Negative"

# HER2-enriched: HER2+ and ER-, PR-
df_preds$Subtype[data_test$HER2_preTrt == 1 & data_test$ER_preTrt == 0 & data_test$PR_preTrt == 0] <- "HER2+"

# Luminal B: HER2+ and ER+ or PR+
df_preds$Subtype[data_test$HER2_preTrt == 1 & (data_test$ER_preTrt == 1 | data_test$PR_preTrt == 1)] <- "Luminal B"

# Luminal A: HER2- and ER+ or PR+
df_preds$Subtype[data_test$HER2_preTrt == 0 & (data_test$ER_preTrt == 1 | data_test$PR_preTrt == 1)] <- "Luminal A"

df_preds$Subtype <- factor(df_preds$Subtype, levels = c("Luminal A", "Luminal B", "HER2+", "Triple Negative", "Other"))


g <- ggplot(df_preds, aes(x = SS_Rlearner, y = MS_Rlearner, color = Subtype)) +
  geom_point(alpha = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = expression(paste("Study-specific R-learner ", hat(tau)(x))), 
       y = expression(paste("Multi-study R-learner ", hat(tau)(x))),
       color = "Subtype", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 15),   # Adjust x-axis tick mark size
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),  # Legend title size
        legend.text = element_text(size = 13))
g

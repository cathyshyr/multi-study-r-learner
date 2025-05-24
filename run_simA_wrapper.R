args <- commandArgs(trailingOnly = TRUE)
array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

get_arg <- function(arg_name) {
  val <- sub(paste0("--", arg_name, "="), "", grep(paste0("--", arg_name, "="), args, value = TRUE))
  if (length(val) == 0) stop(paste("Missing argument:", arg_name))
  return(as.numeric(val))
}

rep_start <- get_arg("rep-start")
rep_end   <- get_arg("rep-end")
a_index   <- get_arg("a")
m_index   <- get_arg("m")
t_index   <- get_arg("t")

cat("Running reps", rep_start, "to", rep_end, "for a =", a_index, "m =", m_index, "t =", t_index, "\n")

source("ovarianDataSimulation_ScenarioA.R")

# Replace outer loop in sim_multiple
nreps = rep_end - rep_start + 1
rep_ids = rep_start:rep_end
sim_multiple <- function(n, p, K, edat_train, edat_test, covariate_cols, perturb_tau, perturb_mu0, perturb_a, k_folds, alpha, numRE, nreps){
  repeat {
    # Coefficients for generating mu0(x) = E[Y(0)|X=x]
    beta_mu0 <- rnorm(p + 1, mean = 0, sd = 1)
    # Coefficients for generating e(x) = E[A|X=x]
    beta_a <- rnorm(p + 1, mean = 0, sd = 1)
    # Coefficients for generating tau(x) = E[Y(1)-Y(0)|X=x]
    beta_tau <- rnorm(p + 1, mean = 0, sd = 1)
    # Coefficients for generating p(k|x) = P(S=k|X=x)
    beta_p_k <- vector("list", length = K - 2)
    for(k in 1:(K - 2)){
      beta_p_k[[k]] <- rnorm(p + 1, mean = 0, sd = 1)   
      # beta_p_k[[k]][sample(2:(p + 1), round(p * (4/5)))] <- 0
    }
    
    # p := P(S|X = x) ascertainment probability
    X_int <- cbind(1, rbind(edat_train, edat_test))
    denom_mat <- sapply(1:length(beta_p_k), function(x){
      exp(as.matrix(X_int) %*% beta_p_k[[x]])
    })
    denom <- apply(denom_mat, 1, sum)
    
    p_vec_dat <- data.frame(sapply(1:length(beta_p_k), function(x){
      exp(as.matrix(X_int) %*% beta_p_k[[x]])/(1 + denom)
    }))
    names(p_vec_dat) <- sapply(2:(K - 1), function(x) paste("p", x, sep = ""))
    p1 <- 1 - apply(p_vec_dat, 1, sum)
    data_train <- data.frame(cbind(edat_train, p1 = p1[1:nrow(edat_train)], p_vec_dat[1:nrow(edat_train), ]))
    data_train <- data_train %>% 
      rowwise %>%
      dplyr::mutate(S = sample(1:(K - 1), size = 1, replace = TRUE, 
                               prob = c_across(matches('^p\\d+$')))) %>% 
      ungroup
    data_train <- as.data.frame(data_train[, c(covariate_cols, "S")])
    p_mat <- cbind(p1, p_vec_dat)
    
    # Indices of random effects
    cols_re_tau <- cols_re_mu0 <-cols_re_a <- vector("list", length = K)
    for(k in 1:(K)){
      cols_re_tau[[k]] <- sample(1:(p + 1), numRE)
      
      cols_re_mu0[[k]] <- sample(1:(p + 1), numRE)
      
      cols_re_a[[k]] <- sample(1:(p + 1), numRE)
    }
    
    data_test <- edat_test
    data_test$S <- K
    data <- rbind(data_train, data_test)
    data$A <- NA
    
    # Check if every level of S has more than 100 observations
    tab <- table(data$S)
    if (all(tab > 100)) break
  }
  
  registerDoMC(cores = 48)
  results = foreach (j = rep_start:rep_end, .combine = rbind) %dopar% {
    print(paste("Perturb_a =", perturb_a, "Perturb_mu0 =", perturb_mu0, "Perturb_tau =", perturb_tau, "Iteration =", j))
    tryCatch({
      sim_each(n = n, p = p, K = K, data = data, covariate_cols = covariate_cols, 
               cols_re_tau = cols_re_tau, cols_re_mu0 = cols_re_mu0, cols_re_a = cols_re_a,
               beta_tau = beta_tau, beta_mu0 = beta_mu0, beta_a = beta_a, 
               p_mat = p_mat,
               perturb_tau = perturb_tau, perturb_mu0 = perturb_mu0, perturb_a = perturb_a,
               k_folds = k_folds, alpha = alpha, iteration = j)
    }, error = function(e) {
      message(paste("Error in iteration", j, ":", conditionMessage(e)))
      return(NA)  # or return(NA), or a row of NAs to keep dimensions consistent
    })
  }
}

perturb_a <- c(0, 1, 2)
perturb_mu0 <- c(0, 1, 2)
perturb_tau <- c(0, 0.5, 1)

# Select based on index (assume passed in as 1-based)
perturb_a_val   <- perturb_a[a_index]
perturb_mu0_val <- perturb_mu0[m_index]
perturb_tau_val <- perturb_tau[t_index]

k_folds <- 3
alpha <- 0.9

# Number of random effects
numRE <- 10

results <- sim_multiple(
  n = n, p = p, K = K,
  edat_train = edat_train, edat_test = edat_test,
  covariate_cols = covariate_cols,
  perturb_tau = perturb_tau_val,
  perturb_mu0 = perturb_mu0_val,
  perturb_a = perturb_a_val,
  k_folds = k_folds, alpha = alpha,
  numRE = numRE,
  nreps = rep_end - rep_start + 1
)

save(results, file = paste0(
  "Simulation_A",
  "_a", a_index,
  "_m", m_index,
  "_t", t_index,
  "_task", array_task_id,
  ".RData"
))

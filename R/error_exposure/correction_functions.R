# ------------------------------------------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS - EXPOSURE
# ------------------------------------------------------------------------------
# This file contains functions for implementing the following ME correction
# methods: 1) Conditional Scoring, 2) Regression Calibration,
# 3) SIMEX, 4) IV, 5) MIME (Webb-Vargas, 2015)

# --------------------------------------------------------
#            NAIVE AND IDEAL APPROACH
# --------------------------------------------------------

srf_ideal <- function(data) {
  
  # Use erf() function to estimate contrast
  results <- with(data, srf(y = Y, a = A, x = data.frame(W, X), delta = 1))
  return(results)
  
}

srf_naive <- function(data, sl.lib, family = gaussian()) {
  
  #' Naively estimates ERF contrast, using mismeasured W in place of true measurements
  
  # Use srf() function to est ATE
  results <- with(data, srf(y = Y, a = A.star, x = data.frame(W, X), delta = 1))
  return(results)
  
}

# ---------------------------------------------------------
#                REGRESSION CALIBRATION
# ---------------------------------------------------------

srf_rc <- function(data, nboot = 100) {
  
  #' Regression calibration
  #' 
  #' INPUTS: 
  #' - data from gen_data()
  #' 
  #' OUPUTS:
  #' - A list containing rcal EST and CI
  
  val_data <- data[which(data$val.idx == 1),]
  
  # fit calibration model
  A.hat <- predict(lm(A ~ A.star + W + X, data = val_data), newdata = data)
  
  # fit model with calibrated exposures
  rc_mod <- with(data, srf(a = A.hat, y = Y, x = data.frame(W, X), delta = 1))
  
  Mhat <- rc_mod$EST
  boot <- rep(NA, nboot)
  
  for (b in 1:nboot) {
    
    boot_idx <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    boot_data <- data[boot_idx,]
    val_data <- data[which(boot_data$val.idx == 1),]
    boot_data$A.tilde <- predict(lm(A ~ A.star + W + X, data = val_data), newdata = boot_data)
    
    boot_mod <- with(boot_data, srf(a = A.tilde, y = Y, x = data.frame(W, X), delta = 1))
    boot[b] <- boot_mod$EST
    
  }
  
  # bootstrap to obtain confidence interval
  CI_hat <- quantile(boot, probs = c(0.025, 0.975))
  
  return(list(EST = Mhat, CI = CI_hat))
  
}

# ---------------------------------------------------------
#                          SIMEX
# ---------------------------------------------------------

srf_simex <- function(data, nboot = 50, k = 3, lambda = seq(0, 2, by = 0.2)) {
  
  sig_u_hat <- with(data, sd(A.star[which(val.idx == 1)] - A[which(val.idx == 1)]))
  
  l.vals <- lapply(lambda, function(lam, data, sig_u, nboot){
    
    A.mat <- replicate(nboot, data$A.star + sqrt(lam)*rnorm(nrow(data), 0, sig_u))
    
    vals <- apply(A.mat, 2, function(A.tmp, Y, W, Z) {
      
      # Estimate srf
      results <- srf(a = A.tmp, y = Y, x = data.frame(W, Z), delta = 1)
      Mhat <- results$EST
      Vhat <- results$VAR
      
      return(c(EST = Mhat, VAR = Vhat))
      
    }, W = data$W, Z = data$X, Y = data$Y)
    
    mu.vals <- vals[1,] 
    sig.vals <- vals[2,]
    
    s.hat <- var(mu.vals)
    sig.hat <- mean(sig.vals)
    
    return(list(est = mean(mu.vals), var = sig.hat - s.hat))
    
  }, data = data, nboot = nboot, sig_u = sig_u_hat)
  
  Psi <- do.call(c, lapply(l.vals, function(o) o$est))
  Phi <- do.call(c, lapply(l.vals, function(o) o$var))
  
  Mhat <- c(predict(mgcv::gam(Psi ~ s(lambda, k = k), data = data.frame(Psi = Psi, lambda = lambda)), 
                  newdata = data.frame(lambda = -1)))
  Vhat <- c(predict(mgcv::gam(Phi ~ s(lambda, k = k), data = data.frame(Phi = Phi, lambda = lambda)), 
                  newdata = data.frame(lambda = -1)))
  
  CI_hat <- c(qnorm(.025)*sqrt(Vhat) + Mhat, qnorm(.975)*sqrt(Vhat) + Mhat)
  
  return(list(EST = Mhat, CI = CI_hat))
  
}

# --------------------------------------------------------
#                 INSTRUMENTAL VARIABLES
# --------------------------------------------------------

srf_iv <- function(data) {
  
  #' Implements instrumental variables correction 
  #'
  #' INPUTS: 
  #' - data: A dataframe created with the generate_data() function
  #' 
  #' OUTPUTS:
  #' - A list containing the estimate and confidence interval
  #' 
  #' REQUIRED PACKAGES:
  #' - AER
  
  # Fit the IV model (first and second stage) with AER package
  iv_mod <- ivreg(Y ~ X + W + A.star | V + X + W, data=data)
  
  return(list(EST = iv_mod$coefficients['A.star'], # point estimate
              CI = confint(iv_mod)['A.star',])) # confidence interval
  
}

# --------------------------------------------------------
#               MULTIPLE IMPUTATION
# --------------------------------------------------------

srf_mime <- function(data, m = 50) {
  
  #' Performs multiple imputation for ME correction, following Josey et 
  #' al. (2023) but with a validation set
  
  # Turn into missing data problem explicitly by recasting the unobserved X
  # values as missing so that we can use mice
  data_imp <- data %>% mutate(A = replace(A, val.idx == 0, NA)) %>% dplyr::select(-val.idx)
  
  # run the (congenial) MI procedure
  invisible(capture.output(imps <- mice(data_imp, method = 'norm.boot', m = m)))
  
  # loop through imputed datasets, estimate PS 
  vals <- sapply(1:m, function(d, imps) {
    
    # fit propensity score model with d-th dataset
    curr_data <- complete(imps, d)
    
    # Record ERF estimate and its SE
    results <- with(curr_data, srf(y = Y, a = A, x = data.frame(W, X), delta = 1))
    est <- results$EST
    var <- results$VAR
    
    return(c(est = est, var = var))
    
  }, imps = imps)
  
  # grab ests and their SEs
  ests <- vals[1,]
  vars <- vals[2,]
  
  # compute SE est 
  Mhat <- mean(ests) 
  bw_var <- sum((ests - Mhat)^2)/(length(ests) - 1) # var bw ests
  wi_var <- mean(vars) # within variance
  Vhat <- wi_var + (1 + (1/length(ests))) * bw_var
  CI_hat <- c(qnorm(.025)*sqrt(Vhat) + Mhat, qnorm(.975)*sqrt(Vhat) + Mhat)
  
  return(list(EST = Mhat, CI = CI_hat))
  
}

# --------------------------------------------------------
#                 CONTROL VARIATES
# --------------------------------------------------------

srf_cv <- function(data) {
  
  #' Given a dataframe data (containing a validation data indicator), implements
  #' the control variates method to obtain ERF estimate
  #'
  #' INPUTS:
  #' - a dataframe from the gen_data() function
  #' OUTPUTS:
  #' - a list containing the ERF estimate, variance estimate and 95% CI
  
  demean <- function(vec) {
    return(vec-mean(vec, na.rm = T))
  }
  
  # Keep track of val data
  data_val <- data %>% filter(val.idx == 1)
  
  # Step 1: estimate ERF in validation data
  psi <- srf(a = data_val$A, y = data_val$Y, delta = 1,
                     x = subset(data_val, select = c("W","X"))) 
  
  tau_hat_val <- psi$EST
  v_hat <- psi$VAR
  varphi_val <- psi$EIF
  
  # Step 2: estimate control variates
  phi1 <- srf(a = data$A.star, y = data$Y, delta = 1,
              x = subset(data, select = c("W","X")))
  
  phi2 <- srf(a = data_val$A.star, y = data_val$Y, delta = 1,
              x = subset(data_val, select = c("W","X")))
  
  tau_ep_main <- phi1$EST
  tau_ep_val <- phi2$EST
  
  # Step 3: estimate Gamma and V
  phi_main <- phi1$EIF
  phi_val <- phi2$EIF
  
  # Get sample sizes
  n_main <- length(phi_main)
  n_val <- length(phi_val)
  
  # Estimate Gamma
  IF <- cbind(demean(phi_val), demean(varphi_val))
  IF <- IF[complete.cases(IF),]
  gamma_hat <- (1 - n_val/n_main)*cov(IF)[1,2]
  
  # Estimate V
  V_hat <- (1 - n_val/n_main)*mean(demean(phi_main)^2, na.rm = TRUE)
  
  ## Step 4: subtract off (may need to flip the subtraction sign)
  tau_cv <- unname(tau_hat_val - (gamma_hat/V_hat)*(tau_ep_val - tau_ep_main))
  
  # Get variance estimate and 95% CI
  var_hat <- (v_hat*n_val - gamma_hat^2/V_hat)/n_val
  ci_low <- tau_cv - qnorm(0.975)*sqrt(var_hat)
  ci_high <- tau_cv + qnorm(0.975)*sqrt(var_hat)
  
  # Return the ERF est and associated variance
  return(list(EST = tau_cv, CI = c(ci_low,ci_high)))
  
}
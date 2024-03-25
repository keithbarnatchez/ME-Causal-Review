# ------------------------------------------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS - EXPOSURE
# ------------------------------------------------------------------------------
# This file contains functions for implementing the following ME correction
# methods: 1) Conditional Scoring, 2) Regression Calibration,
# 3) SIMEX, 4) IV, 5) MIME (Webb-Vargas, 2015)

# --------------------------------------------------------
#            NAIVE AND IDEAL APPROACH
# --------------------------------------------------------

erf_ideal <- function(data) {
  
  # Use erf() function to estimate contrast
  results <- with(data, erf(y = Y, a = A, x = data.frame(W, X), a1 = 1, a0 = 0))
  return(results)
  
}

erf_naive <- function(data, sl.lib, family = gaussian()) {
  
  #' Naively estimates ERF contrast, using mismeasured W in place of true measurements
  
  # Use erf() function to est ATE
  results <- with(data, erf(y = Y, a = A.star, x = data.frame(W, X), a1 = 1, a0 = 0,))
  return(results)
  
}

# ---------------------------------------------------------
#                CONDITIONAL SCORE METHOD
# ---------------------------------------------------------

csme_implement <- function(data) {
  
  Y <- data$Y
  A.star <- data$A.star
  ipw_hat <- data$ipw_hat
  X <- as.matrix(data %>% select(-c(A.star, Y, ipw_hat)))
  
  #data dimensions
  n <- nrow(X)
  p <- ncol(X)
  
  # not confident about the denominator to this estimator
  condexp <- function(X, alpha, beta, delta, sig2_u, sig2_e) {
    (c(X %*% alpha) + delta*c(X %*% beta)) /
      (1 + c(X %*% beta)^2*c(sig2_u / sig2_e))
  }
  
  condvar <- function(X, beta, sig2_u, sig2_e) {
    sig2_e / (1 + c(X %*% beta)^2*c(sig2_u / sig2_e))
  }
  
  function(theta, sig2_u, a0, a1){
    
    # model parameters
    gamma <- theta[1:p]
    alpha <- theta[(p + 1):(2*p)]
    beta <- theta[(2*p + 1):(3*p)]
    sig2_e <- theta[(3*p + 1)]
    eta <- theta[(3*p + 2)]
    
    # gps model uncertainty (questionable whether this is sufficient?)
    gps_eqn <- crossprod(X, A.star - c(X %*% gamma))
    
    # conditional score statistics
    delta <- A.star + c(sig2_u/sig2_e)*Y*c(X %*% beta)
    m_A <- condexp(X = X, alpha = alpha, beta = beta, 
                   delta = delta, sig2_u = sig2_u, sig2_e = sig2_e)
    v_A <- condvar(X = X, beta = beta, sig2_u = sig2_u, sig2_e = sig2_e)
    scl <- (Y - m_A)^2 / v_A
    
    # outcome model estimating equations
    ols_eqn1 <- crossprod(X, ipw_hat*(Y - m_A))
    ols_eqn2 <- crossprod(delta*X, ipw_hat*(Y - m_A))
    disp_eqn <- crossprod(ipw_hat, sig2_e - sig2_e*scl)
    
    # prediction estimating equation
    pred_eqn <- t(rep(a1 - a0, n)) %*% c(X %*% beta) - eta
    
    c(gps_eqn, ols_eqn1, ols_eqn2, disp_eqn, pred_eqn)   
    
  }
  
}

erf_csme <- function(data) {
  
  #' Function for implementing a modified version of the Blette (2022) CSME
  #' approach
  #'
  #' INPUTS:
  #' - data: Simulation vars from gen_data()
  sig_u_hat <- with(data, sd(A.star[which(val.idx == 1)] - A[which(val.idx == 1)]))
  
  # fit GPS model/get weights
  denom_mod <- lm(A.star ~ W + X, data = data)
  p_denom <- predict(denom_mod, type = 'response')
  dens_denom <- dnorm(data$A.star, p_denom, sigma(denom_mod))
  num_mod <- lm(A.star ~ 1, data = data)
  p_num <- predict(num_mod, type = 'response')
  dens_num <- dnorm(data$A.star, p_num, sigma(denom_mod))
  data$ipw_hat <- dens_num / dens_denom
  
  # set up data frame for m_estimate
  mdat <- with(data, data.frame("Y" = Y, "A.star" = A.star, 
                                "ipw_hat" = ipw_hat, "X" = cbind(1, W, X)))
  
  # initial predictions
  outmod <- lm(Y ~ A.star*X + A.star*W, data = data)
  ipwmod <- lm(Y ~ A.star, weights = ipw_hat, data = data)
  start <- unname(c(coef(denom_mod), coef(outmod), sigma(outmod)^2, coef(ipwmod)[2]))
  
  results_csme <- m_estimate(estFUN = csme_implement, data = mdat,
                             inner_args = list(sig2_u = (sig_u_hat)^2, a0 = 0, a1 = 1),
                             root_control = setup_root_control(start = start))
  
  idx <- length(results_csme@estimates)
  
  # get results of interest
  Mhat <- coef(results_csme)[idx]
  Vhat <- vcov(results_csme)[idx, idx]
  
  CI_hat <- c(Mhat - 1.96*sqrt(Vhat), Mhat + 1.96*sqrt(Vhat))
  
  return(list(EST = Mhat, CI = CI_hat))
  
}

# ---------------------------------------------------------
#                REGRESSION CALIBRATION
# ---------------------------------------------------------

erf_rc <- function(data, nboot = 100) {
  
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
  rc_mod <- with(data, erf(a = A.hat, y = Y, x = data.frame(W, X), a1 = 1, a0 = 0))
  
  Mhat <- rc_mod$EST
  boot <- rep(NA, nboot)
  
  for (b in 1:nboot) {
    
    boot_idx <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    boot_data <- data[boot_idx,]
    val_data <- data[which(boot_data$val.idx == 1),]
    boot_data$A.tilde <- predict(lm(A ~ A.star + W + X, data = val_data), newdata = boot_data)
    
    boot_mod <- with(boot_data, erf(a = A.tilde, y = Y, x = data.frame(W, X), a1 = 1, a0 = 0))
    boot[b] <- boot_mod$EST
    
  }
  
  # bootstrap to obtain confidence interval
  CI_hat <- quantile(boot, probs = c(0.025, 0.975))
  
  return(list(EST = Mhat, CI = CI_hat))
  
}

# ---------------------------------------------------------
#                          SIMEX
# ---------------------------------------------------------

erf_simex <- function(data, nboot = 50, k = 3, lambda = seq(0, 2, by = 0.2)) {
  
  sig_u_hat <- with(data, sd(A.star[which(val.idx == 1)] - A[which(val.idx == 1)]))
  
  l.vals <- lapply(lambda, function(lam, data, sig_u, nboot){
    
    A.mat <- replicate(nboot, data$A.star + sqrt(lam)*rnorm(nrow(data), 0, sig_u))
    
    vals <- apply(A.mat, 2, function(A.tmp, Y, W, Z) {
      
      # Estimate ERF
      results <- erf(a = A.tmp, y = Y, x = data.frame(W, Z), a0 = 0, a1 = 1)
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

erf_iv <- function(data) {
  
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

erf_mime <- function(data, m = 50) {
  
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
    results <- with(curr_data, erf(y = Y, a = A, x = data.frame(W, X), a0 = 0, a1 = 1))
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

erf_cv <- function(data) {
  
  #' Given a dataframe data (containing a validation data indicator), implements
  #' the control variates method to obtain ERF estimate
  #'
  #' INPUTS:
  #' - a dataframe from the gen_data() function
  #' OUTPUTS:
  #' - a list containing the ERF estimate, variance estimate and 95% CI
  
  demean <- function(vec) {
    return(vec-mean(vec))
  }
  
  # Keep track of val data
  data_val <- data %>% filter(val.idx == 1)
  
  # Step 1: estimate ERF in validation data
  tau_val_mod <- erf(a = data_val$A, y = data_val$Y, 
                     x = subset(data_val, select = c("W","X")),
                     a0 = 0, a1 = 1, sl.lib = c("SL.glm")) 
  
  tau_hat_val <- tau_val_mod$EST
  v_hat <- tau_val_mod$VAR
  varphi_val <- tau_val_mod$EIF
  
  # Step 2: estimate control variates
  psi1 <- erf(a = data$A.star, y = data$Y, 
              x = subset(data, select = c("W","X")),
              a0 = 0, a1 = 1, sl.lib = c("SL.glm"))
  
  psi2 <- erf(a = data_val$A.star, y = data_val$Y, 
              x = subset(data_val, select = c("W","X")),
              a0 = 0, a1 = 1, sl.lib = c("SL.glm"))
  
  cv_mods <- list(psi1 = psi1, psi2 = psi2)
  
  tau_ep_val <- psi1$EST
  tau_ep_main <- psi2$EST
  
  # Step 3: estimate Gamma and V
  phi_main <- psi1$EIF
  phi_val <- psi2$EIF
  
  # Get sample sizes
  n_main <- length(phi_main)
  n_val <- length(phi_val)
  
  # Estimate Gamma
  gamma_hat <- (1 - n_val/n_main)*(1/n_val)*cov(cbind(demean(phi_val), demean(varphi_val)))[1,2]
  
  # Estimate V
  V_hat <- (1 - n_val/n_main)*(1/n_val)*mean(demean(phi_main)^2)
  
  ## Step 4: subtract off (may need to flip the subtraction sign)
  tau_cv <- unname(tau_hat_val - (gamma_hat/V_hat)*(tau_ep_main - tau_ep_val))
  
  # Get variance estimate and 95% CI
  var_hat <- v_hat - gamma_hat^2/V_hat
  ci_low <- tau_cv - qnorm(0.975)*sqrt(var_hat)
  ci_high <- tau_cv + qnorm(0.975)*sqrt(var_hat)
  
  # Return the ERF est and associated variance
  return(list(EST = tau_cv, CI = c(ci_low,ci_high)))
  
}

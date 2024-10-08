# ------------------------------------------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS - CONFOUNDER
# ------------------------------------------------------------------------------
# This file contains functions for implementing the following ME correction
# methods: 1) Propensity Score Calibration, 2) Regression Calibration,
# 3) SIMEX, 4) IV, 5) MIME (Webb-Vargas, 2015)

# --------------------------------------------------------
#            NAIVE AND IDEAL APPROACH
# --------------------------------------------------------

srf_ideal <- function(data) {
  
  #' Computes ATE under ideal conditions (using X, correctly specified model)
  #'
  #' Returns ATE estimate
  
  # Estimate ATE
  mod <- with(data, srf(a = A, y = Y, x = data.frame(W, X), delta = 1))
  
  return(mod)
  
}

srf_naive <- function(data) {
  
  #' Computes ATE when naively using W in place of X in estimating propensity
  #' scores, but with otherwise correctly-specified model
  #' 
  #' Returns ATE estimates
  
  # Estimate ATE
  mod <- with(data, srf(a = A, y = Y, x = data.frame(W.star, X), delta = 1))
  
  return(mod)
  
}

# ---------------------------------------------------------
#                PROPENSITY SCORE CALIBRATION
# ---------------------------------------------------------

psc_implement <- function(data) {
  
  #' Fit propensity score model in main data so it can be projected to GS measure
  #' via the model fitted with validation data
  
  # Estimate error-prone proopensity score
  data$ep_ps <- predict(glm(A ~ W.star + X, data = data, family = binomial()), type = 'response') # getting ep pscores in main data
  
  # Fit propensity score models in validation data
  val_data <- data[which(data$val.idx == 1),]
  gs_ps_mod <- glm(A ~ W + X, data = val_data, family = binomial()) # gold standard
  val_data$gs_ps <- predict(gs_ps_mod, type = "response") # predict gs pscores in val data

  # Fit model relating gold-standard PS to error-prone PS
  ps_rel_model <- lm(gs_ps ~ A + ep_ps, data = val_data)

  # fit error-prone outcome model
  ps_epe_model <- lm(Y ~ A + ep_ps, data = data)
  
  return(list(ps_rel_model, ps_epe_model))
  
}

srf_psc <- function(data, nboot = 100) {
  
  #' Outer function for the propensity score calibration method. Calls
  #' psc_implement to yield ATE estimate, and psc_bootstrap to obtain con
    
  # Get weights
  psc_mods <- psc_implement(data)
  
  # Estimate ATE
  B_A <- psc_mods[[2]]$coefficients['A']
  B_X <- psc_mods[[2]]$coefficients['ep_ps']
  L_A <- psc_mods[[1]]$coefficients['A']
  L_X <- psc_mods[[1]]$coefficients['ep_ps']
  EST_hat <- B_A - L_A*B_X/L_X
  
  # Bootsrap
  boot <- rep(NA, nboot)
  
  for (b in 1:nboot) {
    
    boot_idx <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    boot_mod <- psc_implement(data[boot_idx,])
    B_A <- boot_mod[[2]]$coefficients['A']
    B_X <- boot_mod[[2]]$coefficients['ep_ps']
    L_A <- boot_mod[[1]]$coefficients['A']
    L_X <- boot_mod[[1]]$coefficients['ep_ps']
    boot[b] <- B_A - L_A*B_X/L_X
    
  }
  
  # Bootstrap to obtain confidence interval
  CI_hat <- quantile(boot, probs = c(0.025, 0.975))
  
  return(list(EST = EST_hat, CI = CI_hat))
  
}

# ---------------------------------------------------------
#                REGRESSION CALIBRATION
# ---------------------------------------------------------

srf_rc <- function(data, nboot = 100) {
  
  #' Standard Regression Calibration
  
  val_data <- data[which(data$val.idx == 1),]
  W.hat <- predict(lm(W ~ W.star + A + X, data = val_data), newdata = data)
  
  # fit error-prone model
  rc_mod <- with(data, srf(a = A, y = Y, x = data.frame(W.hat, X), delta = 1))
  
  EST_hat <- rc_mod$EST
  boot <- rep(NA, nboot)

  for (b in 1:nboot) {

    boot_idx <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
    boot_data <- data[boot_idx,]
    val_data <- data[which(boot_data$val.idx == 1),]
    W.tilde <- predict(lm(W ~ W.star + A + X, data = val_data), newdata = boot_data)
    
    boot_mod <- with(boot_data, srf(a = A, y = Y, x = data.frame(W.tilde, X), delta = 1))
    boot[b] <- boot_mod$EST

  }

  # Bootstrap to obtain confidence interval
  CI_hat <- quantile(boot, probs = c(0.025, 0.975))
  
  return(list(EST = EST_hat, CI = CI_hat))
  
}

# ---------------------------------------------------------
#                          SIMEX
# ---------------------------------------------------------

srf_simex <- function(data, nboot = 50, k = 3, lambda = seq(0, 2, by = 0.2)) {
  
  #' INPUTS:
  #' - data: simulation data
  #' - degree: degree of jackknife extrapolation
  #' - lambda: measurement error incrementation
  #' - nboot: number of bootstrap resamples
  
  sig_u_hat <- with(data, sd(W.star[which(val.idx == 1)] - W[which(val.idx == 1)]))
  
  l.vals <- lapply(lambda, function(lam, data, sig_u, nboot, method){
    
    W.mat <- replicate(nboot, data$W.star + sqrt(lam)*rnorm(nrow(data), 0, sig_u))
    
    vals <- apply(W.mat, 2, function(W.tmp, A, Z, Y, method) {
      
      # Estimate ATE
      simex_mod <- srf(a = A, y = Y, x = data.frame(W.tmp, Z), delta = 1)
      
      EST_hat <- simex_mod$EST
      Vhat <- simex_mod$VAR
      
      return(c(EST = EST_hat, VAR = Vhat))
      
    }, A = data$A, Z = data$X, Y = data$Y, method = method)
    
    mu.vals <- vals[1,] 
    sig.vals <- vals[2,]
    
    s.hat <- var(mu.vals)
    sig.hat <- mean(sig.vals)
    
    return(list(est = mean(mu.vals), var = sig.hat - s.hat))
    
  }, data = data, nboot = nboot, sig_u = sig_u_hat, method = method)
  
  Psi <- do.call(c, lapply(l.vals, function(o) o$est))
  Phi <- do.call(c, lapply(l.vals, function(o) o$var))
  
  EST_hat <- predict(mgcv::gam(Psi ~ s(lambda, k = k), data = data.frame(Psi = Psi, lambda = lambda)), 
                     newdata = data.frame(lambda = -1))
  Vhat <- predict(mgcv::gam(Phi ~ s(lambda, k = k), data = data.frame(Phi = Phi, lambda = lambda)), 
                  newdata = data.frame(lambda = -1))
  
  CI_hat <- c(qnorm(.025)*sqrt(Vhat) + EST_hat, qnorm(.975)*sqrt(Vhat) + EST_hat)
  
  return(list(EST = EST_hat, CI = CI_hat))
  
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
  #' REQUIRED PACKAGES:
  #' - AER

  # Fit the IV model (first and second stage) with AER package
  iv_mod <- ivreg(Y ~ A + W.star + X | A + V + X, data = data)
  
  return(list(EST = iv_mod$coefficients['A'], # point estimate
              CI = confint(iv_mod)['A',])) # confidence interval
  
}

# --------------------------------------------------------
#               MULTIPLE IMPUTATION
# --------------------------------------------------------

srf_mime <- function(data, m = 50) {
  
  #' Performs multiple imputation for ME correction, following Webb-Vargas et 
  #' al. (2015) with the one modification that we assume there is an internal,
  #' not external, validation sample. Makes use of the 'mice' package to perform
  #' multiple imputations
  
  # Turn into missing data problem explicitly by recasting the unobserved X
  # values as missing so that we can use mice
  data_imp <- data %>% mutate(W = replace(W, val.idx == 0, NA)) %>% dplyr::select(-val.idx)
  
  # Run the (congenial) MI procedure
  invisible(capture.output(imps <- mice(data_imp, method = 'norm.boot', m = m)))
  
  # Loop through imputed datasets, estimate PS 
  vals <- sapply(1:m, function(d, imps, method) {
    
    # fit propensity score model with d-th dataset
    curr_data <- complete(imps, d)
    
    # Record ATE estimate and its SE
    mime_mod <- with(curr_data, srf(a = A, y = Y, x = data.frame(W, X), delta = 1))
    
    est <- mime_mod$EST
    var <- mime_mod$VAR
    
    return(c(est = est, var = var))
    
  }, imps = imps, method = method)
  
  # grab ests and their SEs
  ests <- vals[1,]
  vars <- vals[2,]
  
  # compute SE est 
  EST_hat <- mean(ests) 
  bw_var <- sum((ests - EST_hat)^2)/(length(ests)-1) # var bw ests
  wi_var <- mean(vars) # within variance
  Vhat <- wi_var + (1 + (1/length(ests))) * bw_var
  CI_hat <- c(qnorm(.025)*sqrt(Vhat) + EST_hat, qnorm(.975)*sqrt(Vhat) + EST_hat)
  
  return(list(EST = EST_hat, CI = CI_hat))
  
}

# --------------------------------------------------------
#                 CONTROL VARIATES
# --------------------------------------------------------

srf_cv <- function(data) {
  
  #' Given a dataframe data (containing a validation data indicator), implements
  #' the control variates method to obtain ATE estimate
  #'
  #' INPUTS:
  #' - a dataframe from the gen_data() function
  #' OUTPUTS:
  #' - a list containing the ATE estimate, variance estimate and 95% CI
  
  demean <- function(vec) {
    return(vec-mean(vec, na.rm = T))
  }
  
  # Keep track of val data
  data_val <- data %>% filter(val.idx == 1)
  
  # Step 1: estimate ATE in validation data
  psi <- srf(a = data_val$A, y = data_val$Y, delta = 1,
              x = subset(data_val, select = c("W","X"))) 
  
  tau_hat_val <- psi$EST
  v_hat <- psi$VAR
  varphi_val <- psi$EIF
  
  # Step 2: estimate control variates
  phi1 <- srf(a = data$A, y = data$Y, delta = 1,
               x = subset(data, select = c("W.star","X")))
  
  phi2 <- srf(a = data_val$A, y = data_val$Y, delta = 1,
               x = subset(data_val, select = c("W.star","X"))) 
  
  tau_ep_main <- phi1$EST
  tau_ep_val <- phi2$EST
  
  # Step 3: estimate Gamma and V
  phi_main <- phi1$EIF
  phi_val <- phi2$EIF
  
  # Get sample sizes
  n_main <- length(phi_main)
  n_val <- length(phi_val)
  
  # Estimate Gamma
  IF <- cbind(phi_val, varphi_val)
  IF <- IF[complete.cases(IF),]
  gamma_hat <- (1 - n_val/n_main)*cov(IF)[1,2]
  
  # Estimate V
  V_hat <- (1 - n_val/n_main)*mean(demean(phi_main)^2, na.rm = TRUE)
  
  ## Step 4: subtract off (may need to flip the subtraction sign)
  tau_cv <- tau_hat_val - (gamma_hat/V_hat)*(tau_ep_val - tau_ep_main)
  
  # Get variance estimate and 95% CI
  var_hat <- (v_hat*n_val - gamma_hat^2/V_hat)/n_val
  ci_low <- tau_cv - qnorm(0.975)*sqrt(var_hat)
  ci_high <- tau_cv + qnorm(0.975)*sqrt(var_hat)
  
  # Return the ATE est and associated variance
  return(list(EST = tau_cv, CI = c(ci_low,ci_high)))
  
}
 
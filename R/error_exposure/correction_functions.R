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
    
    # simulate nboot 
    z.mat <- replicate(n.boot, z + sqrt(lam)*rnorm(length(z), 0, sqrt(tau2)))
    
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
    
  }, z = z, y = y, x = x, a0 = a0, a1 = a1, sl.lib = sl.lib)
  
  if (any(lambda == 0)){
    
    Psi <- do.call(c, lapply(l.vals, function(o) o$estimate))[,-which(lambda == 0)]
    Phi <- do.call(c, lapply(l.vals, function(o) o$variance))[,-which(lambda == 0)]
    
  } else {
    
    Psi <- do.call(c, lapply(l.vals, function(o) o$estimate))
    Phi <- do.call(c, lapply(l.vals, function(o) o$variance))
    
  }
  
  # Extrapolation step
  L <- cbind(1, poly(lambda, degree = degree, raw = TRUE))
  chi <- c(1, poly(-1, degree = degree, raw = TRUE))
  estimate <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Psi)
  variance <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Phi)
  se_simex <- sqrt(variance)

  # Compile results into list and return the list
  out <- list(estimate = estimate, 
              CI = c(estimate - 1.96*se_simex,
                     estimate + 1.96*se_simex),
              variance = variance, 
              Psi = Psi, Phi = Phi, lambda = lambda)
  return(out)

}

iv_interaction <- function(Y,A,X1,X2,A.star,Z) {
  #' Given input data/valid instrument Z from gen_data(), outputs the results
  #' of a 2SLS estimation of the assumed DGP model
  #' 
  #'
  #'
  
  # Fit first stage
  Ahat <- lm(A.star ~ X1 + X2 + Z)
  
  # Second stage
  
  
  

}



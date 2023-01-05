ideal_method <- function(Y,A,X1,X2,
                         sl.lib,
                         outcome_family=gaussian()) {
  #' Estimates ATE in ideal scenario: with full access to true values of the
  #' exposure A
  #'

  X <- cbind(1,X1,X2)
  
  # Use erf() function to est ATE
  res <- erf(y = Y, a1 = 1, a0 = 0, a = A, x = X[,-1],
  family = outcome_family, sl.lib = sl.lib)
  
  return(list(
    estimate=res$estimate,
    CI=c(res$estimate - 1.96*sqrt(res$variance),res$estimate + 1.96*sqrt(res$variance)),
    se=sqrt(res$variance)
  ))
}

naive_method <- function(Y,A.star,X1,X2,
                         sl.lib,
                         outcome_family=gaussian()) {
  #' Naively estimates ATE, using mismeasured A.star in place of true measurements
  #' A
  #'
  
  X <- cbind(1,X1,X2)
  
  # Use erf() function to est ATE
  res <- erf(y = Y, a1 = 1, a0 = 0, a = A.star, x = X[,-1],
             family = outcome_family, sl.lib = sl.lib)
  
  return(list(
    estimate=res$estimate,
    CI=c(res$estimate - 1.96*sqrt(res$variance),res$estimate + 1.96*sqrt(res$variance)),
    se=sqrt(res$variance)
  ))
}


rcal <- function(Y,A,X1,X2,A.star,v_idx,sl.lib) {
  #' Regression calibration
  #' 
  #' INPUTS: 
  #' - Y, A, X1, X2, A.star: covariates from gen_data()
  #' - v_idx: Indicator for validation data subset
  #' - sl.lib: SuperLearner libraries
  #' 
  #' OUPUTS:
  #' - A list containing rcal ATE estimate, SE and CI
  
  # Validation data
  A.sub <- A ; A.sub[v_idx==0] <- NA
  
  # Covariates
  X <- cbind(1, X1, X2)
  
  # Cal model
  A.tilde <- pred(a = A.sub, z = A.star, x = X[,-1], sl.lib = sl.lib)
  results_rc <- erf(y = Y, a1 = 1, a0 = 0, a = A.tilde, x = X[,-1],
                    family = gaussian(), sl.lib = sl.lib)
  
  # bias <- results_rc$estimate - true_effect
  se_rc <- sqrt(results_rc$variance)
  # coverage_rc <- (results_rc$estimate - 1.96*se_rc) < true_effect &
  #   (results_rc$estimate + 1.96*se_rc) > true_effect
  
  res <- list(estimate=results_rc$estimate,
              CI=c(results_rc$estimate - 1.96*se_rc, results_rc$estimate + 1.96*se_rc),
              se=se_rc
  )
  return(res)
}

csme_linear <- function(Y,A,X1,X2,A.star,
                        sig2_me,
                        init_guess=1,
                        guess_acc=0.25) {
  #' Function for implementing a modified version of the Blette (2022) CSME
  #' approach
  #'
  #' INPUTS:
  #' - Y,A,X1,X2,A.star: Simulation vars from gen_data()
  #' - true_effect: True ATE (needed for setting initial guess)
  #' - guess_acc: SD of initial guess for m-estimation
  
  ## Conditional Score Measurement Error
  X <- cbind(1,X1,X2)
  
  # Initial guess for causal effect (good but variable intuition)
  guess <- init_guess + rnorm(1, 0, 0.1)

  # Fit regression to use for starting values
  mod <- lm(Y ~ 0 + X + X:A.star)
  
  # Fit GPS model/get weights
  denom_mod <- lm(A.star ~ 0 + X)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(A.star, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(A.star ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(A.star, p_num, summary(denom_mod)$sigma)
  wts <- dens_num / dens_denom
  
  # Fit weighted regression for starting values
  wmod <- lm(Y ~ A.star, weights = wts)
  
  # Set up data frame for m_estimate
  data <- data.frame("Y" = Y, "A.star" = A.star, "wts" = wts, "X" = X)
  
  startvec <- unname(c(mean(A.star), coef(denom_mod), coef(mod), sigma(mod)^2, coef(wmod)[2]))
  results_csme <- m_estimate(estFUN = csme_aipw, data = data,
                             inner_args = list(tau2 = sig2_me, a0 = 0, a1 = 1),
                             root_control = setup_root_control(start = startvec))
  
  idx <- length(results_csme@estimates)
  
  # Get results of interest
  # bias_csme <- coef(results_csme)[idx] - true_effect
  se_csme <- sqrt(vcov(results_csme)[idx, idx])
  # coverage_csme <- (coef(results_csme)[idx] - 1.96*se_csme) < true_effect &
    # (coef(results_csme)[idx] + 1.96*se_csme) > true_effect
  CI <- c(coef(results_csme)[idx] - 1.96*se_csme, coef(results_csme)[idx] + 1.96*se_csme)
  
  return(list(
    estimate=coef(results_csme)[idx],
    CI=CI
  ))
  
}

simex_direct <- function(z, y, x, a0, a1, family = gaussian(),
                  n.boot = 100, degree = 2, mc.cores = 3,
                  tau2, lambda = seq(0.1, 2.1, by = 0.25)) {
  #' Direct simex correction
  #' z = mismeasured exposure
  #' x = covariate data
  #' a1 - a0 = points to be contrasted
  #' degree = polynomial function to extrapolate
  #' mc.cores = faster SIMEX with multicore processing?
  #' n.boot = number of bootstrap resamples
  #' tau2 = estimated measurement error variance
  
  l.vals <- lapply(lambda, function(lam, z, y, x, a0, a1, sl.lib, ...){
    
    z.mat <- replicate(n.boot, z + sqrt(lam)*rnorm(length(z), 0, sqrt(tau2)))
    
    vals <- apply(z.mat, 2, erf, y = y, a1 = a1, a0 = a0,
                  x = x, family = gaussian(), sl.lib = sl.lib)
    
    mu.vals <- do.call(c, lapply(vals, function(o) o$estimate))
    sig.vals <- do.call(c, lapply(vals, function(o) o$variance))
   
    
    s.hat <- var(mu.vals)
    sig.hat <- mean(sig.vals)
    
    return(list(estimate = mean(mu.vals), variance = sig.hat - s.hat))
    
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



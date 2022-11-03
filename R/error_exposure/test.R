### Test RC/DR estimator

rm(list = ls())

## Preliminaries
library(SuperLearner)
library(parallel)
library(abind)
library(geex)

# Code for generating and fitting data
source("~/Github/ME-Causal-Review/R/error_exposure/erf.R")
source("~/Github/ME-Causal-Review/R/error_exposure/rc.R")
source("~/Github/ME-Causal-Review/R/error_exposure/csme.R")
source("~/Github/ME-Causal-Review/R/error_exposure/simex.R")

# simulation arguments
n.sim <- 100
n <- 1500
beta1_true <- 0.7
beta3_true <- -0.7
beta5_true <- 0.4
X1_prob <- 0.5
X2_mean <- 1
true_effect <- beta1_true + beta3_true*X1_prob + beta5_true*X2_mean
sig2_me <- 0.16

# model arguments
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction") # SuperLearner Libraries

start <- Sys.time()

out <- mclapply(1:n.sim, function(i, ...) {

  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, X2_mean, 0.5)
  X <- cbind(1, X1, X2)
  A <- rnorm(n, 2 + 0.9*X1 - 0.6*X2, 1.1)
  Y_mean <- 1.5 + beta1_true*A + 0.9*X1 + beta3_true*A*X1 - 0.6*X2 + beta5_true*A*X2
  Y <- rnorm(n, Y_mean, 0.4)
  A.star <- A + rnorm(n, 0, sqrt(sig2_me))
  A.sub <- rbinom(n, 1, 0.2)*A
  A.sub[A.sub == 0] <- NA
  
  ## Regression Calibration
  
  A.tilde <- pred(a = A.sub, z = A.star, x = X[,-1], sl.lib = sl.lib)
  results_rc <- erf(y = Y, a1 = 1, a0 = 0, a = A.tilde, x = X[,-1],
                    family = gaussian(), sl.lib = sl.lib)
  
  bias_rc <- results_rc$estimate - true_effect
  se_rc <- sqrt(results_rc$variance)
  coverage_rc <- (results_rc$estimate - 1.96*se_rc) < true_effect &
    (results_rc$estimate + 1.96*se_rc) > true_effect
  
  ## Conditional Score Measurement Error
  
  # Initial guess for causal effect (good but variable intuition)
  guess <- true_effect + rnorm(1, 0, 0.1)
  
  # Fit regression to use for starting values
  mod <- lm(Y ~ 0 + X + X:A.star)
  
  denom_mod <- lm(A.star ~ 0 + X)
  p_denom <- predict(denom_mod, type='response')
  dens_denom <- dnorm(A.star, p_denom, summary(denom_mod)$sigma)
  num_mod <- lm(A.star ~ 1)
  p_num <- predict(num_mod, type='response')
  dens_num <- dnorm(A.star, p_num, summary(denom_mod)$sigma)
  wts <- dens_num / dens_denom
  
  # Fit weighted regression for starting values
  wmod <- lm(Y ~ A.star, weights = wts)
  
  data <- data.frame("Y" = Y, "A.star" = A.star, "wts" = wts, "X" = X)
  
  startvec <- unname(c(mean(A.star), coef(denom_mod), coef(mod), sigma(mod)^2, coef(wmod)[2]))
  results_csme <- m_estimate(estFUN = csme_aipw, data = data,
                        inner_args = list(tau2 = sig2_me, a0 = 0, a1 = 1),
                        root_control = setup_root_control(start = startvec))
  
  idx <- length(results_csme@estimates)
  
  bias_csme <- coef(results_csme)[idx] - true_effect
  se_csme <- sqrt(vcov(results_csme)[idx, idx])
  coverage_csme <- (coef(results_csme)[idx] - 1.96*se_csme) < true_effect &
    (coef(results_csme)[idx] + 1.96*se_csme) > true_effect
  
  ## SIMEX
  
  results_simex <- simex(z = A.star, y = Y, x = X[,-1], a0 = 0, a1 = 1, family = gaussian(), 
                         n.boot = 100, degree = 2, mc.cores = 1,  tau2 = sig2_me, 
                         lambda = seq(0.1, 2.1, by = 0.25))
  
  bias_simex <- results_simex$estimate - true_effect
  se_simex <- sqrt(results_simex$variance)
  coverage_simex <- (results_simex$estimate - 1.96*se_simex) < true_effect &
    (results_simex$estimate + 1.96*se_simex) > true_effect
  
  # Bias
  bias <- c(bias_rc, bias_csme, bias_simex)
  
  # standard errors
  se <- c(se_rc, se_csme, se_simex)
  
  # coverage
  coverage <- c()
  
  return(list(bias = bias, se = se, coverage = coverage))
  
}, mc.cores = 25, mc.preschedule = TRUE)

stop <- Sys.time()
stop - start

bias <- rbind(lapply(out, function(lst, ...) lst$bias), along = 3)
se <- rbind(lapply(out, function(lst, ...) lst$se), along = 3)
coverage <- rbind(lapply(out, function(lst, ...) lst$se), along = 3)


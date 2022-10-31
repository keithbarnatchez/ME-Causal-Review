rm(list = ls())
library(geex)
library(rootSolve)
library(parallel)
source("~/Github/ME-Causal-Review/R/error_exposure/csme.R")

# Set parameter values
n.sim <- 100
n <- 1500
beta1_true <- 0.7
beta3_true <- -0.7
beta5_true <- 0.4
X1_prob <- 0.5
X2_mean <- 1
true_effect <- beta1_true + beta3_true*X1_prob + beta5_true*X2_mean
sigma_me <- 0.16

out <- mclapply(1:n.sim, function(i, ...) {
  
  # Generate data
  X1 <- rbinom(n, 1, X1_prob)
  X2 <- rnorm(n, X2_mean, 0.5)
  X <- cbind(1, X1, X2)
  A <- rnorm(n, 2 + 0.9*X1 - 0.6*X2, 1.1)
  Y_mean <- 1.5 + beta1_true*A + 0.9*X1 + beta3_true*A*X1 - 0.6*X2 + beta5_true*A*X2
  Y <- rnorm(n, Y_mean, 0.4)
  A.star <- A + rnorm(n, 0, sqrt(sigma_me))
  
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
  results <- m_estimate(estFUN = csme_aipw, data = data,
                        inner_args = list(tau2 = sigma_me, a = 1),
                        root_control = setup_root_control(start = startvec))
  
  p <- length(results@estimates)
  bias <- coef(results)[p] - true_effect
  se <- sqrt(vcov(results)[p, p])
  coverage <- (coef(results)[p] - 1.96*se) < true_effect &
    (coef(results)[p] + 1.96*se) > true_effect
  
 out = list(bias = bias, se = se, coverage = coverage)
 return(out)
  
}, mc.cores = 25, mc.preschedule = TRUE)

bias_est <- mean(do.call(c, lapply(out, function(lst, ...) lst$bias)))
cov_est <- mean(do.call(c, lapply(out, function(lst, ...) lst$coverage)))


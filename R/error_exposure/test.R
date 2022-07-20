### Test RC/DR estimator

rm(list = ls())

## Preliminaries

library(ranger)
library(earth)
library(glmnet)
library(SuperLearner)
library(parallel)
library(abind)

# Code for generating and fitting data
source("~/Github/ME-Causal-Review/R/error_exposure/gen_data.R")
source("~/Github/ME-Causal-Review/R/error_exposure/erf.R")
source("~/Github/ME-Causal-Review/R/error_exposure/rc.R")

# simulation arguments
n.sim <- 100 # number of simulations
omega <- 1 # outcome sd
sigma <- sqrt(2) # exposure sd
tau <- 1 # error sd
mult <- 5 # for replication - the average number of replicate measurements
prob <- 0.1 # for validation - proportion of data where A is observed
n <- 800 # sample size
bw <- 0.5 # bandwidth

# model arguments
a.vals <- seq(6, 14, by = 0.02)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction", 
            "SL.glmnet", "SL.ranger", "SL.earth") # SuperLearner Libraries

start <- Sys.time()

out <- mclapply(1:n.sim, function(i, ...) {
    
  ## replication
  
  dat_rep <- gen_data(n = n, omega = omega, sigma = sigma, tau = tau, type = "replicate", mult = mult)
  
  # extract replication data
  y <- dat_rep$y
  x <- dat_rep$x
  s <- dat_rep$s
  s.id <- dat_rep$s.id
  x.id <- dat_rep$x.id
  
  # REDUNDANT!
  keep <- which((x.id %in% unique(s.id)))
  if (length(keep)!= 0) {
    y <- y[(x.id %in% keep)]
    x <- x[(x.id %in% keep),]
    x.id <- x.id[(x.id %in% keep)]
  }
  
  a.hat <- blp(s = s, s.id = s.id,x = x, x.id = x.id)$a
  rep_hat <- try(erf(y = y, a = a.hat, x = x, a.vals = a.vals, 
                     family = gaussian(), bw = bw, sl.lib = sl.lib), silent = TRUE)
  
  ## validation
  
  dat_val <- gen_data(n = n, omega = omega, sigma = sigma, tau = tau, type = "validate")
  
  # extract validation data
  y <- dat_val$y
  z <- dat_val$z
  x <- dat_val$x
  a <- dat_val$a
  id <- dat_val$id
  
  # validation subset
  a <- a*rbinom(n, 1, prob)
  a[a == 0] <- NA
  
  a.tilde <- pred(a = a, z = z, x = x, sl.lib = sl.lib)
  val_hat <- try(erf(y = y, a = a.tilde, x = x, a.vals = a.vals,
                     family = gaussian(), bw = bw, sl.lib = sl.lib), silent = TRUE)
  
  ## True True
  
  erc <- predict_example(a = a.vals, x = x)
  
  # estimates
  est <- rbind(erc, if (!inherits(rep_hat, "try-error")) {rep_hat$estimate} else {rep(NA, length(a.vals))},
               if (!inherits(val_hat, "try-error")) {val_hat$estimate} else {rep(NA, length(a.vals))})
  
  # standard errors
  se <- rbind(if (!inherits(rep_hat, "try-error")) {sqrt(rep_hat$variance)} else {rep(NA, length(a.vals))},
              if (!inherits(val_hat, "try-error")) {sqrt(val_hat$variance)} else {rep(NA, length(a.vals))})
  
  return(list(est = est, se = se))
  
}, mc.cores = 8, mc.preschedule = TRUE)

stop <- Sys.time()
stop - start

est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
mu.mat <- matrix(rep(rowMeans(est[1,,]), n.sim), nrow = length(a.vals), ncol = n.sim)

# Estimate
out_est <- t(apply(est, 1, rowMeans, na.rm = T))
colnames(out_est) <- a.vals
rownames(out_est) <- c("True ERF","Replication","Validation")

# Mean Absolute Error
out_bias <- t(apply(est[2:3,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
colnames(out_bias) <- a.vals
rownames(out_bias) <- c("Replication","Validation")

# Mean Squared Error
out_mse <- t(apply(est[2:3,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
colnames(out_mse) <- a.vals
rownames(out_mse) <- c("Replication","Validation")

# Coverage Probability
cp <- list(as.matrix((est[2,,] - 1.96*se[1,,]) < mu.mat & (est[2,,] + 1.96*se[1,,]) > mu.mat),
           as.matrix((est[3,,] - 1.96*se[2,,]) < mu.mat & (est[3,,] + 1.96*se[2,,]) > mu.mat))
out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("Replication","Validation")

# Plots
plot(a.vals, out_est[1,], type = "l", col = "black", lwd = 2,
     main = "Exposure Response Curve", xlab = "Exposure", ylab = "Rate of Event")
lines(a.vals, out_est[2,], type = "l", col = "blue", lwd = 2, lty = 1)
lines(a.vals, out_est[3,], type = "l", col = "red", lwd = 2, lty = 1)
legend("topleft", legend=c("True ERF","Replication","Validation"),
       col= c("black", "blue", "red"), lty = c(1,1,1), lwd=2, cex=0.8)

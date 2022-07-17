### Test DR estimator

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
n.sim <- 100
omega <- 2
sigma <- sqrt(2)
tau <- 1
mult <- 5
n <- 400
bw <- 1
prob <- 0.1

# model arguments
a.vals <- seq(6, 14, by = 0.02)

start <- Sys.time()

out <- mclapply(1:n.sim, function(i, ...) {
  
  dat_rep <- gen_data(n = n, omega = omega, sigma = sigma, tau = tau, type = "replicate", mult = mult)
  dat_val <- gen_data(n = n, omega = omega, sigma = sigma, tau = tau, type = "validate")
    
  ## replicate
  
  # data
  y <- dat_rep$y
  x <- dat_rep$x
  s <- dat_rep$s
  s.id <- dat_rep$s.id
  x.id <- dat_rep$id
  
  ### Redundant!
  keep <- which((x.id %in% unique(s.id)))
  if (length(keep)!= 0) {
    y <- y[(x.id %in% keep)]
    x <- x[(x.id %in% keep),]
    x.id <- x.id[(x.id %in% keep)]
  }
  
  a.hat <- blp(s = s, s.id = s.id,x = x, x.id = x.id)$a
  
  # naive
  rep_hat <- try(erf(y = y, a = a.hat, x = x, a.vals = a.vals, bw = bw), silent = TRUE)
  
  ### validation
  
  # validation subset
  a <- dat_val$a*rbinom(n, 1, prob)
  a[a == 0] <- NA
  z <- dat_val$z
  x <- dat_val$x
  
  a.tilde <- pred(a = a, z = z, x = x,
       sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", 
                  "SL.glmnet", "SL.ranger", "SL.earth"))
  
  # real
  val_hat <- try(erf(y = y, a = a.tilde, x = x, a.vals = a.vals, bw = bw), silent = TRUE)
  
  ### true true
  
  erc <- predict_example(a = a.vals, x = x)
  
  # estimates
  est <- rbind(erc, if (!inherits(rep_hat, "try-error")) {rep_hat$estimate} else {rep(NA, length(a.vals))},
               if (!inherits(val_hat, "try-error")) {val_hat$estimate} else {rep(NA, length(a.vals))})
  
  #standard error
  se <- rbind(if (!inherits(rep_hat, "try-error")) {sqrt(rep_hat$variance)} else {rep(NA, length(a.vals))},
              if (!inherits(val_hat, "try-error")) {sqrt(val_hat$variance)} else {rep(NA, length(a.vals))})
  
  return(list(est = est, se = se))
  
}, mc.cores = 25, mc.preschedule = TRUE)

stop <- Sys.time()
stop - start

est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
mu.mat <- matrix(rep(rowMeans(est[1,,]), n.sim), nrow = length(a.vals), ncol = n.sim)

# coverage probability
cp <- list(as.matrix((est[2,,] - 1.96*se[1,,]) < mu.mat & (est[2,,] + 1.96*se[1,,]) > mu.mat),
           as.matrix((est[3,,] - 1.96*se[2,,]) < mu.mat & (est[3,,] + 1.96*se[2,,]) > mu.mat))

out_est <- t(apply(est, 1, rowMeans, na.rm = T))
colnames(out_est) <- a.vals
rownames(out_est) <- c("True ERF","Replicate","Validation")

out_bias <- t(apply(est[2:5,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
colnames(out_bias) <- a.vals
rownames(out_bias) <- c("True ERF","Replicate","Validation")

out_mse <- t(apply(est[2:5,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
colnames(out_mse) <- a.vals
rownames(out_mse) <- c("True ERF","Replicate","Validation")

out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("True ERF","Replicate","Validation")

plot(a.vals, out_est[1,], type = "l", col = "black", lwd = 2,
     main = "Exposure Response Curve", xlab = "Exposure", ylab = "Rate of Event", 
     ylim = c(0,0.25))
lines(a.vals, out_est[2,], type = "l", col = "blue", lwd = 2, lty = 1)
lines(a.vals, out_est[3,], type = "l", col = "red", lwd = 2, lty = 1)
legend(6, 0.15, legend=c("True ERF","NAIVE","RC","BART","GLM"),
       col=hue_pal()(5), lty = c(1,1,1), lwd=2, cex=0.8)


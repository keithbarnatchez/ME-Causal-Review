
## wrapper for erf.R to allow for SIMEX correction

# z = mismeasured exposure
# x = covariate data
# a1 - a0 = points to be contrasted
# degree = polynomial function to extrapolate
# mc.cores = faster SIMEX with multicore processing?
# n.boot = number of bootstrap resamples
# tau2 = estimated measurement error variance

simex <- function(z, y, x, a0, a1, family = gaussian(),
                  n.boot = 100, degree = 2, mc.cores = 3,
                  tau2, lambda = seq(0.1, 2.1, by = 0.25)) {
  
  l.vals <- mclapply(lambda, function(lam, z, y, x, a0, a1, sl.lib, ...){
    
    z.mat <- replicate(n.boot, z + sqrt(lam)*rnorm(length(z), 0, sqrt(tau2)))
    
    vals <- apply(z.mat, 2, erf, y = y, a1 = a1, a0 = a0,
                  x = x, family = gaussian(), sl.lib = sl.lib)
    
    mu.vals <- do.call(c, lapply(vals, function(o) o$estimate))
    sig.vals <- do.call(c, lapply(vals, function(o) o$variance))
    
    s.hat <- var(mu.vals)
    sig.hat <- mean(sig.vals)
    
    return(list(estimate = mean(mu.vals), variance = sig.hat - s.hat))
    
  }, z = z, y = y, x = x, a0 = a0, a1 = a1, sl.lib = sl.lib, mc.cores = mc.cores)
  
  if (any(lambda == 0)){
    
    Psi <- do.call(c, lapply(l.vals, function(o) o$estimate))[,-which(lambda == 0)]
    Phi <- do.call(c, lapply(l.vals, function(o) o$variance))[,-which(lambda == 0)]
    
  } else {
    
    Psi <- do.call(c, lapply(l.vals, function(o) o$estimate))
    Phi <- do.call(c, lapply(l.vals, function(o) o$variance))
    
  }
  
  L <- cbind(1, poly(lambda, degree = degree, raw = TRUE))
  chi <- c(1, poly(-1, degree = degree, raw = TRUE))
  estimate <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Psi)
  variance <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Phi)
  
  out <- list(estimate = estimate, variance = variance, Psi = Psi, Phi = Phi, lambda = lambda)
  
  return(out)
  
}

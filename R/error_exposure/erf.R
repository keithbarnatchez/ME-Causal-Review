# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
erf <- function(a, y, x, a0, a1, family = gaussian(),
                sl.lib = c("SL.mean", "SL.glm")) {	
  
  n <- length(a)
  weights <- rep(1, times = n) # placeholder until we can incorporate this
  a.vals <- seq(min(a), max(a), length.out = 100)  # placeholder until we can incorporate this
  
  wrap <- sl_est(y = y, a = a, x = x, a.vals = a.vals, 
                 family = family, sl.lib = sl.lib)
  
  psi <- wrap$psi
  int.mat <- wrap$int.mat
  
  # asymptotics
  out <- contrast(a0 = a0, a1 = a1, psi = psi, a = a, weights = weights, 
            se.fit = TRUE, a.vals = a.vals, int.mat = int.mat)
  
  estimate <- out[1]
  variance <- out[2]
  
  # bootstrap
  # boot.idx <- cbind(1:n, replicate(200, sample(x = n, size = n, replace = TRUE)))
  # 
  # out <- apply(boot.idx, 2, function(idx, ...) {
  #   
  #   sapply(a.vals, kern_est, psi = psi[idx], a = a[idx], 
  #          weights = weights[idx], bw = bw, se.fit = FALSE)
  #   
  # })
  # 
  # estimate <- out[,1]
  # variance <- apply(out[,2:ncol(out)], 1, var)
  
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}


# estimate glm outcome model (Keith - adapt code to fit a Poisson outcome as an excercise)
sl_est <- function(a, y, x, a.vals, family = gaussian(), sl.lib = c("SL.mean", "SL.glm"), ...) {
  
  if (is.null(weights))
    weights <- rep(1, nrow(x))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a)
  
  # outcome model
  mumod <- SuperLearner(Y = y, X = xa, family = family, SL.library = sl.lib)
  muhat <- mumod$SL.predict
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    xa.tmp <- data.frame(x, a = a.tmp)
    predict(mumod, newdata = xa.tmp)$pred
    
  })
  
  mhat.vals <- colMeans(muhat.mat)
  mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = a)$y

  # GPS model
  pimod <- SuperLearner(Y = a, X = x, SL.library = sl.lib)
  pimod.vals <- pimod$SL.predict
  pimod.sd <- sd(a - pimod.vals)
  a.std <- c(a - pimod.vals)/pimod.sd
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd # don't forget le Jacobian!
  
  phat.vals <- sapply(a.vals, function(a.tmp, ...){
    
    a.std.tmp <- c(a.tmp - pimod.vals)/pimod.sd
    mean(approx(x = dens$x, y = dens$y, xout = a.std.tmp)$y / pimod.sd, na.rm = TRUE)
    
  })
  
  phat <- predict(smooth.spline(a.vals, phat.vals), x = a)$y
  phat[phat <= 0] <- .Machine$double.eps
  
  # pseudo outcome
  psi <- c(y - muhat)*c(phat/pihat) + mhat
  
  # integration matrix
  mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
  phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)  
  int.mat <- (muhat.mat - mhat.mat)*phat.mat
  
  out <- list(psi = psi, int.mat = int.mat)
  
  return(out)
  
}

# Kernel weighted least squares
contrast <- function(a0, a1, a, psi, bw = 1, weights = NULL, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  if (is.null(weights))
    weights <- rep(1, times = n)
  
  # Regression
  g.std <- cbind(1, a)
  
  mod <- lm(psi ~ a, weights = weights)
  mu <- (a1 - a0)*mod$coefficients[2]
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- mod$fitted.values
    
    # Gaussian Kernel Matrix
    g.vals <- matrix(rep(a.vals, n), byrow = T, nrow = n)
    intfn1.mat <- int.mat
    intfn2.mat <- g.vals * int.mat
    
    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
    U <- solve(crossprod(g.std, weights*g.std))
    V <- cbind(weights*((psi - eta) + int1),
               weights*(a * (psi - eta) + int2))
    sig2 <- t(c(0, a1 - a0)) %*%U %*% crossprod(V) %*% U %*% c(0, a1 - a0)
    
    return(c(mu = mu, sig2 = sig2))
    
    
  } else
    return(mu)
  
}

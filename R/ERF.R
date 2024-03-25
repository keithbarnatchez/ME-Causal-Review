# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
erf <- function(a, y, x, a0, a1, family = gaussian(),
                sl.lib = c("SL.mean", "SL.glm")) {	
  
  n <- length(a)
  a.vals <- seq(min(a), max(a), length.out = 100)  # placeholder until we can incorporate this
  
  wrap <- sl_est(y = y, a = a, x = x, a.vals = a.vals, 
                 family = family, sl.lib = sl.lib)
  
  psi <- wrap$psi
  int.mat <- wrap$int.mat
  
  # asymptotics
  results <- contrast(a0 = a0, a1 = a1, psi = psi, a = a, se.fit = TRUE, 
                      a.vals = a.vals, int.mat = int.mat)
  
  Mhat <- results$mu
  Vhat <- results$sig2
  EIF_hat <- results$eif
  
  # Rely on asymptotics for CI construction
  lower_ci <- Mhat - sqrt(Vhat)*qnorm(0.975)
  upper_ci <- Mhat + sqrt(Vhat)*qnorm(0.975)
  
  return(list(EST = Mhat, CI = c(lower_ci, upper_ci),
              VAR = Vhat, EIF = EIF_hat))
  
  return(out)
  
}


# estimate glm outcome model
sl_est <- function(a, y, x, a.vals, family = gaussian(), sl.lib = c("SL.mean", "SL.glm"), ...) {
  
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
contrast <- function(a0, a1, a, psi, bw = 1, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  # Regression
  mod <- lm(psi ~ a)
  mu <- (a1 - a0)*mod$coefficients[2]
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- mod$fitted.values
    
    # Gaussian Kernel Matrix
    g.std <- cbind(1, a)
    g.vals <- matrix(rep(a.vals, n), byrow = TRUE, nrow = n)
    intfn1.mat <- int.mat
    intfn2.mat <- g.vals * int.mat
    
    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), 
                           byrow = T, nrow = n) *
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n),
                           byrow = T, nrow = n) *
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
    U <- solve(crossprod(g.std))
    V <- n*cbind((psi - eta) + int1, a * (psi - eta) + int2)
    eif <- c(t(c(0, a1 - a0) %*% U %*% t(V)))
    sig2 <- var(eif)/n
    
    return(list(mu = mu, sig2 = sig2, eif = eif))
    
  } else
    return(mu)
  
}

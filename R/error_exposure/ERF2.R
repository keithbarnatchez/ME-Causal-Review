# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
erf <- function(a, y, x, a0, a1, family = gaussian(),
                bw.seq = seq(0.1, 5, by = 0.1),
                sl.lib = c("SL.mean", "SL.glm")) {	
  
  n <- length(a)
  a.vals <- seq(min(a), max(a), length.out = 100)  # placeholder until we can incorporate this
  
  wrap <- sl_est(y = y, a = a, x = x, a.vals = a.vals, 
                 family = family, sl.lib = sl.lib)
  
  psi <- wrap$psi
  int.mat <- wrap$int.mat
  
  # leave-one-out cross-validation to select bandwidth
  w.fn <- function(bw) {
    
    w.vals <- NULL
    
    for (a.val in a.vals) {
      a.std <- (a - a.val) / bw
      k.std <- dnorm(a.std) / bw
      w.vals <- c(w.vals, mean(a.std^2 * k.std) * (dnorm(0) / bw) /
                    (mean(k.std) * mean(a.std^2 * k.std) - mean(a.std * k.std)^2))
    }
    
    return(w.vals / n)
    
  }
  
  hatvals <- function(bw) {
    approx(a.vals, w.fn(bw), xout = a)$y
  }
  
  cts.eff.fn <- function(out, bw) {
    approx(locpoly(a, out, bandwidth = bw), xout = a)$y
  }
  
  # note: choice of bandwidth range depends on specific problem,
  # make sure to inspect plot of risk as function of bandwidth
  risk.fn <- function(h) {
    hats <- hatvals(h)
    mean(((psi - cts.eff.fn(psi, bw = h)) / (1 - hats))^2)
  }
  
  risk.est <- sapply(bw.seq, risk.fn)
  bw <- bw.seq[which.min(risk.est)]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  
  # asymptotics
  results <- lapply(list(a0, a1), contrast, psi = psi, a = a, se.fit = TRUE, 
                    a.vals = a.vals, int.mat = int.mat)
  
  Mhat <- results[[2]]$mu - results[[1]]$mu
  EIF_hat <- results[[2]]$eif - results[[1]]$eif
  Vhat <- var(EIF_hat)/n
  
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
contrast <- function(a.new, a, psi, bw = 1, se.fit = FALSE, a.vals = NULL, int.mat = NULL) {
  
  n <- length(a)
  
  a.std <- (a - a.new) / bw
  g.std <- cbind(1, a.std)
  k.std <- dnorm(a.std) / bw
  mod <- lm(psi ~ a.std, weights = k.std)
  mu <- mod$coefficients[1]
  
  if (se.fit & !is.null(int.mat) & !is.null(a.vals)) {
    
    eta <- mod$fitted.values
    
    # matrices
    k.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = TRUE, nrow = n)
    g.mat <- matrix(rep((a.vals - a.new) / bw, n), byrow = TRUE, nrow = n)
    
    # Gaussian Kernel Matrix
    intfn1.mat <- k.mat * int.mat
    intfn2.mat <- g.mat * k.mat * int.mat
    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n),  byrow = T, nrow = n) *
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
    U <- solve(crossprod(g.std))
    V <- n*cbind((psi - eta) + int1, a * (psi - eta) + int2)
    eif <- c(c(1, 0) %*% U %*% t(V))
    sig2 <- var(eif)/n
    
    return(list(mu = mu, sig2 = sig2, eif = eif))
    
  } else
    return(mu)
  
}

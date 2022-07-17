# wrapper function to fit a nonparametric ERF with measurement error using kernel-weighted regression
erf <- function(a, y, x, family = gaussian(),
                a.vals = seq(min(a), max(a), length.out = 100),
                bw = NULL, bw.seq = seq(0.1, 2, by = 0.1), folds = 5,
                sl.lib = c("SL.mean", "SL.glm")) {	
  
  n <- length(a)
  weights <- rep(1, times = n) # placeholder until we can incorporate this
    
  wrap <- sl_est(y = y, a = a, x = x, a.vals = a.vals, 
                 family = family, sl.lib = sl.lib)
  
  psi <- wrap$psi
  int.mat <- wrap$int.mat
  
  # select bw if null
  if (is.null(bw))
    bw <- cv_bw(a = a, psi = psi, weights = weights, folds = folds, bw.seq = bw.seq)
  
  # asymptotics
  out <- sapply(a.vals, kern_est, psi = psi, a = a, weights = weights,
                bw = bw, se.fit = TRUE, int.mat = int.mat, a.vals = a.vals)
  
  estimate <- out[1,]
  variance <- out[2,]
  
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
  
  names(estimate) <- names(variance) <- a.vals
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
    
    a.tmp <- c(a.tmp - pimod.vals)/pimod.sd
    mean(approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd, na.rm = TRUE)
    
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
kern_est <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE, int.mat = NULL, a.vals = NULL) {
  
  n <- length(a)
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
    
    eta <- c(g.std %*% b)
    
    # Gaussian Kernel Matrix
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.vals * kern.mat * int.mat
    
    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}

# k-fold cross validation to select bw
cv_bw <- function(a, psi, weights = NULL, folds = 5, bw.seq = seq(0.1, 2, by = 0.1)) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  idx <- sample(x = folds, size = n, replace = TRUE)
  
  cv.mat <- sapply(bw.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a[idx == k], kern_est, psi = psi.sub[idx != k], a = a.sub[idx != k], 
                      weights = weights[idx != k], bw = h, se.fit = FALSE)
      cv.vec[k] <- mean((psi.sub[idx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  bw <- bw.seq[which.min(cv.err)]
  
  return(bw)
  
}


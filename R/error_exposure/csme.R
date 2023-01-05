# alpha: regression coefficients
# beta: regression coefficients interacting with delta
# sig2: outcome model variance

# X: covariate matrix
# Y: outcome
# A.star: error prone exposure
# wts: Inverse probability weights fit with glm conditional on X
# tau2: measurement error variance

csme_aipw <- function(data){
  
  A.star <- data$A.star
  Y <- data$Y
  wts <- data$wts
  X <- as.matrix(data[,4:ncol(data)])
  
  #data dimensions
  n <- nrow(X)
  p <- ncol(X)
  
  # not confident about the denominator to this estimator
  condexp <- function(X, alpha, beta, delta, tau2, sig2) {
    (X%*%alpha + delta*(X%*%beta)) /
      (1 + (X%*%beta)^2*c(tau2 / sig2))
  }
  
  condvar <- function(X, beta, tau2, sig2) {
    sig2 / (1 + (X%*%beta)^2*c(tau2 / sig2))
  }
  
  function(theta, tau2, a0, a1){
  
    # model parameters
    gamma <- theta[2:(p + 1)]
    alpha <- theta[(p + 2):(2*p + 1)]
    beta <- theta[(2*p + 2):(3*p + 1)]
    sig2 <- theta[(3*p + 2)]
    eta <- theta[(3*p + 3)]
    
    # gps model uncertainty (questionable whether this is sufficient?)
    gps_num_eqn <- A.star - theta[1]
    gps_denom_eqn <- crossprod(X, A.star - c(X %*% gamma))
    
    # conditional score statistics
    delta <- A.star + (tau2/sig2)*Y*c(X%*%beta)
    m_A <- condexp(X = X, alpha = alpha, beta = beta, 
                    delta = delta, tau2 = tau2, sig2 = sig2)
    v_A <- condvar(X = X, beta = beta, tau2 = tau2, sig2 = sig2)
    
    # outcome model estimating equations
    ols_eqn1 <- crossprod(X, wts*(Y - m_A))
    ols_eqn2 <- crossprod(delta*X, wts*(Y - m_A))
    scl <- (Y - m_A) / sqrt(v_A)
    disp_eqn <- sig2 - sig2*crossprod(scl, wts*scl)

    # prediction estimating equation
    pred_eqn <- (a1 - a0) * c(X %*% beta) - eta
    
    c(gps_num_eqn, gps_denom_eqn, ols_eqn1, ols_eqn2, disp_eqn, pred_eqn)   
    
  }
  
}

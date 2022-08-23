outcome_me <- function(A, S, X, Y, ...) {
  
  n1 <- sum(S)
  n0 <- sum(1 - S)
  n <- n1 + n0
  m <- ncol(X)
  
  # full calibration weights
  X <- X %>% mutate_if(is.numeric, scale)
  x.mat <- model.matrix(~ ., data = data.frame(x))
  astar <- c(a_x - mean(a_x))/var(a_x)
  astar2 <- c((a_x - mean(a_x))^2/var(a_x) - 1)
  mod <- calibrate(cmat = cbind(1, x.mat*astar, astar2), 
                   target = c(n, rep(0, ncol(x.mat) + 1)))
  
  # try direct optimization
  fit_out <- try( calibrate(constraint = A,
                            target = b,
                            base_weights = base_weights,
                            coefs_init = coefs_init,
                            distance = distance,
                            optim_ctrl = optim_ctrl),
                  silent = TRUE )
  
  
  if (!inherits(fit_out, "try-error")) {
    
    weights <- fit_out$weights
    converged <- fit_out$converged
    coefs <- fit_out$coefs
    
  } else {
    stop("optimization failed")
  }
  
  if (!converged)
    warning("model failed to converge")
  
  tau <- sum(weights*(2*Z - 1)*Y)/sum(weights*Z)
  
  U <- matrix(0, ncol = 5*m, nrow = 5*m)
  v <- rep(0, times = 5*m + 1)
  meat <- matrix(0, ncol = 5*m + 1, nrow = 5*m + 1)
  
  for (i in 1:n) {
    
    U[1:(4*m),1:(4*m)] <- U[1:(4*m),1:(4*m)] - weights[i] * A[i,] %*% t(A[i,])
    U[1:m,(4*m + 1):(5*m)] <- U[1:m,(4*m + 1):(5*m)] - diag(S[i], m, m)
    U[(m + 1):(2*m),(4*m + 1):(5*m)] <- U[(m + 1):(2*m),(4*m + 1):(5*m)] - diag(S[i], m, m)
    U[(2*m + 1):(3*m),(4*m + 1):(5*m)] <- U[(2*m + 1):(3*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    U[(3*m + 1):(4*m),(4*m + 1):(5*m)] <- U[(3*m + 1):(4*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    U[(4*m + 1):(5*m),(4*m + 1):(5*m)] <- U[(4*m + 1):(5*m),(4*m + 1):(5*m)] - diag(1 - S[i], m, m)
    
    v[1:(4*m)] <- v[1:(4*m)] - weights[i]*(2*Z[i] - 1)*(Y[i] - Z[i]*tau)*A[i,]
    v[5*m + 1] <- v[5*m + 1] - weights[i]*Z[i]
    
    meat <- meat +  tcrossprod(esteq_fusion(X = X[i,], Y = Y[i], Z = Z[i], S = S[i], 
                                            p = weights[i], q = base_weights[i], 
                                            theta = theta, tau = tau))
    
  }
  
  invbread <- matrix(0, nrow = 5*m + 1, ncol = 5*m + 1)
  invbread[1:(5*m),1:(5*m)] <- U
  invbread[5*m + 1, ] <- v
  
}


# general calibration function
calibrate <- function(cmat, target, base_weights = NULL, coefs_init = NULL,
                      optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  fn <- match.fun(lagrange_ent)
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, nrow(cmat))
  } else if (length(base_weights) != nrow(cmat)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(cmat)) 
  } else if (length(coefs_init) != ncol(cmat)) {
    stop("length(coefs_init) != ncol(cmat)")
  }
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS", hessian = TRUE,
                      cmat = cmat, base_weights = base_weights,
                      target = target, control = optim_ctrl)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  weights <- c( base_weights*exp(-cmat %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}

# entropy objective function
lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}
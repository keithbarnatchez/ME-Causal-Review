# error_exposure/data_functions.R

# First define functions for quickly computing common functions
logit <- function(x) {
  return( log( x/(1-x) ) )
}

expit <- function(x) {
  return( exp(x) / (1+exp(x)) )
}
 
meas_model <- function(A, sig_u) {
  
  #' Generates error-prone measurements of exposure A, with variance sig_u
  #' Working with just classical measurement error for now, but there is option
  #' to alter code for different processes
  #' 
  #' INPUTS:
  #' - A: true exposure variable
  #' - sig_u: standard deviation of meas error
  #' 
  #' OUTPUTS: 
  #' - Error-prone exposurement measurement variable
  
  return(A + rnorm(length(A), 0, sig_u))
  
}

trt_model <- function(W, X, V, aw = .25, ax = .25,  av = .25, a0 = 0) {
  
  #' Generates treatment model for T with error-prone exposure X and 
  #' properly-measured exposure Z
  
  mu <- a0 + aw*W + ax*X + av*V
  A <- rnorm(length(X),mean = mu,sd=1) 
  return(A)
  
}

out_model <- function(A, W, X, sig_e = 1, binary = FALSE, ba = 1, bw = 1,
                          baw = 0.2, bx = 1, bax = 0.2, b0 = 0) {
  
  #' Generate Y from N(mu,sig_e) where mu is a linear function of T, X and Z
  #' INPUTS:
  #' - A, W, X: True exposure and confounders
  #' - ba, bw, bx, b0: coefficients on A, W, X, and intercept in outcome model
  #' - sig_e: variance of outcome cond on covariates (if continuous)
  #' 
  #' OUTPUTS:
  #' - a vector containing simulated values of the outcome Y, given the specified
  #'   dgp (governed by the user-supplied parameters)
  
  mu <- b0 + ba*A + bw*W + bx*X + baw*A*W + bax*A*X
  
  if (!binary) { # if continuous outcome
    return(rnorm(length(A), mean = mu, sd = sig_e))
  }  else { #if binary outcome
    testit::assert((ba < 1) & (0 < boa) ) # make sure effect of a is < 1
    p <- (1 - ba)*expit(b0 + bw*W + bx*X) + ba*A
    return(rbinom(length(X), size=1, prob = p)  )
  }
  
}

generate_covariates <- function(n, rho, psi) {
  
  #' Generates X, Z, and an instrument V from a multivariate normal distribution
  #' It is assumed each of X, Z and V are marginally N(0,1), and that 
  #' corr(X,Z)=rho, corr(X,V)=psi, and that corr(Z,V) = 0
  #'  
  #' INPUTS:
  #' - rho: correlation between X and Z
  #' - psi: correlation between X and V
  #' 
  #' OUTPUTS:
  #' - Matrix (X, Z, V) of simulated values following a MVN dist as specified
  #'  by user
  
  return(mvtnorm::rmvnorm(n = n, sigma = matrix(c(1 ,  rho,  psi,
                                              rho,  1 ,   0,
                                              psi,  0 ,   1),
                                            nrow = 3, byrow = T)))
  
}

generate_data <- function(n, sig_e = 1, sig_u = 0.1,
                          rho = 0.5, psi = -0.25,
                          aw = 0.5, ax = -0.5, av = 0.25, a0 = 0,
                          ba = 1, bw = -1, bx = 0.5, b0 = 0,
                          baw = 0.25, bax = -0.25,
                          v_share = 0.1, binary = FALSE) {
  
  #' Generates dataset with outcome variable y, error-prone exposure X with 
  #' measurements W, binary treatment of interest T, and confounding variable Z.
  #' For now, I'm assuming homoskedasticity in the outcome and measurement error
  #' models
  #' 
  #' INPUTS:
  #' - n: sample size
  #' - sig_e: outcome model variance
  #' - sig_u: measurement error variance
  #' - rho: correlation b/w error-prone exposure and error-free covariate 
  #' - psi: correlation b/w error-prone exposure and instrumental variable V
  #' - ax, az, a0: coefficients on W, X and intercept in PS model
  #' - ba, bw, bx, b0: coefficients on A, W, X and intercept in outcome model
  #' 
  #' OUTPUTS:
  #' - simulated dataset
  #' 
  
  # Simulate covariates
  covariates <- generate_covariates(n, rho, psi)
  W <- covariates[,1] ; X <- covariates[,2] ; V <- covariates[,3]
  
  # Simulate treatment process
  A <- trt_model(W = W, X = X, V = V, aw = aw, ax = ax, av = av, a0 = a0)
  
  # Simulate outcome 
  Y <- out_model(A = A, W = W, X = X, 
                 ba = ba, bw = bw, baw = baw,
                 bx = bx, bax = bax, b0 = b0, 
                 sig_e = sig_e, binary = binary)
  
  # Simulate error-prone measurements
  A.star <- meas_model(A, sig_u) 
  
  # Get validation data index
  val.idx <- rbinom(n, size = 1, prob = v_share)
  
  out_data <- data.frame(Y = Y, A = A, A.star = A.star,
                         W = W, X = X, V = V, 
                         val.idx = val.idx)
  
  return(out_data)
  
}
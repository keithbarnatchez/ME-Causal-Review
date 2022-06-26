logit <- function(x) {
  return( log( x/(1-x) ) )
}

expit <- function(x) {
  return( exp(x) / (1+exp(x)) )
}

meas_model <- function(X,sig_u) {
  #' Generates error-prone measurements of exposure X, with variance sig_u
  #' Working with just classical measurement error for now, but there is option
  #' to alter code for different processes
  #' 
  #' INPUTS:
  #' - X: true exposure variable
  #' - sig_u: standard deviation of meas error
  #' 
  #' OUTPUTS: 
  #' - Error-prone exposurement measurement variable
  
  return( X + rnorm( length(X),0,sig_u ) )
}

tmt_model <- function(X,Z,
                      ax = .5,
                      az = .5,
                      a0 = 0) {
  #' Generates treatment model for T with error-prone exposure X and 
  #' properly-measured exposure Z
  
  p_tmt <- expit(a0 + ax*X + az*Z)
  A <- rbinom(length(X),size=1,prob=p_tmt) 
  return(A)
}

outcome_model <- function(A,X,Z,
                          ba = 1,
                          bx = 1,
                          bz = 1,
                          b0 = 0,
                          sig_e = 1,
                          binary=0) {
  #' Generate Y from N(mu,sig_e) where mu is a linear function of T, X and Z
  
  mu <- b0 + ba*A + bx*X + bz*Z
  if (binary==0) { # if continuous outcome
    return( rnorm(length(X), mean = mu, sd=sig_e))
  }
  else { #if binary outcome
    return( rbinom(length(X),size=1,prob=expit(mu))  )
  }
}

generate_covariates <- function(n, rho, psi) {
  #' Generate X and Z from multivariate normal dist with specified covariance
  #' rho
  
  return( rmvnorm(n, sigma = matrix(c(1,  rho, psi,
                                      rho,  1,   0,
                                      psi,  0,   1),nrow=3,byrow=T)) )
}

generate_data <- function(n,
                       sig_e=1,
                       sig_u=0.1,
                       rho=0.5, psi=0.5,
                       ax=0.5,az=0.5,a0=0,
                       ba=1,bx=1,bz=1,b0=0,
                       v_share=0.1,
                       binary=0) {
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
  #' - ax, az, a0: coefficients on x, z and intercept in PS model
  #' - ba, bx, bz, b0: coefficients on a, x, z and intercept in outcome model
  #' 
  #' OUTPUTS:
  #' - simulated dataset
  #' 
  
  # Simulate covariates
  covariates <- generate_covariates(n,rho,psi)
  X <- covariates[,1] ; Z <- covariates[,2] ; V <- covariates[,3]
  
  # Simulate treatment process
  A <- tmt_model(X,Z,ax,az,a0)
  
  # Simulate outcome 
  Y <- outcome_model(A,X,Z,ba,bx,bx,b0,sig_e,binary)
  
  # Simulate error-prone measurements
  W <- meas_model(X, sig_u) 
  
  # Get validation data index
  v_share <- 0.1 # share in validation data
  v_idx <- rbinom(n,size=1,prob=v_share)
    
  out_data <- data.frame(Y=Y,X=X,A=A,W=W,Z=Z,V=V,v_idx=v_idx)
  return(out_data)
}
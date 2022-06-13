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
  T <- rbinom(length(X),size=1,prob=p_tmt) 
  return(T)
}

outcome_model <- function(T,X,Z,
                          bt = 1,
                          bx = 1,
                          bz = 1,
                          b0 = 0,
                          sig_e = 1) {
  #' Generate Y from N(mu,sig_e) where mu is a linear function of T, X and Z
  
  mu <- b0 + bt*T + bx*X + bz*Z
  
  return( rnorm(length(X), mean = mu, sd=sig_e))
  
}

generate_covariates <- function(n, rho) {
  #' Generate X and Z from multivariate normal dist with specified covariance
  #' rho
  
  return( rmvnorm(n, sigma = matrix(c(1,rho,rho,1),nrow=2,byrow=T)) )
}

generate_data <- function(n,
                       sig_e=1,
                       sig_u=0.1,
                       rho=0.5,
                       ax=0.5,az=0.5,a0=0,
                       bt=1,bx=1,bz=1,b0=0) {
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
  #' - ax, az, a0: coefficients on x, z and intercept in PS model
  #' - bt, bx, bz, b0: coefficients on t, x, z and intercept in outcome model
  #' 
  #' OUTPUTS:
  #' - simulated dataset
  
  # Simulate covariates
  covariates <- generate_covariates(n,rho)
  X <- covariates[,1] ; Z <- covariates[,2]
  
  # Simulate treatment process
  T <- tmt_model(X,Z,ax,az,a0)
  
  # Simulate outcome 
  Y <- outcome_model(T,X,Z,bt,bx,bx,b0,sig_e)
  
  # Simulate error-prone measurements
  W <- meas_model(X, sig_u) 
    
  out_data <- data.frame(Y=Y,X=X,T=T,W=W,Z=Z)
  
}



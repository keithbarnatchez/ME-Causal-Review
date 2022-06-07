# ----------------------
rm(list=ls())
library(mvtnorm)
1
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
                      bx = .5,
                      bz = .5) {
  #' Generates treatment model for T with error-prone exposure X and 
  #' properly-measured exposure Z
  
  p_tmt <- expit(bx*X + bz*Z)
  T <- rbinom(length(X),size=1,prob=p_tmt) 
  return(T)
}

outcome_model <- function(T,X,Z,
                          bt = 1,
                          bx = 1,
                          bz = 1,
                          b0 = 0,
                          sig_e = 1) {
  
  mu <- b0 + bt*T + bx*X + bz*Z
  
  return( rnorm(length(X), mean = mu, sd=sig_e))
  
}

generate_covariates <- function(n, rho) {
  return( rmvnorm(n, sigma = matrix(c(1,rho,rho,1),nrow=2,byrow=T)) )
}

generate_data <- function(n,
                       sig_e=1,
                       sig_u=1.1,
                       rho=0.5) {
  #' Generates dataset with outcome variable y, error-prone exposure X with 
  #' measurements W, binary treatment of interest T, and confounding variable Z.
  #' For now, I'm assuming homoskedasticity in the outcome and measurement error
  #' models
  #' 
  #' INPUTS:
  #' - n: sample size
  #' - sig_e: outcome model variance
  #' - sig_u: measurement error variance
  #' - rho: correlation b/w error-prone exposure and 
  #' 
  #' OUTPUTS:
  #' 
  
  # Simulate covariates
  covariates <- generate_covariates(n,rho)
  X <- covariates[,1] ; Z <- covariates[,2]
  
  # Simulate treatment process
  T <- tmt_model(X,Z)
  
  # Simulate outcome 
  Y <- outcome_model(T,X,Z)
  
  # Simulate error-prone measurements
  W <- meas_model(X, sig_u) 
    
  out_data <- data.frame(Y=Y,X=X,T=T,W=W,Z=Z)
  
}



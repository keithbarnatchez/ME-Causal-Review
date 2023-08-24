gen_corr_norm_bern <- function(n,
                               rho,
                               p=0.5) {
  #' Given a sample size b and desired correlation rho, returns a std normal
  #' rv Y and mean p bernoulli rv B with correlation rho
  
  # Find 1-pth quantile of X1
  x0 <- qnorm(1-p)
  
  # Generate X1, X2 to be ind std normal
  X1 <- rnorm(n) ; X2 <- rnorm(n)
  
  # Define r for setting correlation
  r <- rho*sqrt(p*(1-p))/dnorm(x0)
  
  # Define Y 
  Y <- r*X1 + sqrt(1-r^2)*X2
  
  # Define B
  B <- ifelse(X1>x0,1,0)   
  
  return(cbind(B,Y))
  
}


gen_data <- function(n=1500,vshare=0.1,
                     b0=0,bA=1,bX1=0.7,bX2=-0.7, bA_X1=0.2, bA_X2=0.1,
                     muX1=0.5, muX2=1, muZ=1,
                     a0=2,aX1=0.9,aX2=-0.6,aZ=0.5,
                     sigU=0.3, sigA=1,sigE=1) {
  
  #' Code for simulating data with mismeasured exposure A
  #'
  #' Variables:
  #' A: True exposure values
  #' X1: BINARY confounder
  #' X2: Continuous confonder
  #' Z: continous instrument
  #' Y: outcome
  #' Astar: error-prone exposure measurements
  #'
  #' DGP:
  #' A ~ norm(a0 + aX1(X1) + aX2(X2) + aZ(Z), sigA^2)
  #' Y ~ norm(b0 + bA(A) + bX1(X1) + bX2(X2) +bAX1(snfksd))
  #' Astar ~ norm(A, sigU^2+sigE^2)
  
  # Simulate confounders
  # X1 <- rbinom(n,1,muX1) ; X2 <- rnorm(n,muX2,1)
  X <- gen_corr_norm_bern(n,rho=0.5,p=muX1)
  X1 <- X[,1] ; X2 <- X[,2]+muX2
  
  # Simulate instrument
  Z <- rnorm(n,muZ,1)
  
  # Simulate exposure
  A <- rnorm(n,a0 + aX1*X1 + aX2*X2 + aZ*Z,sigA)
  
  # Simulate outcome
  Y <- rnorm(n, 
             b0 + bA*A + bX1*X1 + bX2*X2 + bA_X1*A*X1 + bA_X2*A*X2,
             sigE)
  
  # Simulate exposure measurements
  Astar <- A + rnorm(n,mean = 0,sd = sigU)
  
  # Get the validation idx
  n_v <- floor(n*vshare)
  v_idx <- sample( c(rep(1,n_v),rep(0,n-n_v)) )
  
  # Form the df
  data <- data.frame(A=A,Y=Y,X1=X1,X2=X2,Z=Z,Astar=Astar,v_idx=v_idx)
  
  return(data)
}
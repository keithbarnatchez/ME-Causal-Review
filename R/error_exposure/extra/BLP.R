### Replicate Measurements

# s = replicate of mismeasured variable a
# s.id = id of replicate measures; unique(s.id) should be the length of a
# x = covariate measurements for unique s.id values
# x.id = id measures for x (unique(s.id) == x.id, for example)

blp <- function(s, s.id, x = NULL, x.id = NULL) {
  
  wts <- c(unname(table(s.id)))
  
  # initialize exposures
  tmp <- aggregate(s, by = list(s.id), mean)
  id <- c(tmp[,1])
  z <- c(tmp[,2])
  
  if (!is.null(x.id)) {
    
    if(sum(duplicated(x.id)) != 0)
      stop("duplicate x.id detected")
    
    if (nrow(x) != length(x.id))
      stop("x provided does not match length of unique x.id")
    
    x <- x[order(x.id),]
    id.check <- x.id[order(x.id)]
    
    if (!identical(id, id.check))
      stop("levels of z.id provided does not match unique values of x.id")
    
  }
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  
  ord <- order(s.id)
  z_s <- rep(z, wts)[order(ord)]
  
  if (!is.null(x)) {
    
    x <- as.matrix(x)
    p <- ncol(x)
    
    mu_z <- sum(wts*z)/m
    mu_x <- colMeans(x)
    mat_x <- matrix(rep(mu_x, n), nrow = n, byrow = TRUE)
    nu <- m - sum(wts^2)/m
    
    tau2 <- sum((s - z_s)^2)/(m - n)
    sigma2 <- c(sum(wts*(z - mu_z)^2) - (n - 1)*tau2)/nu
    psi <- c(crossprod(as.matrix(wts*(z - mu_z)), as.matrix(x - mat_x)))/nu
    Omega <- cov(x)
    
    phi <- c(sigma2, psi)
    V <- matrix(NA, p + 1, p + 1)
    V[2:(p+1),] <- V[,2:(p+1)] <- psi 
    V[2:(p+1),2:(p+1)] <- Omega
    
    a <- sapply(1:n, function(i, ...) {
      
      V[1,1] <- sigma2 + tau2/wts[i]
      star <- c(z[i] - mu_z, t(x[i,]) - mu_x)
      out <- c(mu_z + t(phi)%*%solve(V)%*%star)
      
      return(out)
      
    })
    
    cvar <- sapply(1:n, function(i, ...) {
      
      V[1,1] <- sigma2 + tau2/wts[i]
      out <- c(sigma2 - t(phi)%*%solve(V)%*%phi)
      
      return(out)
      
    })
    
  } else {
    
    mu_z <- sum(wts*z)/m
    nu <- m - sum(wts^2)/m
    tau2 <- sum((s - z_s)^2)/(m - n)
    sigma2 <- (sum(wts*(z - mu_z)^2) - (n - 1)*tau2)/nu
    
    a <- sapply(1:n, function(i, ...) {
      
      out <- c(mu_z + sigma2/(sigma2 + tau2/wts[i])*c(z[i] - mu_z))
      
      return(out)
      
    })
    
    cvar <- sapply(1:n, function(i, ...) {
      
      out <- c(sigma2 - sigma2^2/(sigma2 + tau2/wts[i]))
      return(out)
      
    })
    
  }
  
  return(data.frame(id = id, a = a, cvar = cvar))
  
}

# same as blp() but for multi-valued s; dim(z) = m x p

multi_blp <- function(s, s.id, x = NULL, x.id = NULL) {
  
  wts <- c(unname(table(s.id)))
  
  # initialize exposures
  tmp <- aggregate(s, by = list(s.id), mean)
  id <- c(tmp[,1])
  z <- as.matrix(tmp[,-1])
  
  if (!is.null(x.id)) {
    
    if (nrow(x) != length(id))
      stop("x provided does not match length of unique s.id")
    
    if(sum(duplicated(id)) != 0)
      stop("duplicate id detected")
    
    x <- x[order(x.id),]
    id.check <- x.id[order(x.id)]
    
    if (!identical(id, id.check))
      stop("id provided does not match unique values of s.id")
    
  }
  
  # dimensions
  m <- length(z.id)
  n <- length(id)
  
  ord <- order(z.id)
  z_tmp <- z[rep(1:nrow(z), wts),]
  z_s <- z_tmp[order(ord),]
  
  p <- ncol(z)
  q <- ncol(x)
  
  if (!is.null(x)) {
    
    mu_z <- colSums(wts*z)/m
    mu_x <- colMeans(x)
    
    mat_z <- matrix(rep(mu_z, n), nrow = n, byrow = TRUE)
    mat_x <- matrix(rep(mu_x, n), nrow = n, byrow = TRUE)
    nu <- m - sum(wts^2)/m
    
    Omega <- as.matrix(crossprod(as.matrix(s - z_s)))/(m - n)
    Sigma <- as.matrix(crossprod(c(wts)*(z - mat_z), (z - muMat_z)) - (n - 1)*Omega)/nu
    Psi <- as.matrix(crossprod(c(wts)*(z - mat_z), as.matrix(x - muMat_x)))/nu
    Tau <- cov(x)
    
    Phi <- cbind(Sigma, Psi)
    V <- matrix(NA, p + q, p + q)
    V[(p+1):(p+q),1:p] <- t(Psi)
    V[1:p,(p+1):(p+q)] <- Psi
    V[(p+1):(p+q),(p+1):(p+q)] <- Tau
    
    a <- t(sapply(1:n, function(i, ...) {
      
      V[1:p,1:p] <- Sigma + Omega/wts[i]
      star <- c(t(z[i,]) - mu_z, t(x[i,]) - mu_x)
      out <- c(mu_z + Phi %*% chol2inv(chol(V)) %*% star)
      
      return(out)
      
    }))
    
  } else {
    
    mu_z <- colSums(wts*z)/m
    mat_z <- matrix(rep(mu_z, n), nrow = n, byrow = TRUE)
    nu <- m - sum(wts^2)/m
    
    Omega <- crossprod(as.matrix(s - z_s))/(m - n)
    Sigma <- as.matrix(crossprod(wts*(z - mat_z), (z - mat_z)) - (n - 1)*Omega)/nu
    
    a <- t(sapply(1:n, function(i, ...) {
      
      V <- Sigma + Omega/wts[i]
      out <- c(mu_z + Sigma %*% solve(V) %*% c(t(z[i,]) - mu_z))
      
      return(out)
      
    }))
    
  }
  
  colnames(a) <- colnames(z)
  return(data.frame(id = id, a))
  
}
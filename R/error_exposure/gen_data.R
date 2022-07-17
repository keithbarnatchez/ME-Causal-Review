## Inputs

# n = number of clusters
# mult = number of grids equals n*mult
# sigma = sd of the true exposure
# tau = exposure error
# type = "replicate" or "validate"
# mult = if type == "replicate" identifies  average # of replications

## Output

# a = true exposure
# z = error prone exposure
# s = replicate error prone exposure
# y = outcome
# id = id for a and y
# s.id = replicate id for s, same levels as id
# x = confounders

gen_data <- function(n = c(400, 800), omega = 2, sigma = sqrt(2), tau = 1, type = c("replicate", "validate"), mult = c(5, 10)) {
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  x <- cbind(x1, x2, x3, x4)
  
  mu_gps <- 10 + 0.5*x[,1] - 0.5*x[,2] - 0.5*x[,3] + 0.5*x[,4]

  a <- rnorm(n, mu_gps, sigma)
  id <- 1:n
  
  if (type == "replicate") {
    
    if (mult == 10) {
      s.id <- rep(id, rep(c(2,4,6,8,12,14,16,18), each = n/8))
    } else if (mult == 5){
      s.id <- rep(id, rep(c(1,2,3,4,6,7,8,9), each = n/8))
    }
    
    stab <- table(s.id)
    a_s <- rep(a, stab)  
    s <- rnorm(mult*n, a_s, tau)
    
  } else {
    
    z <- rnorm(n, a*(-0.25*x1 + 0.75*x2 + -0.75*x3 + 0.25*x4),  tau)
    
  }

  mu_out <- 2 - 0.5*x[,1] - 0.25*x[,2] + 0.25*x[,3] + 0.5*x[,4] +
    0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 0.25*(a - 10)*x[,1]
  
  y <- rnorm(n, mu_out) # Keith - try simulating and fitting a poisson model
  
  # create simulation dataset
  
  if (type == "replicate")
    sim <- list(a = a, s = s, y = y, x = x, id = id, s.id = s.id)
  else
    sim <- list(a = a, z = z, y = y, x = x, id = id)
    
    
  return(sim)
  
}

# get the true sample ERF
predict_example <- function(a.vals, x) {
  
  out <- rep(NA, length(a.vals))
  
  for(i in 1:length(a.vals)) {
    
    a.vec <- rep(a.vals[i],nrow(x))
    mu_out <- 2 + x %*% c(-0.5,-0.25,0.25,0.5) + 0.25*(a.vec - 10) - 
      0.75*cos(pi*(a.vec - 6)/4) - 0.25*(a.vec - 10)*x[,1]
    out[i] <- mean(mu_out)
    
  }
  
  return(out)
  
}


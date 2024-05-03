# ---------------------------------------------------------
# UTILITY FUNCTIONS
# ---------------------------------------------------------

ipw <- function(a, y, x, pihat = NULL, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction")) {
  
  if (is.null(pihat)) { # if aren't pre-supplying propensity scores
    # Get predicted propensity scores
    pimod <- SuperLearner(Y = a, X = x, SL.library = sl.lib, family = binomial())
    pihat <- pimod$SL.predict
  }
  
  xmat <- model.matrix(~ ., data = x)
  m <- ncol(xmat)
  n <- nrow(xmat)
  
  # estimate
  ipw_hat <- ifelse(a == 1, 1/pihat, 1/(1 - pihat))
  ATE_hat <- sum((2*a - 1)*ipw_hat*y)/sum(ipw_hat*a)
  
  U <- matrix(0, ncol = m, nrow = m)
  v <- rep(0, times = m + 1)
  meat <- matrix(0, ncol = m + 1, nrow = m + 1)
  
  esteq <- function(a, x, y, pihat, ipw, ATE) {
    
    eq1 <- (a - pihat)*x
    eq2 <- a*ipw*(y - ATE) - (1 - a)*ipw*y
    
    eq <- c(eq1, eq2) 
    return(eq)
    
  }
  
  for (i in 1:n) {
    
    U[1:m,1:m] <- U[1:m,1:m] - pihat[i] * tcrossprod(xmat[i,])
    v[1:m] <- v[1:m] - a[i] * (ipw_hat[i] - 1) * (y[i] - a[i]*ATE_hat) * xmat[i,] - 
      (1 - a[i]) * (ipw_hat[i] - 1) * y[i] * xmat[i,]
    v[m + 1] <- v[m + 1] - a[i]*ipw_hat[i]
    meat <- meat + tcrossprod(esteq(a = a[i], x = xmat[i,], y = y[i],
                                    ipw = ipw_hat[i], 
                                    pihat[i], ATE = ATE_hat))
    
  }
  
  invbread <- matrix(0, nrow = m + 1, ncol = m + 1)
  invbread[1:m,1:m] <- U
  invbread[m + 1,] <- v
  
  bread <- try(solve(invbread), silent = TRUE)
  
  if (inherits(bread, "try-error"))
    stop("inversion of \"bread\" matrix failed")
  
  sandwich <- bread %*% meat %*% t(bread)
  Vhat <- sandwich[m + 1, m + 1]
  
  # Rely on asymptotics for CI construction
  lower_ci <- ATE_hat - sqrt(Vhat)*qnorm(0.975)
  upper_ci <- ATE_hat + sqrt(Vhat)*qnorm(0.975)
  
  return(list(ATE = ATE_hat, VAR = Vhat, CI = c(lower_ci, upper_ci)))
  
}

aipw <- function(a, y, x, pihat = NULL, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction")) {
  
  #' Main function for implementing AIPW estimator.
  #' 
  #' INPUTS: 
  #' - data: A dataframe created by the generate_data() function
  #' - methods: Character vector with libraries to feed into SuperLearner. By
  #'            default, we just use regular glm (no non-parametric stuff)
  #'
  #' Outputs:
  #' - Estimate of the ATE
  
  # Outcome model
  outcome_model <- SuperLearner(Y = y, X = data.frame(a = a, x),
                                SL.library = sl.lib, family = gaussian())
  
  # Predict cond on A = 0, A = 1
  rhs0 <- data.frame(a = 0, x) # scenario where all tmt vals are 0
  rhs1 <- data.frame(a = 1, x) # scenario where all tmt vals are 1
  muhat0 <- predict(outcome_model, newdata = rhs0)$pred # E[Y|X,A=0]
  muhat1 <- predict(outcome_model, newdata = rhs1)$pred # E[Y|X,A=1]
  
  if (is.null(pihat)) { # if aren't pre-supplying propensity scores
    # Get predicted propensity scores
    pimod <- SuperLearner(Y = a, X = x, SL.library = sl.lib, family = binomial())
    pihat <- pimod$SL.predict
  }
  
  # Calculate the AIPW estimate and return the result
  
  # Compute relevant terms for estimator
  term1 <- a*y/pihat - (1 - a)*y/(1 - pihat)
  term2 <- (a - pihat)/(pihat*(1 - pihat))
  term3 <- (1 - pihat)*muhat1 + pihat*muhat0
  
  # ATE estimator
  ATE_hat <- mean(term1 - term2*term3)
  
  # Sandwich estimator
  EIF_hat <- (term1 - term2*term3) - ATE_hat
  Vhat <- sum(EIF_hat^2)/(length(EIF_hat)^2)
  
  # Rely on asymptotics for CI construction
  lower_ci <- ATE_hat - sqrt(Vhat)*qnorm(0.975)
  upper_ci <- ATE_hat + sqrt(Vhat)*qnorm(0.975)
  
  return(list(EST = ATE_hat, VAR = Vhat, 
              CI = c(lower_ci, upper_ci),
              EIF = EIF_hat))
  
}

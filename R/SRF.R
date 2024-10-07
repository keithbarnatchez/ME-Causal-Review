# Stochastic intervention estimator
#' @param y length n vector of observations
#' @param x n x d matrix of covariates
#' @param a length n vector of exposures
#' @param delta exposure shift
#' @param sl.lib SuperLearner libraries for fitting nuisance models
#' @param ... Extra arguments for osqp solver
#' @export

srf <- function(y, x, a, delta, family = gaussian(),
                sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction"), ...) {
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a = a)
  
  # outcome model
  mumod <- SuperLearner(Y = y, X = xa, family = family, SL.library = sl.lib)
  muhat <- mumod$SL.predict
  mutilde <- c(predict(mumod, newdata = data.frame(x, a = a + delta))$pred)
  
  # GPS model
  pimod <- SuperLearner(Y = a, X = x, SL.library = "SL.glm")
  pimod.vals <- pimod$SL.predict
  pimod.sd <- sd(a - pimod.vals)
  
  a.std <- c(a - pimod.vals)/pimod.sd
  a.shift <- c((a - delta) - pimod.vals)/pimod.sd
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / pimod.sd # don't forget le Jacobian!
  pitilde <- approx(x = dens$x, y = dens$y, xout = a.shift)$y / pimod.sd # don't forget le Jacobian!
  
  # Point Estimates
  psi <- c(y - muhat)*(pitilde/pihat) + mutilde
  mu <- mean(psi, na.rm = TRUE)
  theta <- mean(psi - y, na.rm = TRUE)
  
  # Efficient IF
  mu_eif <- psi - mu
  theta_eif <- mu_eif - (y - mean(y))
  
  # variances (needs work)
  omega2 <- var(theta_eif, na.rm = TRUE)/n
  
  # Rely on asymptotics for CI construction
  lower_ci <- theta - sqrt(omega2)*qnorm(0.975)
  upper_ci <- theta + sqrt(omega2)*qnorm(0.975)
  
  return(list(EST = theta,VAR = omega2,
              CI = c(lower_ci, upper_ci), 
              EIF = theta_eif))
    
  
}

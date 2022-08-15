
csme <- function(y, astar, w, beta, weights){
  
  delta <- function(y, astar, eta, sigma_me, sigma_ep) {
    astar + sigma_me*y*eta/sigma_ep
  }
  

  
}

eefun_csme_aipw <- function(data, model) {
  Y <- data$Y
  Astar <- data$Astar
  L2 <- data$L2
  Lmat <- model.matrix(model, data = data)
  sw <- data$sw
  delta <- function(beta1, beta3, sigma_ep) {
    Astar + sigma_me*(beta1 + beta3*L2)*Y / sigma_ep
  }
  condexp <- function(beta0, beta1, beta2, beta3, sigma_ep) {
    (beta0 + (beta1 + beta3*L2)*delta(beta1, beta3, sigma_ep) +
       beta2*L2) /
      (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
  }
  condvar <- function(beta1, beta3, sigma_ep) {
    sigma_ep / (1 + ((beta1 + beta3*L2)^2)*sigma_me / sigma_ep)[[1]]
  }
  function(theta) {
    p  <- length(theta)
    p1 <- length(coef(model))
    rho <- Lmat %*% theta[1:p1]
    
    score_eqns <- apply(Lmat, 2, function(x) sum((Astar - rho) * x))
    
    c(score_eqns,
      Astar - theta[p-7],
      sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                      theta[p-2])),
      sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                      theta[p-2]))*L2,
      sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                      theta[p-2]))*
        delta(theta[p-5], theta[p-3], theta[p-2]),
      sw*(Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                      theta[p-2]))*
        L2*delta(theta[p-5], theta[p-3], theta[p-2]),
      sw*(theta[p-2] - theta[p-2]*
            (Y - condexp(theta[p-6], theta[p-5], theta[p-4], theta[p-3],
                         theta[p-2]))^2 /
            condvar(theta[p-5], theta[p-3], theta[p-2])),
      L2 - theta[p-1],
      theta[p-5] + theta[p-3]*theta[p-1] - theta[p]
    )
  }
}

failed <- TRUE
j <- 1

while(failed == TRUE & j < 6) {
  
  failed <- FALSE
  startvec <- c(coef(denom_mod), mean(Astar), coef(mod), sigma(mod)^2,
                mean(L2), guess)*(j == 1) +
    c(coef(denom_mod), mean(Astar), rnorm(4, 0, j/5), sigma(mod)^2,
      mean(L2), guess)*(j > 1)
  results_csme_aipw <-
    tryCatch(m_estimate(estFUN = eefun_csme_aipw, data = data,
                        outer_args = list(denom_mod),
                        root_control =
                          setup_root_control(start = startvec)),
             error = function(e) { failed <<- TRUE})
  if (failed == FALSE) {
    if (abs(coef(results_csme_aipw)[11]) > 3) { failed <- TRUE }
  }
  j <- j + 1
  
}

bias_aipw_ps <- coef(results_csme_aipw)[11] - true_effect
se_aipw_ps <- sqrt(vcov(results_csme_aipw)[11, 11])
coverage_aipw_ps <-
  1*(coef(results_csme_aipw)[11] - 1.96*se_aipw_ps < true_effect &
       coef(results_csme_aipw)[11] + 1.96*se_aipw_ps > true_effect)
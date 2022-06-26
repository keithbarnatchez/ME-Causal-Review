# ------------------------------------------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS
# ------------------------------------------------------------------------------
# This file contains functions for implementing the following ME correction
# methods: 1) propensity score calibration, 2) indirect SIMEX (Kyle et al. 2015),
# 3) IV, 4) 
#

# -------------------------------------------
# PROPENSITY SCORE CALIBRATION
# -------------------------------------------
psc_implement <- function(data) {
  #' Implements the Sturmer propensity score calibration approach
  #' INPUTS:
  #' - data: Analysis dataset obtained from the generate_data() function. 
  #'         
  #' - v_idx: Logical vector marking observations that are part of validation
  #'          dataset
  #' OUTPUTS:
  #' - Dataset containing IPTW weights obtained via PSC
  
  # Fit propensity score models in validation data
  val_data <- data[which(data$v_idx==1),]
  ep_ps_mod <- glm(A ~ W + Z, data=val_data,family='binomial') # error-prone
  gs_ps_mod <- glm(A ~ X + Z, data=val_data,family='binomial') # gold standard
  val_data$ep_ps <- predict(ep_ps_mod, type='response')
  val_data$gs_ps <- predict(gs_ps_mod, type='response')
  
  # Fit model relating gold-standard PS to error-prone PS
  ps_rel_model <- glm(gs_ps ~ ep_ps + A + Z, data=val_data, family=gaussian(link=),
                      link)
  
  # Fit propensity score model in main data so it can be projected to GS measure
  # via the model fitted with validation data
  data$ep_ps <- predict(glm(A ~ W + Z, data=data,family='binomial'),
                        type='response')
  data$gs_ps_hat <- predict(ps_rel_model,newdata=data,type='response')
  
  # Project X and W
  cal_mod <- lm(X ~ W, data=val_data)
  data$xcal <- predict(cal_mod,newdata=data)
  data$XX <- data$X ; data$X <- data$xcal
  data$ps2 <- predict(gs_ps_mod,newdata=data,type='response')
  data$X <- data$XX
  
  # Get weights for IPTW
  data$iptw <- ifelse(data$A==1,
                      1/data$gs_ps_hat,
                      1/(1-data$gs_ps_hat))
  data$iptw2 <- ifelse(data$A==1,
                      1/data$ps2,
                      1/(1-data$ps2))
  
  # thinking about returning altered df but for now returning ipw
  return(data$iptw)
}

psc <- function(data) {
  
  # Get weights
  data$iptw <- psc_implement(data)
  
  # Estimate ATE
  ATE <- lm(Y ~ A, weights = data$iptw, data=data)$coefficients[2]
  return(ATE)
}

# -------------------------------------------
#  SIMEX
# -------------------------------------------

simex_indirect_implement <- function(data) {
  #' Implements the indirect SIMEX adjustment described in Kyle et al. (2016)
  #' 
  
  # Estimate ME variance
  sig_u_hat <- sd(data$X[which(data$v_idx==1)] - data$W[which(data$v_idx==1)])
  
  # Run SIMEX on PS coefficient
  ps_model <- glm(A ~ W + Z, data=data, family='binomial', x=TRUE)
  simex_model <- simex(model=ps_model,
                       SIMEXvariable = "W",
                       measurement.error = sig_u_hat)
  
  # Get weights
  e_hat <- predict(simex_model,type='response')
  w_hat <- ifelse(data$A==1,
                  1/e_hat,
                  1/(1-e_hat))
  
  # Estimate ATE
  ATE_mod <- lm(Y ~ A, weights=data$iptw,data=data)$coefficients[2]
  
  return(ATE_mod)
}

simex_bootstrap <- function(data, nboot=1e3) {
    #' Obtains a bootstrap confidence interval for the SIMEX-adjusted ATE estimate
    #' Runs nboot instances of simex_indirect_implement, each with a new
    #' bootstrapped dataset
    #' INPUTS:
    #' - data: The dataset to be bootstrapped
    #' - nboot: Number of bootstrap iterations
    #' OUTPUTS:
    #' - A vector containing the lower and upper bound of a 95% CI
  
    # Current iteration
    ests <- rep(NA,nboot)
    for (b in 1:nboot) {
      boot_idx <- sample(1:nrow(data),replace=T)
      curr_mod <-  simex_indirect_implement(data[boot_idx,])
      ests[b] <- curr_mod$coefficients[2]
    }
    return(quantile(ests, probs=c(0.025,0.975)))
}

simex_indirect <- function(data, bs=T) {
  #' Implement the indirect SIMEX correction procedure described in Kyle et al.
  #' (2016)
  #'
  #' INPUTS:
  #' - data: A dataframe created with the generate_data() function
  #' - bs: Logical that = TRUE if want to bootstrap
  #' 
  #' OUTPUTS:
  #' - Estimate of ATE
  
  # Get weights for IPTW
  ATE_mod <- simex_indirect_implement(data)
  ATE_est <- ATE_mod$coefficients[2]
  
  if (bs==TRUE) {
    CI_est <- simex_bootstrap()
  }
}

iv_confounder <- function(data) {
  #' Implements instrumental variables correction 
  #'
  #'
  #'

  # Fit the IV model (first and second stage) with AER package
  iv_mod <- ivreg(Y ~ W + Z + A | Z + A + V, data=data)
  
  return(list(iv_mod$coefficients[4], # point estimate
         confint(iv_mod)[4,])) # confidence interval
}

ate_ideal <- function(data) {
  #' Computes ATE under ideal conditions (using X, correctly specified model)
  #' 
  #' Returns ATE estimate
  ps_model <- glm(A ~ X + Z, data=data, family='binomial')
  e_hat <- predict(ps_model,type='response')
  w_hat <- ifelse(data$A==1,
                  1/e_hat,
                  1/(1-e_hat))
  # Estimate ATE
  ATE_mod <- lm(Y ~ A, weights=w_hat,data=data)$coefficients[2]
  results <- list(ATE_mod$coefficients[2],
                  confint(ATE_mod)[2,])
  return(results)
}

ate_naive <- function(data) {
  #' Computes ATE when naively using W in place of X in estimating propensity
  #' scores
  #' 
  #' Returns ATE estimates
  ps_model <- glm(A ~ W + Z, data=data, family='binomial')
  e_hat <- predict(ps_model,type='response')
  w_hat <- ifelse(data$A==1,
                  1/e_hat,
                  1/(1-e_hat))
  # Estimate ATE
  ATE_mod <- lm(Y ~ A, weights=w_hat,data=data)
  results <- list(ATE_mod$coefficients[2],
                  confint(ATE_mod)[2,])
  
  return(results)
}

mime <- function(data,m=20) {
  #' Performs multiple imputation for ME correction, following Webb-Vargas et 
  #' al. (2015) with the one modification that we assume there is an internal,
  #' not external, validation sample. Makes use of the 'mice' package to perform
  #' multiple imputations
  #' 
  
  # Turn into missing data problem explicitly by recasting the unobserved X
  # values as missing so that we can use mice
  data_imp <- data %>% mutate(X = replace(X,v_idx==0,NA))
  
  # Run the (congenial) MI procedure
  imps <- mice(data_imp,method='norm.boot',m=m)
  
  # Loop through imputed datasets, estimate PS 
  vals <- sapply(1:m, function(d, imps, ...) {
    
    # fit propensity score model with d-th dataset
    curr_data <- complete(imps,d)
    ps_hat <- predict(glm(A ~ X + Z,family='binomial',data=curr_data),
                      type='response')
    w_hat <- ifelse(curr_data$A==1,
                    1/ps_hat,
                    1/(1-ps_hat))
    
    # Record ATE estimate and its SE
    ATE <- lm(Y ~ A, data=curr_data,weights=w_hat)
    beta <- ATE$coefficients[2] 
    se <- sqrt(diag(vcov(ATE))[2]) 
    
    return(c(beta = beta, se = se))
    
  }, imps = imps)
  
  betas <- vals[1,]
  ses <- vals[2,]
  
  return(mean(betas))
  
}

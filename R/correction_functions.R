# -------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS
# -------------------------------------------


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
  ep_ps_mod <- glm(T ~ W + Z, data=val_data,family='binomial') # error-prone
  gs_ps_mod <- glm(T ~ X + Z, data=val_data,family='binomial') # gold standard
  val_data$ep_ps <- predict(ep_ps_mod, type='response')
  val_data$gs_ps <- predict(gs_ps_mod, type='response')
  
  # Fit model relating gold-standard PS to error-prone PS
  ps_rel_model <- lm(gs_ps ~ ep_ps + T + Z, data=val_data)
  
  # Fit propensity score model in main data so it can be projected to GS measure
  # via the model fitted with validation data
  data$ep_ps <- predict(glm(T ~ W + Z, data=data,family='binomial'),
                        type='response')
  data$gs_ps_hat <- predict(ps_rel_model,newdata=data)
  
  # Project X and W
  cal_mod <- lm(X ~ W, data=val_data)
  data$xcal <- predict(cal_mod,newdata=data)
  data$XX <- data$X ; data$X <- data$xcal
  data$ps2 <- predict(gs_ps_mod,newdata=data,type='response')
  data$X <- data$XX
  
  # Get weights for IPTW
  data$iptw <- ifelse(data$T==1,
                      1/data$gs_ps_hat,
                      1/(1-data$gs_ps_hat))
  data$iptw2 <- ifelse(data$T==1,
                      1/data$ps2,
                      1/(1-data$ps2))
  
  # thinking about returning altered df but for now returning ipw
  return(data$iptw)
}

psc <- function(data) {
  
  # Get weights
  data$iptw <- psc_implement(data)
  
  # Estimate ATE
  ATE <- lm(Y ~ T, weights = data$iptw, data=data)$coefficients[2]
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
  ps_model <- glm(T ~ W + Z, data=data, family='binomial', x=TRUE)
  simex_model <- simex(model=ps_model,
                       SIMEXvariable = "W",
                       measurement.error = sig_u_hat)
  
  # Get weights
  e_hat <- predict(simex_model,type='response')
  w_hat <- ifelse(data$T==1,
                  1/e_hat,
                  1/(1-e_hat))
  
  return(w_hat)
}

simex_indirect <- function(data) {
  #' Implement the indirect SIMEX correction procedure described in Kyle et al.
  #' (2016)
  #'
  #' INPUTS:
  #' - data: A dataframe created with the generate_data() function
  #' 
  #' OUTPUTS:
  #' - Estimate of ATE
  
  # Get weights for IPTW
  data$iptw <- simex_indirect_implement(data)
  
  # Estimate ATE
  ATE <- lm(Y ~ T, weights=data$iptw,data=data)$coefficients[2]
  
}

iv_confounder <- function(data) {
  #' Implements instrumental variables correction 
  #'
  #'
  #'
  
  data$xhat <- lm(W ~ Z + T + V, data=data)$fitted.values
  stage2reg <- lm(Y ~ xhat + Z + T, data=data)
  
  return(stage2reg$coefficients[4])
}

ate_ideal <- function(data) {
  #' Computes ATE under ideal conditions (using X, correctly specified model)
  #' 
  #' Returns ATE estimate
  ps_model <- glm(T ~ X + Z, data=data, family='binomial')
  e_hat <- predict(ps_model,type='response')
  w_hat <- ifelse(data$T==1,
                  1/e_hat,
                  1/(1-e_hat))
  # Estimate ATE
  ATE <- lm(Y ~ T, weights=w_hat,data=data)$coefficients[2]
  
  return(ATE)

}

ate_naive <- function(data) {
  #' Computes ATE when naively using W in place of X in estimating propensity
  #' scores
  #' 
  #' Returns ATE estimates
  ps_model <- glm(T ~ W + Z, data=data, family='binomial')
  e_hat <- predict(ps_model,type='response')
  w_hat <- ifelse(data$T==1,
                  1/e_hat,
                  1/(1-e_hat))
  # Estimate ATE
  ATE <- lm(Y ~ T, weights=w_hat,data=data)$coefficients[2]
  
  return(ATE)
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
  betas <- rep(NA,m)
  ses   <- rep(NA,m)
  for (d in 1:m) {
    # fit propensity score model with d-th dataset
    curr_data <- complete(imps,d)
    ps_hat <- predict(glm(T ~ X + Z,family='binomial',data=curr_data),
                      type='response')
    w_hat <- ifelse(data$T==1,
                    1/ps_hat,
                    1/(1-ps_hat))
    
    # Record ATE estimate and its SE
    ATE <- lm(Y ~ T, data=curr_data,weights=w_hat)
    betas[d] <- ATE$coefficients[2] 
    ses[d] <- sqrt(diag(vcov(ATE))[2]) 
  }
  return(mean(betas))
}

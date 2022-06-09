# -------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS
# -------------------------------------------


# -------------------------------------------
# PROPENSITY SCORE CALIBRATION
# -------------------------------------------
psc <- function(data,v_idx) {
  #' Implements the Sturmer propensity score calibration approach
  #' INPUTS:
  #' - data: Analysis dataset obtained from the generate_data() function. 
  #'         
  #' - v_idx: Logical vector marking observations that are part of validation
  #'          dataset
  #' OUTPUTS:
  #' - Dataset containing IPTW weights obtained via PSC
  
  # Fit propensity score models in validation data
  val_data <- data[which(v_idx==1),]
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

# -------------------------------------------
#  SIMEX
# -------------------------------------------

simex_indirect <- function(data, v_idx) {
  #' Implements the indirect SIMEX adjustment described in Kyle et al. (2016)
  #' 
  
  # Estimate ME variance
  sig_u_hat <- sd(data$X[which(v_idx==1)] - data$W[which(v_idx==1)])
  
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







n <- 10000
v_share <- 0.1
v_idx <- rbinom(n,size=1,prob=v_share)

data <- generate_data(n,sig_u = 0.2)
data_iptw <- psc(data,v_idx)

simex_weights <- simex_indirect(data,v_idx)


lm(Y ~ T + Z + W, data=data_iptw)
lm(Y ~ T + Z + X, data=data_iptw)
lm(Y ~ T+0, data=data_iptw, weights=data_iptw$iptw)
lm(Y ~ T+0, data=data_iptw, weights=data_iptw$iptw2)
lm(Y ~ T+0, data=data_iptw, weights=simex_weights)



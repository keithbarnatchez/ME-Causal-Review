#-------

psc <- function(data,v_idx) {
  
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
  
  return(data)
}

n <- 10000
v_share <- 0.1
v_idx <- rbinom(n,size=1,prob=v_share)

data <- generate_data(n)
data_iptw <- psc(data,v_idx)


lm(Y ~ T + Z + W, data=data_iptw)
lm(Y ~ T, data=data_iptw, weights=data_iptw$iptw)
lm(Y ~ T, data=data_iptw, weights=data_iptw$iptw2)

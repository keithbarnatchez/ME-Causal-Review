# ------------------------------------------------------------------------------
# CAUSAL ME CORRECTION FUNCTIONS
# ------------------------------------------------------------------------------
# This file contains functions for implementing the following ME correction
# methods: 1) propensity score calibration, 2) indirect SIMEX (Kyle et al. 2015),
# 3) IV, 4) 
#
#
#
# ---------------------------------------------------------
# utility functions to supplement other calculations

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

simex_predictions <- function(simex_model) {
  
  # Get log odds
  lodds <- model.matrix(simex_model$model)%*%simex_model$coefficients
  
  return(expit(lodds)) # return predicted probabilities
}

aipw <- function(data,pihat=NULL,methods = 'SL.glm') {
  #' Main function for implementing AIPW estimator. Calls aipw_calc for the 
  #' actual calculation after calculating all necessary objects
  #' 
  #' INPUTS: 
  #' - data: A dataframe created by the generate_data() function
  #' - methods: Character vector with libraries to feed into SuperLearner. By
  #'            default, we just use regular glm (no non-parametric stuff)
  #'
  #' Outputs:
  #' - Estimate of the ATE
  
  # Outcome model
  y <- data$Y
  X <- data %>% select(W,Z,A)
  outcome_model <- SuperLearner(Y=y,X=X,SL.library = methods, family=gaussian)
  
  # Predict cond on A=0, A=1
  rhs0 <- X %>% mutate(A=0) # scenario where all tmt vals are 0
  rhs1 <- X %>% mutate(A=1) # scenario where all tmt vals are 1
  yhat0 <- predict(outcome_model,newdata=rhs0,X=X,Y=y)$pred # E[Y|X,A=0]
  yhat1 <- predict(outcome_model,newdata=rhs1,X=X,Y=y)$pred # E[Y|X,A=1]
  
  if (is.null(pihat)) { # if aren't pre-supplying propensity scores
    # Get predicted propensity scores
    Xtilde <- X %>% select(-A)
    tmt_model <- SuperLearner(Y=X$A, X=Xtilde, SL.library=methods,family=binomial())
    pihat <- predict(tmt_model, type='response')$pred
  }
  
  # Calculate the AIPW estimate and return the result
  return(  aipw_calc(y,data$A,pihat,yhat0,yhat1)  )
}

aipw_calc <- function(y,a,pihat,yhat0,yhat1) {
  #' Implements the augmented inv probability estimator (see eq3 of Glynn &
  #' Quinn, Political Analysis 2010)
  #' INPUTS: 
  #' - y: outcome
  #' - a: exposure
  #' - pihat: predicted propensity scores
  #' - yhat0: E(Y|A=0,Z)
  #' - yhat1: E(Y|A=1,Z)
  #'
  #' OUTPUTS:
  #' - Estimate of ATE
  
  # Compute relevant terms for estimator
  term1 <- a*y/pihat - (1-a)*y/(1-pihat)
  term2 <- (a-pihat)/(pihat*(1-pihat))
  term3 <- (1-pihat)*yhat1 + pihat*yhat0
  
  # ATE estimator
  ATE_hat <- mean(term1 - term2*term3)
  
  # Sandwich estimator
  Vhat <- (term1 - term2*term3) - ATE_hat
  Vhat <- sum(Vhat^2)/(length(Vhat)^2)
  
  # Rely on asymptotics for CI construction
  lower_ci <- ATE_hat - sqrt(Vhat)*qnorm(0.975)
  upper_ci <- ATE_hat + sqrt(Vhat)*qnorm(0.975)
  
  return(list(ATE=ATE_hat,
              CI=c(lower_ci,upper_ci)) )
  
}


# ---------------------------------------------------------
# ---------------------------------------------------------
#                PROPENSITY SCORE CALIBRATION
# ---------------------------------------------------------
psc_implement <- function(data) {
  #' Implements a modified Sturmer propensity score calibration approach, where
  #' the ATE is estimated via IPTW. As such, a beta regression model is assumed
  #' between the gold standard and naive propensity scores to ensure all predicted
  #' weights are in (0,1)
  #'  
  #' INPUTS:
  #' - data: Analysis dataset obtained from the generate_data() function
  #'         
  #' - v_idx: Logical vector marking observations that are part of validation
  #'          dataset
  #' OUTPUTS:
  #' - Dataset containing IPTW weights obtained via PSC
  
  # Fit propensity score models in validation data
  val_data <- data[which(data$v_idx==1),]
  ep_ps_mod <- glm(A ~ W + Z, data=val_data,family='binomial') # error-prone
  gs_ps_mod <- glm(A ~ X + Z, data=val_data,family='binomial') # gold standard
  val_data$ep_ps <- predict(ep_ps_mod, type='response') # predict ep pscores in val data
  val_data$gs_ps <- predict(gs_ps_mod, type='response') # predict gs pscores in val data
  
  # Fit model relating gold-standard PS to error-prone PS
  # ps_rel_model <- glm(gs_ps ~ ep_ps + A + Z, data=val_data, family=gaussian(link=),
  #                     link)
  # Note: using beta regressions at the moment to respect 0/1 bounds and to 
  # prevent issues with weights being negative 
  n_val <- nrow(val_data)
  val_data$gs_ps <- (val_data$gs_ps*(n_val-1) + 0.5)/n_val # smithson (2006) transformation to ensure (0,1) not [0,1]
  ps_rel_model <- betareg(gs_ps ~ ep_ps + A, data=val_data)
  
  # Fit propensity score model in main data so it can be projected to GS measure
  # via the model fitted with validation data
  data$ep_ps <- predict(glm(A ~ W + Z, data=data,family='binomial'),
                        type='response') # getting ep pscores in main data
  # Use the GS-EP pscore model to predict GS pscores in the main data
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
  # mean(data$Y * data$A * data$iptw) - mean(data$Y * (1-data$A) * data$iptw) 
  outcome_mod <- lm(Y ~ A, weights=data$iptw,data=data)
  return(ATE_mod)
}

psc_implement_reg <- function(data) {
  # Fit propensity score models in validation data
  val_data <- data[which(data$v_idx==1),]
  ep_ps_mod <- glm(A ~ W + Z, data=val_data,family='binomial') # error-prone
  gs_ps_mod <- glm(A ~ X + Z, data=val_data,family='binomial') # gold standard
  val_data$ep_ps <- predict(ep_ps_mod, type='response') # predict ep pscores in val data
  val_data$gs_ps <- predict(gs_ps_mod, type='response') # predict gs pscores in val data
  
  # Fit model relating gold-standard PS to error-prone PS
  # ps_rel_model <- glm(gs_ps ~ ep_ps + A + Z, data=val_data, family=gaussian(link=),
  #                     link)
  ps_rel_model <- lm(gs_ps ~ ep_ps + A, data=val_data)
  
  # Fit propensity score model in main data so it can be projected to GS measure
  # via the model fitted with validation data
  data$ep_ps <- predict(glm(A ~ W + Z, data=data,family='binomial'),
                        type='response') # getting ep pscores in main data
  # Use the GS-EP pscore model to predict GS pscores in the main data
  data$gs_ps_hat <- predict(ps_rel_model,newdata=data,type='response') 
  
  # fit the final model
  ATE_mod <- lm(Y ~ A + gs_ps_hat, data=data)
  return(list(ps_rel_model,ATE_mod))
}

psc_bootstrap <- function(data,nboot,iptw) {
  #' Giving a specified number of bootstrap iterations, generates a bootstrap
  #' confidence interval of the PSC ATE estimate
  #' INPUTS:
  #' - data: Analysis dataset obtained from the generate_data() function
  #' - nboot: Bootstrap iterations
  #' OUTPUTS:
  #' - A vector containing the CI lower and upper bounds

  # Current iteration
  ests <- rep(NA,nboot)
  for (b in 1:nboot) {
    boot_idx <- sample(1:nrow(data),replace=T)
    
    if (iptw==1) {
      curr_mod <-  psc_implement(data[boot_idx,])
      ests[b] <- curr_mod$coefficients[2]
    }
    else {
      curr_mods <-  psc_implement_reg(data[boot_idx,])
      B_A<- curr_mods[[2]]$coefficients['A'] ; B_X <- curr_mods[[2]]$coefficients['gs_ps_hat']
      L_A <- curr_mods[[1]]$coefficients['A'] ; L_X <- curr_mods[[1]]$coefficients['ep_ps']

      ests[b] <-  B_A - L_A*B_X/L_X
    }
  }
  return(quantile(ests, probs=c(0.025,0.975)))
}

psc <- function(data,nboot=100,iptw=1) {
  #' Outer function for the propensity score calibration method. Calls
  #' psc_implement to yield ATE estimate, and psc_bootstrap to obtain con
  #'
  #' INPUTS
  #'
  
  # Get weights
  if (iptw==1) {ATE_mod <- psc_implement(data)}
  else {psc_mods <- psc_implement_reg(data)}
  
  # Estimate ATE
  if (iptw==1) {
    ATE <- ATE_mod$coefficients['A']
  }
  else {
    B_A<- psc_mods[[2]]$coefficients['A']
    B_X <- psc_mods[[2]]$coefficients['gs_ps_hat']
    L_A <- psc_mods[[1]]$coefficients['A']
    L_X <- psc_mods[[1]]$coefficients['ep_ps']
    ATE <- B_A - L_A*B_X/L_X
  }
  
  # Bootstrap to obtain confidence interval
  
  CI_est <- psc_bootstrap(data,nboot,iptw)

  
  results <- list(ATE=ATE,
                  CI =CI_est)
  return(results)
}

# ---------------------------------------------------------
#                          SIMEX
# ---------------------------------------------------------

simex_indirect_implement <- function(data) {
  #' Implements the indirect SIMEX adjustment described in Kyle et al. (2016)
  #' Called within the simex_indirect() function
  #' 
  #' Returns a model object for the ATE
  
  # Estimate ME variance via the validation data (i.e. sd of X-W in val data)
  # Note: should be fine to do it this way since we're assuming non-dif/classical
  # but may need to explore extensions
  sig_u_hat <- sd(data$X[which(data$v_idx==1)] - data$W[which(data$v_idx==1)]) 
  
  # Run SIMEX on PS coefficient
  ps_model <- glm(A ~ W + Z, data=data, family='binomial', x=TRUE)
  simex_model <- simex(model=ps_model,
                       SIMEXvariable = "W",
                       measurement.error = sig_u_hat)
  
  # Get weights
  # e_hat <- predict(simex_model,type='response')
  e_hat <- simex_predictions(simex_model)
  w_hat <- ifelse(data$A==1,
                  1/e_hat,
                  1/(1-e_hat))
  data$w_hat <- w_hat
  
  # Estimate ATE
  ATE_mod <- lm(Y ~ A, weights=data$w_hat,data=data)
  
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

simex_indirect <- function(data, nboot=0) {
  #' Implement the indirect SIMEX correction procedure described in Kyle et al.
  #' (2016)
  #'
  #' INPUTS:
  #' - data: A dataframe created with the generate_data() function
  #' - bs: Logical that = TRUE if want to bootstrap
  #' 
  #' OUTPUTS:
  #' - Estimate of ATE
  #' REQUIRED PACKAGES:
  #' - simex
  
  # Get weights for IPTW
  ATE_mod <- simex_indirect_implement(data)
  ATE_est <- ATE_mod$coefficients[2]
  print('ATE est:')
  print(ATE_est)
  
  if (nboot>0) {
    CI_est <- simex_bootstrap(data,nboot)
  }
  
  return(list(ATE_mod$coefficients[2],
              confint(ATE_mod)))
}

# --------------------------------------------------------
#                 Instrumental variables
# --------------------------------------------------------

iv_confounder <- function(data) {
  #' Implements instrumental variables correction 
  #'
  #' INPUTS: 
  #' - data: A dataframe created with the generate_data() function
  #' 
  #' OUTPUTS:
  #' - A list containing the ATE estimate and confidence interval
  #' 
  #' REQUIRED PACKAGES:
  #' - AER

  # Fit the IV model (first and second stage) with AER package
  iv_mod <- ivreg(Y ~ W + Z + A | Z + A + V, data=data)
  
  return(list(iv_mod$coefficients[4], # point estimate
         confint(iv_mod)[4,])) # confidence interval
}

# --------------------------------------------------------
#            Naive and ideal approaches
# --------------------------------------------------------

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
  ATE_mod <- lm(Y ~ A, weights=w_hat,data=data)
  results <- list(ATE_mod$coefficients[2],
                  confint(ATE_mod)[2,])
  return(results)
}

ate_naive <- function(data) {
  #' Computes ATE when naively using W in place of X in estimating propensity
  #' scores, but with otherwise correctly-specified model
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

# --------------------------------------------------------
#               Multiple imputation
# --------------------------------------------------------

mime <- function(data,m=20) {
  #' Performs multiple imputation for ME correction, following Webb-Vargas et 
  #' al. (2015) with the one modification that we assume there is an internal,
  #' not external, validation sample. Makes use of the 'mice' package to perform
  #' multiple imputations
  #' 
  
  # Turn into missing data problem explicitly by recasting the unobserved X
  # values as missing so that we can use mice
  data_imp <- data %>% mutate(X = replace(X,v_idx==0,NA)) %>%
                       select(-v_idx)
  
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
  
  # grab ests and their SEs
  betas <- vals[1,]
  ses <- vals[2,]
  
  # compute SE est 
  ovr_mean <- mean(betas) 
  bw_var   <- sum( (betas - ovr_mean)^2 )/(length(betas)-1) # var bw ests
  wi_var   <- mean(ses) # within variance
  SE_est <- wi_var + ((1 + (1/length(betas))) * bw_var )
  CI <- c(-qnorm(.975)*SE_est + ovr_mean,
          qnorm(.975)*SE_est + ovr_mean)
  
  return(list(ATE=ovr_mean,
              CI=CI))
  
}
 
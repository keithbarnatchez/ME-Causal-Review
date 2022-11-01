# exposure_correction_functions.R
################################################################################

iv_exposure <- function(data) {
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
  iv_mod <- ivreg(Y ~ Z + X + W | Z + X + V, data=data)
  
  return(list(iv_mod$coefficients['W'], # point estimate
              confint(iv_mod)['W',])) # confidence interval
}

drf_ideal <- function(data,surface='linear') {
  #' Function for computing ERF under "ideal" scenario (we have measurements)
  #' of true exposure A
  #'
  if (surface=='linear') {
    erf_mod <- lm(Y ~ A + X + Z)
  }
  return(erf_mod)
}



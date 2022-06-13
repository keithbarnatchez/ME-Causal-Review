# --------------------------
# CAUSAL-ME OUTPUT FUNCTIONS
# --------------------------
#
# - percent bias
# - rmse
# - effect estimate
# - ATE variance
#
#
#


calc_stats <- function(methods,
                       data,
                       v_idx,
                       u,t,
                       nsim) {
  #' For a given iteration, calculates the following stats of interest:
  #' - bias, variance, CI coverage
  #' These stats are used to calculate the operating characteristics (across 
  #' all iterations) described in get_results.
  #' 
  #' INPUTS:
  #' - methods:
  #' - data:
  #' - v_idx: 
  #' - u,t:
  #' - nsim: Number of iterations
  #' 
  #' OUTPUTS:
  #' A dataframe with the following variables:
  #'     -
  #'     -
  #'     -
  #'     -
  #'     -
  #'     -
  
  # SIMEX
  if ('simex' %in% methods) {
    ATE <- simex_indirect(data,v_idx)
    
    # Get bias
    bias <- ATE-t
    temprow <- data.frame()
  }
  
  # PSC
  if ('psc' %in% methods) {
    
  }
  
  
}



get_results <- function(methods,
                        sig_u_grid,
                        bt_grid,
                        nsim=100) {
  
  #' Function for computing results of chosen ME-correction approaches, with
  #' varying data parameters
  #'
  #' INPUTS:
  #' - methods: Vector of strings for desired ME-correction methods to implement
  #'    - psc: Propensity score calibration (Sturmer 2005)
  #'    - simex_ind: Indirect Simex (Kyle et al. 2016)
  #'    - iv: Instrumental variables
  #'  - sig_u_grid: Vector of values for measurement error variance
  #'  
  #' OUTPUTS:
  #' - A dataframe
  #' 
  #' 
  
  for (u in sig_u_grid) { # loop over ME variances
  for (t in bt_grid) {
    # Initialize dataframe for each stat of interest
    stats_df <- data.frame(bias=double(),                          ,
                         method=character())
    for (s in 1:nsim) { # loop over simulation iterations
      
    }
  } 
} # sig_u_grid
} # bt_grid

# bias_plot <- function(results)
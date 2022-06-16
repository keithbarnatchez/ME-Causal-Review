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

calc_stats <- function(data,methods,sim_stats,t,s) {
  #' Calculates stats of interest for specified correction methods and updates 
  #' existing dataframe containing stats across simulation iterations
  #'
  #' INPUTS:
  #' - data:
  #' - methods:
  #' - simstats:
  #' - t:
  #' - s:
  #' OUTPUTS:
  #' - updated simstats dataframe
  #' 
  
  # Ideal and naive
  ATE_ideal <- ate_ideal(data) ; ATE_naive <- ate_naive(data)
  bias_ideal <- ATE_ideal-t ; bias_naive <- ATE_naive-t
  sim_stats <- rbind(sim_stats,
                     data.frame(bias=bias_ideal,ATE=ATE_ideal,method='Ideal',iteration=s))
  sim_stats <- rbind(sim_stats,
                    data.frame(bias=bias_naive,ATE=ATE_naive,method='Naive',iteration=s))
                                    
  # SIMEX
  if ('simex_ind' %in% methods) {
    # Get ATE estimate
    ATE <- simex_indirect(data)
    
    # Get bias
    bias <- ATE-t
    temprow <- data.frame(bias=bias,ATE=ATE,method='simex_ind',iteration=s)
    sim_stats <- rbind(sim_stats,temprow)
  }
  
  # PSC
  if ('psc' %in% methods) {
    # Get ATE estimate
    ATE <- psc(data)
    
    bias <- ATE-t
    temprow <- data.frame(bias=bias,ATE=ATE,method='psc',iteration=s)
    sim_stats <- rbind(sim_stats,temprow)
  }
  
  # IV 
  if ('iv' %in% methods) {
    # Get ATE estimate
    ATE <- iv_confounder(data)
    
    bias <- ATE-t
    temprow <- data.frame(bias=bias,ATE=ATE,method='iv',iteration=s)
    sim_stats <- rbind(sim_stats,temprow)
  }
  
  return(sim_stats)
}



get_stats_table <- function(methods,
                            u,t,n,
                            nsim) {
  #' For a given iteration, calculates the following stats of interest:
  #' - bias, variance, CI coverage
  #' These stats are used to calculate the operating characteristics (across 
  #' all iterations) described in get_results.
  #' 
  #' INPUTS:
  #' - methods:
  #' - v_idx: 
  #' - u,tn,:
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
  
  # set up dataframe for calculating operating characteristics by group
  sim_stats <- data.frame(bias=double(),ATE=double(),method=character(),it=double())
  for (s in 1:nsim) {
    
    # simulate data for current iteration
    data <- generate_data(n,sig_u=u,bt=t)
    
    # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
    sim_stats <- calc_stats(data,methods,sim_stats,t,s)
  } # for s in 1:nsim
  
  # Compute avg operating characteristics from stats of interest 
  final_stats <- sim_stats %>% group_by(method) %>%
    summarize(bias = mean(bias),
              mse = mean(bias^2),
              ATE = mean(ATE))
  final_stats$u=u ; final_stats$t=t ; final_stats$n=n
  return(final_stats)
} # calc_stats

get_results <- function(methods,
                        sig_u_grid,
                        bt_grid,
                        n_grid,
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
  
  # Initialize dataframe for each stat of interest
  op_chars <- data.frame(bias=double(),                          
                         method=character(),
                         mse=double(),
                         ATE=double())
  for (u in sig_u_grid) { # loop over ME variances
  for (t in bt_grid) { # loop over tmt effect sizes
  for (n in n_grid) {
  
    # Get operating characteristics for current grid point
    stats_df <- get_stats_table(methods,u,t,n,nsim)
    op_chars <- rbind(op_chars,stats_df)
    
  } # sig_u_grid
  } # bt_grid
  } # n_grid
  return(op_chars)
} # get_results 

bias_plot <- function(results) {
  
  results %>% ggplot(aes(x=bias))
  
}
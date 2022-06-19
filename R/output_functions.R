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

calc_stats <- function(data,methods,t,s) {
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
  stats <- data.frame(bias=bias_ideal,ATE=ATE_ideal,method='Ideal',iteration=s)
  stats <- rbind(stats, data.frame(bias=bias_naive,ATE=ATE_naive,method='Naive',iteration=s))
                                    
  # SIMEX
  if ('simex_ind' %in% methods) {
    # Get ATE estimate
    ATE <- simex_indirect(data)
    
    # Get bias
    bias <- ATE-t
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE,method='simex_ind',iteration=s))
  }
  
  # PSC
  if ('psc' %in% methods) {
    # Get ATE estimate
    ATE <- psc(data)
    bias <- ATE-t
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE,method='psc',iteration=s))
  }
  
  # IV 
  if ('iv' %in% methods) {
    # Get ATE estimate
    ATE <- iv_confounder(data)
    bias <- ATE-t
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE,method='iv',iteration=s))
  }
  
  return(stats)
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
  sim_stats_list <- lapply(1:nsim, function(s, n, u, t, ...) {
    
    # simulate data for current iteration
    data <- generate_data(n,sig_u=u,bt=t)
    
    # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
    return(calc_stats(data,methods,t,s))
    
  }, n = n, u = u, t = t) # for s in 1:nsim
  
  sim_stats <- do.call(rbind, sim_stats_list)
  
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
  scen_df <- expand.grid(sig_u = sig_u_grid, bt = bt_grid, n = n_grid,
                         KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, ...) {
    
    u <- scen$sig_u
    t <- scen$bt
    n <- scen$n
    
    # Get operating characteristics for current grid point
    op_tmp <- get_stats_table(methods,u,t,n,nsim)
    return(op_tmp)
    
  }, methods = methods, nsim = nsim)
    
  op_chars <- do.call(rbind, op_list)
    
  return(op_chars)
} # get_results 

bias_plot <- function(results) {
  
  results %>% ggplot(aes(x=bias))
  
}
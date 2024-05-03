# --------------------------
# CAUSAL-ME OUTPUT FUNCTIONS
# --------------------------
#
# - percent bias
# - rmse
# - effect estimate
# - EST variance

calc_stats <- function(data, methods, true_effect, sim) {
  
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
  EST_ideal <- srf_ideal(data) ; EST_naive <- srf_naive(data)
  CI_ideal <- EST_ideal[[3]] ; CI_naive <- EST_naive[[3]]
  
  # Bias
  bias_ideal <- EST_ideal[[1]] - true_effect ; bias_naive <- EST_naive[[1]] - true_effect
  
  # CI coverage
  ci_cov_ideal <- ifelse( (CI_ideal[1] <= true_effect) & (CI_ideal[2] >= true_effect) , 1, 0)
  ci_cov_naive <- ifelse( (CI_naive[1] <= true_effect) & (CI_naive[2] >= true_effect) , 1, 0)
  
  # Power
  pow_ideal <- ifelse(CI_ideal[1] > 0 | CI_ideal[2] < 0, 1, 0) 
  pow_naive <- ifelse(CI_naive[1] > 0 | CI_naive[2] < 0, 1, 0) 
  
  stats <- data.frame(bias = bias_ideal, EST = EST_ideal[[1]], 
                      ci_cov = ci_cov_ideal, pow = pow_ideal, 
                      method = 'Ideal', iteration = sim)
  stats <- rbind(stats, data.frame(bias = bias_naive, EST = EST_naive[[1]],
                                   ci_cov = ci_cov_naive, pow = pow_naive,
                                   method = 'Naive', iteration = sim))
                         
  # Propensity Score Calibration
  if ('psc' %in% methods) {
    
    EST <- srf_psc(data) 
    bias <- EST[[1]] - true_effect ; CI <- EST[[2]]
    ci_cov <- ifelse((CI[1] <= true_effect) & (true_effect <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, EST=EST[[1]], 
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'PSC', iteration = sim)) 
  
  }
  
  # Regression Calibration
  if ('rc' %in% methods) {
    
    EST <- srf_rc(data) 
    bias <- EST[[1]] - true_effect ; CI <- EST[[2]]
    ci_cov <- ifelse((CI[1] <= true_effect) & (true_effect <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, EST = EST[[1]], 
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'RC', iteration = sim)) 
    
  }
  
  # IV 
  if ('iv' %in% methods) {
    
    EST <- srf_iv(data)
    bias <- EST[[1]] - true_effect ; CI <- EST[[2]]
    ci_cov <- ifelse((CI[1] <= true_effect) & (true_effect <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias = bias, EST = EST[[1]],
                                     ci_cov = ci_cov, pow = pow, 
                                     method = 'IV', iteration = sim))
    
  }
  
  # SIMEX
  if ('simex' %in% methods) {
    
    EST <- srf_simex(data)
    bias <- EST[[1]] - true_effect ; CI <- EST[[2]]
    ci_cov <- ifelse((CI[1] <= true_effect) & (true_effect <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]] > 0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, EST = EST[[1]],
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'SIMEX', iteration = sim))
    
  }
  
  # MIME 
  if ('mime' %in% methods) {
    
    EST <- srf_mime(data)
    bias <- EST[[1]] - true_effect ; CI <- EST[[2]]
    ci_cov <- ifelse((CI[1] <= true_effect) & (true_effect <= CI[2]), 1, 0)
    pow <- ifelse(CI[[1]] > 0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, EST = EST[[1]],
                                     ci_cov = ci_cov, pow = pow, 
                                     method = 'MIME', iteration = sim))
    
  }
  
  # CV
  if ('cv' %in% methods) {
    
    EST <- srf_cv(data) 
    bias <- EST[[1]] - true_effect ; CI <- EST[[2]]
    ci_cov <- ifelse((CI[1] <= true_effect) & (true_effect <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, EST=EST[[1]], 
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'CV', iteration = sim)) 
    
  }
  
  return(stats)
  
}

get_results <- function(methods,
                        sig_u_grid, # me variance
                        ba_grid, 
                        aw_grid,
                        mis_grid, 
                        n_grid,     # n
                        rho_grid,   # correlation between confoundersß
                        nsim = 100,
                        mc.cores = 1) {
  
  #' Function for computing results of chosen ME-correction approaches, with
  #' varying data parameters
  #'
  #' INPUTS:
  #' - methods: Vector of strings for desired ME-correction methods to implement
  #'    - psc: Propensity score calibration (Sturmer 2005)
  #'    - simex_ind: Indirect Simex (Kyle et al. 2016)
  #'    - iv: Instrumental variables
  #'  - sig_u_grid: Vector of values for measurement error variance
  #'  - rho_grid: Vector of values specifying correlation between confounders
  #'  - nsim: number of simulations per grid point
  #'  
  #' OUTPUTS:
  #' - A dataframe containing calculated operating characteristics for 
  #'   each method and parameter combination (e.g. a given row will contain info
  #'   on the avg percent bias, mse, 95% CI coverage, etc for a specific correction
  #'   method and specific values of sig_u, ba, n, etc)
  
  # Initialize dataframe for each stat of interest
  scen_df <- expand.grid(sig_u = sig_u_grid, 
                         ba = ba_grid,
                         aw = aw_grid,
                         mis = mis_grid,
                         n = n_grid,
                         rho = rho_grid,
                         KEEP.OUT.ATTRS = TRUE, 
                         stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, mc.cores) {
    
    sig_u <- scen$sig_u # me variance
    ba <- scen$ba # effect of exposure on outcome
    aw <- scen$aw # effect of EP confounder on A
    mis <- scen$mis # effect of EP confounder on Y
    n <- scen$n # sample size
    rho <- scen$rho # rhoary outcome indicator
    
    # Keep track of progress
    print(paste(scen))
    
    # Get operating characteristics for current grid point
    sim_stats_list <- mclapply(1:nsim, function(sim, n, sig_u, aw, ba, mis, rho, methods) {
      
      # simulate data for current iteration
      data <- generate_data(n = n, sig_u = sig_u, rho = rho, aw = aw, mis = mis, ba = ba)
      
      # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
      return(calc_stats(data, methods = methods, true_effect = ba, sim = sim))
      
    }, n = n, sig_u = sig_u, rho = rho, aw = aw, ba = ba, mis = mis, methods = methods, mc.cores = mc.cores) # for s in 1:nsim
    
    sim_stats <- do.call(rbind, sim_stats_list)
    
    # Compute avg operating characteristics from stats of interest 
    final_stats <- sim_stats %>% group_by(method) %>%
      summarize(bias = 100*mean(bias/abs(ba), na.rm = T),
                rmse = sqrt(mean(bias^2, na.rm = T)),
                estimate = mean(EST, na.rm = T),
                ci_cov = 100*mean(ci_cov, na.rm = T),
                power = 100*mean(pow, na.rm = T))
    
    # make note of the grid point
    final_stats$n = n; final_stats$sig_u = sig_u; 
    final_stats$rho = rho; final_stats$aw = aw; 
    final_stats$ba = ba; final_stats$mis = mis;

    return(final_stats)
    
  }, methods = methods, nsim = nsim, mc.cores = mc.cores)
  
  # Save the op_chars table after ach grid point is done so we don't need to
  # wait for the entire run to be done (can plot grid points as they come in)  
  op_chars <- do.call(rbind, op_list)
    
  return(op_chars)
  
} # get_results 


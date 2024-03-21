# --------------------------
# CAUSAL-ME OUTPUT FUNCTIONS
# --------------------------
#
# - percent bias
# - rmse
# - effect estimate
# - ATE variance

calc_stats <- function(data, methods, a, s) {
  
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
  CI_ideal <- ATE_ideal[[3]] ; CI_naive <- ATE_naive[[3]]
  
  # Bias
  bias_ideal <- ATE_ideal[[1]] - a ; bias_naive <- ATE_naive[[1]] - a
  
  # CI coverage
  ci_cov_ideal <- ifelse( (CI_ideal[1] <= a) & (CI_ideal[2] >= a) , 1, 0)
  ci_cov_naive <- ifelse( (CI_naive[1] <= a) & (CI_naive[2] >= a) , 1, 0)
  
  # Power
  pow_ideal <- ifelse(CI_ideal[1] > 0 | CI_ideal[2] < 0, 1, 0) 
  pow_naive <- ifelse(CI_naive[1] > 0 | CI_naive[2] < 0, 1, 0) 
  
  stats <- data.frame(bias = bias_ideal, ATE = ATE_ideal[[1]], 
                      ci_cov = ci_cov_ideal, pow = pow_ideal, 
                      method = 'Ideal', iteration = s)
  stats <- rbind(stats, data.frame(bias = bias_naive, ATE = ATE_naive[[1]],
                                   ci_cov = ci_cov_naive, pow = pow_naive,
                                   method = 'Naive', iteration = s))
                         
  # Propensity Score Calibration
  if ('psc' %in% methods) {
    
    ATE <- ate_psc(data) 
    bias <- ATE[[1]] - a ; CI <- ATE[[2]]
    ci_cov <- ifelse((CI[1] <= a) & (a <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, ATE=ATE[[1]], 
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'PSC', iteration = s)) 
  
  }
  
  # Regression Calibration
  if ('rc' %in% methods) {
    
    ATE <- ate_rc(data, method = "aipw") 
    bias <- ATE[[1]] - a ; CI <- ATE[[2]]
    ci_cov <- ifelse((CI[1] <= a) & (a <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, ATE = ATE[[1]], 
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'RC', iteration = s)) 
    
  }
  
  # IV 
  if ('iv' %in% methods) {
    
    ATE <- ate_iv(data)
    bias <- ATE[[1]] - a ; CI <- ATE[[2]]
    ci_cov <- ifelse((CI[1] <= a) & (a <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias = bias, ATE = ATE[[1]],
                                     ci_cov = ci_cov, pow = pow, 
                                     method = 'IV', iteration = s))
    
  }
  
  # SIMEX
  if ('simex' %in% methods) {
    
    ATE <- ate_simex(data, method = "aipw")
    bias <- ATE[[1]] - a ; CI <- ATE[[2]]
    ci_cov <- ifelse((CI[1] <= a) & (a <= CI[2]), 1 , 0)
    pow <- ifelse(CI[[1]] > 0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, ATE = ATE[[1]],
                                     ci_cov = ci_cov, pow = pow,
                                     method = 'SIMEX', iteration = s))
    
  }
  
  # MIME 
  if ('mime' %in% methods) {
    
    ATE <- ate_mime(data, method = "aipw")
    bias <- ATE[[1]] - a ; CI <- ATE[[2]]
    ci_cov <- ifelse((CI[1] <= a) & (a <= CI[2]), 1, 0)
    pow <- ifelse(CI[[1]] > 0 | CI[[2]] < 0, 1, 0) # power
    stats <- rbind(stats, data.frame(bias = bias, ATE = ATE[[1]],
                                     ci_cov = ci_cov, pow = pow, 
                                     method = 'MIME', iteration = s))
    
  }
  
  return(stats)
  
}

get_results <- function(methods,
                        sig_u_grid, # me variance
                        aw_grid,    # trt effect
                        ba_grid,    # trt effect
                        bw_grid, 
                        n_grid,     # n
                        bin_grid,   # whether outcome is binary or not
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
  #'  - bin_grid: Vector of values specifying whether to use binary outcome. At 
  #'              most length 2 (either 0, 1 or (0,1) )
  #'  - nsim: number of simulations per grid point
  #'  
  #' OUTPUTS:
  #' - A dataframe containing calculated operating characteristics for 
  #'   each method and parameter combination (e.g. a given row will contain info
  #'   on the avg percent bias, mse, 95% CI coverage, etc for a specific correction
  #'   method and specific values of sig_u, ba, n, etc)
  
  # Initialize dataframe for each stat of interest
  scen_df <- expand.grid(sig_u = sig_u_grid, 
                         aw = aw_grid,
                         ba = ba_grid,
                         bw = bw_grid,
                         n = n_grid,
                         bin = bin_grid,
                         KEEP.OUT.ATTRS = TRUE, 
                         stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, mc.cores) {
    
    sig_u <- scen$sig_u # me variance
    aw <- scen$aw # effect of EP confounder on A
    ba <- scen$ba # effect of exposure on outcome
    bw <- scen$bw # effect of EP confounder on Y
    n <- scen$n # sample size
    bin <- scen$bin # binary outcome indicator
    
    # Keep track of progress
    print(paste(scen))
    
    # Get operating characteristics for current grid point
    sim_stats_list <- mclapply(1:nsim, function(s, n, sig_u, aw, ba, bw, bin) {
      
      # simulate data for current iteration
      data <- generate_data(n = n, sig_u = sig_u, binary = bin, aw = aw, ba = ba, bw = bw)
      
      # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
      return(calc_stats(data, methods, a = ba, s))
      
    }, n = n, sig_u = sig_u, bin = bin, aw = aw, ba = ba, bw = bw, mc.cores = mc.cores) # for s in 1:nsim
    
    sim_stats <- do.call(rbind, sim_stats_list)
    
    # Compute avg operating characteristics from stats of interest 
    final_stats <- sim_stats %>% group_by(method) %>%
      summarize(bias = 100*mean(bias/abs(ba), na.rm = T),
                rmse = sqrt(mean(bias^2, na.rm = T)),
                estimate = mean(ATE, na.rm = T),
                ci_cov = 100*mean(ci_cov, na.rm = T),
                power = 100*mean(pow, na.rm = T))
    
    # make note of the grid point
    final_stats$n = n; final_stats$sig_u = sig_u ; 
    final_stats$bin = bin; final_stats$aw = aw; 
    final_stats$ba = ba; final_stats$bw = bw

    return(final_stats)
    
  }, methods = methods, nsim = nsim, mc.cores = mc.cores)
  
  # Save the op_chars table after ach grid point is done so we don't need to
  # wait for the entire run to be done (can plot grid points as they come in)  
  op_chars <- do.call(rbind, op_list)
    
  return(op_chars)
  
} # get_results 


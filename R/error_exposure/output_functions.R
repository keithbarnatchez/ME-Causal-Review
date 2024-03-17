# ----------------------------------------------------
# CAUSAL-ME OUTPUT FUNCTIONS for exposure ME
# ----------------------------------------------------
#
# - percent bias
# - rmse
# - effect estimate
# - ERF variance


stats_getter <- function(results, true_effect) {
  
  #' Called within calc_stats() on each selected method
  #' Retrieves list of results of interest
  
  bias <- (results$EST - true_effect)/true_effect
  rmse <- sqrt(mean(bias^2))
  ci_cov <- ifelse((results$CI[1] <= true_effect) & (results$CI[2] >= true_effect), 1, 0)
  power <- ifelse((results$CI[1] <= 0) & (results$CI[2] >= 0), 0, 1)
  
  return(list(bias=bias, rmse = rmse, ci_cov = ci_cov, 
              pow = power, estimate = results$EST))
  
}

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
  
  sl.lib <- c("SL.mean", "SL.glm")
  stats <- c()
  
  # Extract variables to feed into correction functions
  Y <- data$Y ; W <- data$W ; X <- data$X ; V <- data$V
  A <- data$A ; A.star <- data$A.star ; v_idx <- data$v_idx
    
  # Ideal
  res_ideal <- erf_ideal(data)
  ideal_stats <- stats_getter(res_ideal, true_effect = a)
  stats <- data.frame(bias=ideal_stats$bias,
                      est = ideal_stats$estimate,
                      ci_cov = ideal_stats$ci_cov,
                      pow = ideal_stats$pow,
                      rmse = ideal_stats$rmse,
                      method = 'Ideal', iteration = s)
  
  # Naive
  res_naive <- erf_naive(data)
  naive_stats <- stats_getter(res_naive, true_effect = a)
  stats <- rbind(stats, data.frame(bias = naive_stats$bias,
                                   est = naive_stats$estimate,
                                   ci_cov = naive_stats$ci_cov,
                                   pow = naive_stats$pow,
                                   rmse = naive_stats$rmse,
                                   method = 'Naive', iteration = s))

  # Regression Calibration
  if ('rc' %in% methods) {

    res_rcal <- erf_rc(data)
    rcal_stats <- stats_getter(res_rcal, true_effect = a)
    stats <- rbind(stats, data.frame(bias = rcal_stats$bias,
                                     est = rcal_stats$estimate,
                                     ci_cov = rcal_stats$ci_cov,
                                     pow = rcal_stats$pow,
                                     rmse = rcal_stats$rmse,
                                     method = 'RC', iteration = s))
    
  }
  
  # CSME
  if ('csme' %in% methods) {
    
    res_csme <- erf_csme(data)
    csme_stats <- stats_getter(res_csme, true_effect = a)
    stats <- rbind(stats, data.frame(bias = csme_stats$bias,
                                     est = csme_stats$estimate,
                                     ci_cov = csme_stats$ci_cov,
                                     pow = csme_stats$pow,
                                     rmse = csme_stats$rmse,
                                     method = 'CSME',iteration=s))
    
  }
  
  # IV
  if ('iv' %in% methods) {
    
    res_iv <- erf_iv(data)
    iv_stats <- stats_getter(res_iv, true_effect = a)
    stats <- rbind(stats, data.frame(bias = iv_stats$bias,
                                     est = iv_stats$estimate,
                                     ci_cov = iv_stats$ci_cov,
                                     pow = iv_stats$pow,
                                     rmse = iv_stats$rmse,
                                     method = 'IV', iteration = s))
    
  }
  
  # SIMEX
  if ('simex' %in% methods) {
    
    res_simex <- erf_simex(data)
    simex_stats <- stats_getter(res_simex, true_effect = a)
    stats <- rbind(stats, data.frame(bias = simex_stats$bias,
                                     est = simex_stats$estimate,
                                     ci_cov = simex_stats$ci_cov,
                                     pow = simex_stats$pow,
                                     rmse = simex_stats$rmse,
                                     method = 'SIMEX', iteration = s))

  }
  
  # MIME
  if ('mime' %in% methods) {

    res_mime <- erf_mime(data)
    mime_stats <- stats_getter(res_mime, true_effect = a)
    stats <- rbind(stats, data.frame(bias = mime_stats$bias,
                                     est = mime_stats$estimate,
                                     ci_cov = mime_stats$ci_cov,
                                     pow = mime_stats$pow,
                                     rmse = mime_stats$rmse,
                                     method = 'MIME', iteration = s))
    
  }
  
  return(stats)
  
}

get_results <- function(methods,
                        sig_u_grid, # me variance
                        ba_grid,    # main tmt effect
                        n_grid,     # n
                        bin_grid,   # interactions
                        nsim = 100,
                        mc.cores = 1) {
  
  #' Function for computing results of chosen ME-correction approaches, with
  #' varying data parameters
  #'
  #' INPUTS:
  #'  - methods: Vector of strings for desired ME-correction methods to implement
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
                         ba = ba_grid, 
                         n = n_grid, 
                         bin = bin_grid,
                         KEEP.OUT.ATTRS = TRUE, 
                         stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, mc.cores) {
    
    sig_u <- scen$sig_u # me variance
    ba <- scen$ba # effect of exposure on outcome
    n <- scen$n # sample size
    bin <- scen$bin # binary outcome indicator
    
    # set up dataframe for calculating operating characteristics by group
    sim_stats_list <- mclapply(1:nsim, function(s, n, sig_u, ba, bin, methods, ...) {
      
      # Keep track of progress
      print(paste('On iteration', s, 'of', nsim))
      
      # simulate data for current iteration
      data <- generate_data(n = n, sig_u = sig_u, binary = bin, ba = ba)
      
      # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
      return(calc_stats(data = data, methods = methods, a = ba, s = s))
      
    }, n = n, sig_u = sig_u, ba = ba, bin = bin, methods = methods, mc.cores = mc.cores) # for s in 1:nsim
    
    sim_stats <- do.call(rbind, sim_stats_list)
    
    # Compute avg operating characteristics from stats of interest 
    final_stats <- sim_stats %>% group_by(method) %>%
      summarize(bias = 100*mean(bias/a, na.rm = T),
                mse = mean(bias^2, na.rm = T),
                est = mean(est, na.rm = T),
                ci_cov = 100*mean(ci_cov, na.rm = T),
                power = 100*mean(pow, na.rm = T))
    
    # make note of the grid point
    final_stats$n = n; final_stats$sig_u = sig_u ; 
    final_stats$bin = bin; final_stats$ba = ba; 
    
    return(final_stats)
    
  }, methods = methods, nsim = nsim, mc.cores = mc.cores)
  
  # Save the op_chars table after each grid point is done so we don't need to
  # wait for the entire run to be done (can plot grid points as they come in)  
  op_chars <- do.call(rbind, op_list)
  
  return(op_chars)
  
} # get_results 



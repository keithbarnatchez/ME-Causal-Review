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
  
  bias <- results$EST - true_effect
  ci_cov <- ifelse((results$CI[1] <= true_effect) & (results$CI[2] >= true_effect), 1, 0)
  power <- ifelse((results$CI[1] <= 0) & (results$CI[2] >= 0), 0, 1)
  
  return(list(bias = bias, ci_cov = ci_cov, 
              pow = power, est = results$EST))
  
}

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
  
  # Extract variables to feed into correction functions
  Y <- data$Y ; W <- data$W ; X <- data$X ; V <- data$V
  A <- data$A ; A.star <- data$A.star ; v_idx <- data$v_idx
    
  # Ideal
  res_ideal <- srf_ideal(data)
  ideal_stats <- stats_getter(res_ideal, true_effect = true_effect)
  stats <- data.frame(bias = ideal_stats$bias,
                      est = ideal_stats$est,
                      ci_cov = ideal_stats$ci_cov,
                      pow = ideal_stats$pow,
                      method = 'Ideal',
                      iteration = sim)
  
  # Naive
  res_naive <- srf_naive(data)
  naive_stats <- stats_getter(res_naive, true_effect = true_effect)
  stats <- rbind(stats, data.frame(bias = naive_stats$bias,
                                   est = naive_stats$est,
                                   ci_cov = naive_stats$ci_cov,
                                   pow = naive_stats$pow,
                                   method = 'Naive',
                                   iteration = sim))

  # Regression Calibration
  if ('rc' %in% methods) {

    res_rcal <- srf_rc(data)
    rcal_stats <- stats_getter(res_rcal, true_effect = true_effect)
    stats <- rbind(stats, data.frame(bias = rcal_stats$bias,
                                     est = rcal_stats$est,
                                     ci_cov = rcal_stats$ci_cov,
                                     pow = rcal_stats$pow,
                                     method = 'RC', 
                                     iteration = sim))
    
  }
  
  # CSME
  if ('csme' %in% methods) {
    
    res_csme <- srf_csme(data)
    csme_stats <- stats_getter(res_csme, true_effect = true_effect)
    stats <- rbind(stats, data.frame(bias = csme_stats$bias,
                                     est = csme_stats$est,
                                     ci_cov = csme_stats$ci_cov,
                                     pow = csme_stats$pow,
                                     method = 'CSME',
                                     iteration = sim))
    
  }
  
  # IV
  if ('iv' %in% methods) {
    
    res_iv <- srf_iv(data)
    iv_stats <- stats_getter(res_iv, true_effect = true_effect)
    stats <- rbind(stats, data.frame(bias = iv_stats$bias,
                                     est = iv_stats$est,
                                     ci_cov = iv_stats$ci_cov,
                                     pow = iv_stats$pow,
                                     method = 'IV', iteration = sim))
    
  }
  
  # SIMEX
  if ('simex' %in% methods) {
    
    res_simex <- srf_simex(data)
    simex_stats <- stats_getter(res_simex, true_effect = true_effect)
    stats <- rbind(stats, data.frame(bias = simex_stats$bias,
                                     est = simex_stats$est,
                                     ci_cov = simex_stats$ci_cov,
                                     pow = simex_stats$pow,
                                     method = 'SIMEX', 
                                     iteration = sim))

  }
  
  # MIME
  if ('mime' %in% methods) {

    res_mime <- srf_mime(data)
    mime_stats <- stats_getter(res_mime, true_effect = true_effect)
    stats <- rbind(stats, data.frame(bias = mime_stats$bias,
                                     est = mime_stats$est,
                                     ci_cov = mime_stats$ci_cov,
                                     pow = mime_stats$pow,
                                     method = 'MIME', 
                                     iteration = sim))
    
  }
  
  # CV
  if ('cv' %in% methods) {
    
    res_cv <- srf_cv(data)
    cv_stats <- stats_getter(res_cv, true_effect = true_effect)
    stats <- rbind(stats, data.frame(bias = cv_stats$bias,
                                     est = cv_stats$est,
                                     ci_cov = cv_stats$ci_cov,
                                     pow = cv_stats$pow,
                                     method = 'CV', 
                                     iteration = sim))
    
  }
  
  return(stats)
  
}

get_results <- function(methods,
                        sig_u_grid, # me variance
                        ba_grid,    # main tmt effect
                        aw_grid,    # effect of covariate on a
                        bw_grid,
                        mis_grid,
                        n_grid,   # sample size
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
                         aw = aw_grid, 
                         bw = bw_grid,
                         mis = mis_grid,
                         n = n_grid,
                         KEEP.OUT.ATTRS = TRUE, 
                         stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, mc.cores) {
    
    sig_u <- scen$sig_u # me variance
    ba <- scen$ba # effect of exposure on outcome
    aw <- scen$aw # effect size of confounder in gps
    bw <- scen$aw # effect size of confounder in om
    mis <- scen$mis # effect of EP confounder on Y
    n <- scen$n # sample size

    # Keep track of progress
    print(paste(scen))
    
    # set up dataframe for calculating operating characteristics by group
    sim_stats_list <- mclapply(1:nsim, function(sim, n, sig_u, ba, aw, bw, mis, methods, ...) {

      # simulate data for current iteration
      data <- generate_data(n = n, sig_u = sig_u, aw = aw, bw = bw, ba = ba, mis = mis)
      
      # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
      return(calc_stats(data = data, methods = methods, true_effect = ba, sim = sim))
      
    }, n = n, sig_u = sig_u, ba = ba, aw = aw, bw = bw, mis = mis, methods = methods, mc.cores = mc.cores) # for s in 1:nsim
    
    sim_stats <- do.call(rbind, sim_stats_list)
    
    # Compute avg operating characteristics from stats of interest 
    final_stats <- sim_stats %>% group_by(method) %>%
      summarize(bias = 100*mean(bias/abs(ba), na.rm = T),
                rmse = sqrt(mean(bias^2, na.rm = T)),
                estimate = mean(est, na.rm = T),
                ci_cov = 100*mean(ci_cov, na.rm = T),
                power = 100*mean(pow, na.rm = T))
    
    # make note of the grid point
    final_stats$n = n; final_stats$sig_u = sig_u; final_stats$ba = ba;  
    final_stats$aw = aw; final_stats$bw = bw; final_stats$mis = mis;
    
    return(final_stats)
    
  }, methods = methods, nsim = nsim, mc.cores = mc.cores)
  
  # Save the op_chars table after each grid point is done so we don't need to
  # wait for the entire run to be done (can plot grid points as they come in)  
  op_chars <- do.call(rbind, op_list)
  
  return(op_chars)
  
} # get_results 



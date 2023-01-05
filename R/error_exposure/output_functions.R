# ----------------------------------------------------
# CAUSAL-ME OUTPUT FUNCTIONS for exposure ME
# ----------------------------------------------------
#
# - percent bias
# - rmse
# - effect estimate
# - ATE variance
#
#
#

stats_getter <- function(results,true_effect) {
  #' Called within calc_stats() on each selected method
  #' Retrieves list of results of interest
  #'
  bias <- (results$estimate - true_effect)/true_effect
  rmse <- sqrt(mean(bias^2))
  ci_cov <- ifelse((results$CI[1]<=true_effect) & (results$CI[2]>=true_effect), 1, 0)
  power <- ifelse((results$CI[1]<=0) & (results$CI[2]>=0), 0, 1)
  
  return(list(
    bias=bias,rmse=rmse,ci_cov=ci_cov,pow=power,estimate=results$estimate
  ))
  
}

calc_stats <- function(data,methods,a,bax1,bax2,s,
                       sl.lib,
                       outcome_family=gaussian()) {
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
  
  sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
  true_effect <- a + 0.5*bax1 + 1*bax2
  stats <- c()
  
  # Extract variables to feed into correction functions
  Y <- data$Y ; X1 <- data$X1 ; X2 <- data$X2 ; A <- data$A ; Z <- data$Z
  v_idx <- data$v_idx ; A.star <- data$Astar
    
  # Ideal
  res_ideal <- ideal_method(Y,A,X1,X2,
                            sl.lib,
                            outcome_family)
  ideal_stats <- stats_getter(res_ideal,true_effect)
  
  stats <- rbind(stats, data.frame(bias=ideal_stats$bias,
                                   ATE=ideal_stats$estimate,
                                   ci_cov=ideal_stats$ci_cov,
                                   pow=ideal_stats$pow,
                                   rmse=ideal_stats$rmse,
                                   method='Ideal',iteration=s))
  
  # Naive
  res_ideal <- ideal_method(Y,A.star,X1,X2,
                            sl.lib,
                            outcome_family)
  naive_stats <- stats_getter(res_ideal,true_effect)
  stats <- rbind(stats, data.frame(bias=naive_stats$bias,
                                   ATE=naive_stats$estimate,
                                   ci_cov=naive_stats$ci_cov,
                                   pow=naive_stats$pow,
                                   rmse=naive_stats$rmse,
                                   method='Naive',iteration=s))

  # Regression calibration
  if ('rcal' %in% methods) {

    res_rcal <- rcal(Y,A,X1,X2,A.star,v_idx,sl.lib)
    
    # Get stats
    rcal_stats <- stats_getter(res_rcal,true_effect)
    stats <- rbind(stats, data.frame(bias=rcal_stats$bias,
                                     ATE=rcal_stats$estimate,
                                     ci_cov=rcal_stats$ci_cov,
                                     pow=rcal_stats$pow,
                                     rmse=rcal_stats$rmse,
                                     method='Reg. cal.',iteration=s))
  }
  
  # CSME
  if ('csme' %in% methods) {
    sig2_me <- var(data$Astar[v_idx==1] - data$A[v_idx==1])
    res_csme <- csme_linear(Y,A,X1,X2,A.star,sig2_me=sig2_me)
    
    # Get stats
    csme_stats <- stats_getter(res_csme,true_effect)
    stats <- rbind(stats, data.frame(bias=csme_stats$bias,
                                     ATE=csme_stats$estimate,
                                     ci_cov=csme_stats$ci_cov,
                                     pow=csme_stats$pow,
                                     rmse=csme_stats$rmse,
                                     method='CSME',iteration=s))
  }
  
  # SIMEX
  if ('simex' %in% methods) {
    # Estimate ME variance to plug in
    sigUhat <- var(A[v_idx==1] - A.star[v_idx==1])
    
    # Run SIMEX
    res_simex <- simex_direct(A.star,Y,cbind(X1,X2),
                              a0=0,a1=1,
                              degree=2, 
                              tau2=sigUhat)
    
    simex_stats <- stats_getter(res_simex,true_effect)
    stats <- rbind(stats, data.frame(bias=simex_stats$bias,
                                     ATE=simex_stats$estimate,
                                     ci_cov=simex_stats$ci_cov,
                                     pow=simex_stats$pow,
                                     rmse=simex_stats$rmse,
                                     method='SIMEX',iteration=s))

  }
  

  return(stats)
}



get_stats_table <- function(methods,
                            u,a,bax1,bax2,n,
                            nsim) {
  #' For a given iteration, calculates the following stats of interest:
  #' - bias, MSE, CI coverage
  #' These stats are used to calculate the operating characteristics (across 
  #' all iterations) described in get_results.
  #' 
  #' INPUTS:
  #' - methods:
  #' - v_idx: 
  #' - u,a,n,b:
  #' - nsim: Number of iterations
  #' 
  #' OUTPUTS:
  #' A dataframe with the following variables:
  #'     - method: method used (psc, iv, etc)
  #'     - bias: average percent bias
  #'     - mse: mean squared error
  #'     - ATE: average ATE estimate
  #'     - ci_cov: 95% CI coverage (share of CIs containing true ATE)
  #'     - u: ME variance
  #'     - a: true ATE
  #'     - n: sample size
  #'     - b: whether outcome is binary or not
  #'     
  
  # set up dataframe for calculating operating characteristics by group
  sim_stats_list <- lapply(1:nsim, function(s, n, u, a, bax1, bax2, ...) {
    
    # Keep track of progress
    print(paste('On iteration',s,'of',nsim))
    
    # simulate data for current iteration
    data <- gen_data(n,vshare=0.1,
                     b0=0,bA=1,bX1=0.7,bX2=-0.7, bA_X1=bax1, bA_X2=bax2,
                     muX1=0.5, muX2=1, muZ=1,
                     a0=2,aX1=0.9,aX2=-0.6,aZ=0.5,
                     sigU=u, sigA=1,sigE=1)
    
    # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
    return(calc_stats(data,methods,a,bax1,bax2,s))
    
  }, n = n, u = u, a = a,bax1 = bax1, bax2 = bax2) # for s in 1:nsim
  
  sim_stats <- do.call(rbind, sim_stats_list)
  
  # Compute avg operating characteristics from stats of interest 
  final_stats <- sim_stats %>% group_by(method) %>%
    summarize(bias = 100*mean(bias/a,na.rm=T),
              mse = mean(bias^2,na.rm=T),
              ATE = mean(ATE,na.rm=T),
              ci_cov = 100*mean(ci_cov,na.rm=T),
              power=100*mean(pow,na.rm=T))
  
  # make note of the grid point
  final_stats$u=u ; final_stats$a=a ; final_stats$n=n ;
  final_stats$bax1=bax1 ; final_stats$bax2=bax2
  
  return(final_stats)
} # calc_stats

get_results <- function(methods,
                        sig_u_grid, # me variance
                        ba_grid,    # main tmt effect
                        bax1_grid, bax2_grid, # interactions
                        n_grid,     # n
                        nsim=100,
                        save_intermediate=T) {
  
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
  #'  - save_intermediate: 0/1, if 1, will save a copy of op_chars every time
  #'                       a new row is updated
  #'  
  #' OUTPUTS:
  #' - A dataframe containing calculated operating characteristics for 
  #'   each method and parameter combination (e.g. a given row will contain info
  #'   on the avg percent bias, mse, 95% CI coverage, etc for a specific correction
  #'   method and specific values of sig_u, ba, n, etc)
  
  # Initialize dataframe for each stat of interest
  scen_df <- expand.grid(sig_u = sig_u_grid, ba = ba_grid, n = n_grid,
                         bax1 = bax1_grid, bax2 = bax2_grid, 
                         KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, ...) {
    
    u <- scen$sig_u # me variance
    a <- scen$ba # main effect
    bax1 <- scen$bax1 ; bax2 <- scen$bax2
    n <- scen$n # sample size
    
    # Get operating characteristics for current grid point
    op_tmp <- get_stats_table(methods,u,a,bax1,bax2,n,nsim)
    
    return(op_tmp)
    
  }, methods = methods, nsim = nsim)
  
  # Save the op_chars table after each grid point is done so we don't need to
  # wait for the entire run to be done (can plot grid points as they come in)  
  op_chars <- do.call(rbind, op_list)
  
  
  return(op_chars)
} # get_results 



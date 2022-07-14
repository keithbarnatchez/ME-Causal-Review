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

calc_stats <- function(data,methods,a,s) {
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
  CI_ideal <- ATE_ideal[[2]] ; CI_naive <- ATE_naive[[2]]
  
  # Bias
  bias_ideal <- ATE_ideal[[1]]-a ; bias_naive <- ATE_naive[[1]]-a
  
  # CI coverage
  ci_cov_ideal <- ifelse( (CI_ideal[1] <= a) & (CI_ideal[2] >= a) , 1, 0  )
  ci_cov_naive <- ifelse( (CI_naive[1] <= a) & (CI_naive[2] >= a) , 1, 0  )
  
  # Power
  pow_ideal <- ifelse(CI_ideal[[1]]>0 | CI_ideal[[2]] < 0,1,0) 
  pow_naive <- ifelse(CI_naive[[1]]>0 | CI_naive[[2]] < 0,1,0) 
  
  stats <- data.frame(bias=bias_ideal,ATE=ATE_ideal[[1]],ci_cov=ci_cov_ideal,pow=pow_ideal,method='Ideal',iteration=s)
  stats <- rbind(stats, data.frame(bias=bias_naive,ATE=ATE_naive[[1]],ci_cov=ci_cov_naive,pow=pow_naive,method='Naive',iteration=s))
                         
  # SIMEX
  if ('simex_ind' %in% methods) {
    # Get ATE estimate
    ATE <- simex_indirect(data) 
    
    # Get bias
    bias <- ATE[[1]]-a ; CI <- ATE[[2]]
    ci_cov <- ifelse( (CI[1] <= a) & (a <= CI[2]), 1 ,0  )
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE[[1]],
                                     ci_cov=ci_cov,pow=pow,
                                     method='SIMEX',iteration=s))
  }
  
  # PSC using IPTW
  if ('psc' %in% methods) {
    # Get ATE estimate
    ATE <- psc(data) 
    bias <- ATE[[1]]-a ; CI <- ATE[[2]]
    ci_cov <- ifelse( (CI[1] <= a) & (a <= CI[2]), 1 ,0  )
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE[[1]],ci_cov=ci_cov,
                                     pow=pow,method='PSC-IPTW',iteration=s))
  }
  
  if ('psc_reg' %in% methods) {
    ATE <- psc(data,iptw=0) 
    bias <- ATE[[1]]-a ; CI <- ATE[[2]]
    ci_cov <- ifelse( (CI[1] <= a) & (a <= CI[2]), 1 ,0  )
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE[[1]],ci_cov=ci_cov,
                                     pow=pow,method='PSC',iteration=s))  
  }
  
  # IV 
  if ('iv' %in% methods) {
    # Get ATE estimate
    ATE <- iv_confounder(data)
    bias <- ATE[[1]]-a ; CI <- ATE[[2]]
    ci_cov <- ifelse( (CI[1] <= a) & (a <= CI[2]), 1 ,0  )
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE[[1]],ci_cov=ci_cov,
                                     pow=pow,method='IV',iteration=s))
  }
  
  # MIME 
  if ('mime' %in% methods) {
    ATE <- mime(data)
    bias <- ATE[[1]]-a ; CI <- ATE[[2]]
    ci_cov <- ifelse( (CI[1] <= a) & (a <= CI[2]), 1 ,0  )
    pow <- ifelse(CI[[1]]>0 | CI[[2]] < 0,1,0) # power
    stats <- rbind(stats, data.frame(bias=bias,ATE=ATE[[1]],ci_cov=ci_cov,
                                     pow=pow,method='MIME',iteration=s))
  }
  
  return(stats)
}



get_stats_table <- function(methods,
                            u,a,n,
                            rho,psi,ax,b,
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
  sim_stats_list <- lapply(1:nsim, function(s, n, u, a, b, rho, psi, ax, ...) {
    
    # Keep track of progress
    print(paste('On iteration',s,'of',nsim))
    
    # simulate data for current iteration
    data <- generate_data(n,sig_u=u,ba=a,binary=b,
                          rho=rho, psi=psi, ax=ax)

    # Calculate stats of interest (e.g. bias, whether CI covers true param val, etc)
    return(calc_stats(data,methods,a,s))
    
  }, n = n, u = u, a = a, b = b, rho=rho, psi=psi, ax=ax) # for s in 1:nsim
  
  sim_stats <- do.call(rbind, sim_stats_list)

  # Compute avg operating characteristics from stats of interest 
  final_stats <- sim_stats %>% group_by(method) %>%
    summarize(bias = 100*mean(bias/a,na.rm=T),
              mse = mean(bias^2,na.rm=T),
              ATE = mean(ATE,na.rm=T),
              ci_cov = 100*mean(ci_cov,na.rm=T),
              power=100*mean(pow,na.rm=T))
  
  # make note of the grid point
  final_stats$u=u ; final_stats$a=a ; final_stats$n=n ; final_stats$b=b
  final_stats$rho = rho ; final_stats$psi = psi ; final_stats$ax=ax
  
  return(final_stats)
} # calc_stats

get_results <- function(methods,
                        sig_u_grid, # me variance
                        ba_grid,    # tmt effect
                        n_grid,     # n
                        rho_grid,   # corr bw x and z
                        psi_grid,   # corr bw x and instrument v
                        ax_grid,
                        bin_grid,   # whether outcome is binary or not
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
                         rho=rho_grid, psi = psi_grid, ax=ax_grid,
                         bin = bin_grid,
                         KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE)
  
  scen_list <- split(scen_df, seq(nrow(scen_df)))
  
  op_list <- lapply(scen_list, function(scen, methods, nsim, ...) {
    
    u <- scen$sig_u # me variance
    a <- scen$ba # effect size
    n <- scen$n # sample size
    rho <- scen$rho
    psi <- scen$psi
    ax <- scen$ax
    b <- scen$bin # binary outcome indicator
    
    
    # Get operating characteristics for current grid point
    op_tmp <- get_stats_table(methods,u,a,n,rho,psi,ax,b,nsim)
    
    return(op_tmp)
    
  }, methods = methods, nsim = nsim)
  
  # Save the op_chars table after ach grid point is done so we don't need to
  # wait for the entire run to be done (can plot grid points as they come in)  
  op_chars <- do.call(rbind, op_list)
  
    
  return(op_chars)
} # get_results 

bias_plot <- function(results) {
  
  results %>% ggplot(aes(x=bias))
  
}
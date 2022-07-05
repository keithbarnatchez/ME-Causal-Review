# sim_main.R
#
# Main file for performing simulations. User can adjust simulation parameters, 
# and provide vectors corresponding to parameters if they want to perform
# simulations for different combinations of parameters
#
# Example:
#
#
#
#
# Output:
# Running this code will result in two outputs:
# 1) a csv file containing the results of performing get_results() with the 
#    user-specified parameter values/combinations, named
#    "sim_resuls_<datestring>.csv", where datestring is generatef with 
# 2) a text file containing the specified parameter values
# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':','-',Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_',gsub(' ','_',flnm),'.csv',sep = '')

simdir <- '../output/sim_results/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name 
# ------------------------------------------------------------------------------
# Set up simulation parameters
methods <- c('iv','mime','psc') 
sig_u_grid <- c(0.1,0.2,0.3,0.5,0.9) ; bt_grid <- c(1) ; n_grid <- c(5000)
bin_grid <- c(0)
# ------------------------------------------------------------------------------
# Run simulations, store results
op_chars <- get_results(methods,sig_u_grid,bt_grid,n_grid, bin_grid, nsim=100)

# Output the results as a csv
write.csv(op_chars, file='../output/sim_results/')



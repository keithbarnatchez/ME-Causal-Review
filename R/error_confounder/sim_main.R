# sim_main.R
#
# Main file for performing simulations. User can adjust simulation parameters, 
# and provide vectors corresponding to parameters if they want to perform
# simulations for different combinations of parameters
#
# Output:
# Running this code will result in two outputs:
# 1) a csv file containing the results of performing get_results() with the 
#    user-specified parameter values/combinations, named
#    "sim_resuls_<datestring>.csv", where datestring is generated with 
# 2) a text file containing the specified parameter values

# ------------------------------------------------------------------------------
library(mvtnorm) # for simulating data
library(tidyverse)
library(mice) # for multiple imputation
library(AER) # for instrumental variables
library(SuperLearner)
library(mgcv)
library(parallel)

# Code folder is root directory in case you're working interactively
source('~/Github/ME-Causal-Review/R/error_confounder/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_confounder/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_confounder/output_functions.R')
source('~/Github/ME-Causal-Review/R/utility.R')

options(show.error.locations = TRUE)

# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':', '-', Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_', gsub(' ', '_', flnm), '.csv', sep = '')

simdir <- '../../output/sim_results/' # directory
fullpath <- paste(simdir, flnm, sep = '') # construct path to final file name 

# ------------------------------------------------------------------------------
# Set up simulation parameters
methods <- c('psc', 'simex', 'iv', 'mime') 
sig_u_grid <- c(0, 0.25, 0.5, 0.75, 1) # ME variances
ba_grid <- c(1) # treatment effect
aw_grid <- c(0.25) # coef of w in the ps model
bw_grid <- c(-1)
n_grid <- c(1000, 2000) # sample size
bin_grid <- c(FALSE, TRUE) # 0/1 continuous/binary

# ------------------------------------------------------------------------------
# Run simulations, store results
op_chars <- get_results(methods,
                        sig_u_grid,
                        ba_grid,
                        n_grid, 
                        rho_grid,
                        psi_grid,
                        ax_grid,
                        bin_grid, nsim=100)

# Output the results as a csv

write.csv(op_chars, file = fullpath)

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
# .libPaths("~/apps/R_4.1.0")
library(mvtnorm) # for simulating data
library(tidyverse)
library(mice) # for multiple imputation
library(AER) # for instrumental variables
library(SuperLearner)
library(mgcv)
library(AIPW)
library(parallel)

# Code folder is root directory in case you're working interactively
source('~/Github/ME-Causal-Review/R/error_confounder/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_confounder/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_confounder/output_functions.R')
source('~/Github/ME-Causal-Review/R/ATE.R')

options(show.error.locations = TRUE)
# ------------------------------------------------------------------------------
# Set up relevant paths for output
# flnm <- gsub(':', '-', Sys.time()) # filename suffixed with date/time
# flnm <- paste('sim_results_', gsub(' ', '_', flnm), '.csv', sep = '')

flnm <- "confounder_results.csv"
simdir <- '~/Github/ME-Causal-Review/output/sim_results/' # directory
fullpath <- paste(simdir, flnm, sep = '') # construct path to final file name 

# ------------------------------------------------------------------------------
# Set up simulation parameters
methods <- c('psc', 'rc', 'simex', 'iv', 'mime') 
sig_u_grid <- c(0.1, 0.25, 0.5, 0.75, 0.9) # ME variances
ba_grid <- c(1) # treatment effect
aw_grid <- c(0.5, -0.5) # coef of w in the ps model
bw_grid <- c(0.5, -0.5) # coef of w in the outcome model
n_grid <- c(1000, 2000) # sample size
bin_grid <- c(FALSE) # T/F continuous/binary


# ------------------------------------------------------------------------------
# Run simulations, store results
op_chars <- get_results(methods,
                        sig_u_grid,
                        aw_grid,
                        ba_grid,
                        bw_grid,
                        n_grid, 
                        bin_grid, 
                        mc.cores = 20,
                        nsim = 500)

# Output the results as a csv

write.csv(op_chars, file = fullpath)

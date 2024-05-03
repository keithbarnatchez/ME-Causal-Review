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
library(dplyr)
library(parallel)

# Code folder is root directory in case you're working interactively
source('~/Github/ME-Causal-Review/R/error_confounder/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_confounder/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_confounder/output_functions.R')
source('~/Github/ME-Causal-Review/R/SRF.R')

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
methods <- c('rc', 'mime', 'simex', 'iv', 'cv') 
sig_u_grid <- c(0.1, 0.25, 0.5, 0.75, 0.9) # ME variances
ba_grid <- c(1) # treatment effect
aw_grid <- c(0.75) # coef of w in the ps model
mis_grid <- c("base", "ps-mis", "out-mis") # coef of w in the ps model
n_grid <- c(1000) # sample size
rho_grid <- c(0.15, 0.85) # high and low numbers on (-1,1)


# ------------------------------------------------------------------------------
# Run simulations, store results
op_chars <- get_results(methods,
                        sig_u_grid = sig_u_grid,
                        ba_grid = ba_grid,
                        aw_grid = aw_grid,
                        mis_grid = mis_grid,
                        n_grid = n_grid, 
                        rho_grid = rho_grid, 
                        mc.cores = 24,
                        nsim = 500)

# Output the results as a csv

write.csv(op_chars, file = fullpath)

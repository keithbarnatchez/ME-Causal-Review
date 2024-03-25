rm(list=ls())

library(mvtnorm) # for simulating data
library(tidyverse)
library(mice) # for multiple imputation
library(AER) # for instrumental variables
library(SuperLearner)
library(mgcv)
library(geex)
library(parallel)

# Code folders
source('~/Github/ME-Causal-Review/R/error_exposure/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/output_functions.R')
source('~/Github/ME-Causal-Review/R/ERF.R')

#-------------------------------------------------
# Test correction functions 
#-------------------------------------------------

# Generate the data
data <- generate_data(n = 1000, sig_u = 0.5, binary = FALSE, ba = 1) 

# Get IPW from different approaches
ideal <- erf_ideal(data)
naive <- erf_naive(data)
csme <- erf_csme(data) # conditional score method
rc <- erf_rc(data) # propensity score calibration 
iv <- erf_iv(data) # IV
mime <- erf_mime(data) # multiple imputation
simex <- erf_simex(data) # simulation extrapolation
cv <- erf_cv(data) # simulation extrapolation

# ------------------------------------------------
# Test the get_results() function
# ------------------------------------------------

methods <- c('iv', 'cv') 
sig_u_grid <- c(0.1,0.3,0.5,0.9) 
ba_grid <- c(-1)
n_grid <- c(2000)
bin_grid <- c(FALSE)

op_chars <- get_results(methods,
                        sig_u_grid,
                        ba_grid,
                        n_grid, 
                        bin_grid, 
                        mc.cores = 1,
                        nsim = 100)

# ------------------------------
# Plot 
# ------------------------------

op_chars %>% ggplot(aes(x = sig_u, y = bias, color = method)) + geom_line() + theme_bw()

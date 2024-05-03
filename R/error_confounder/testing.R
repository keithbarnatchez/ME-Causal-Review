rm(list=ls())
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
source('~/Github/ME-Causal-Review/R/SRF.R')

#-------------------------------------------------
# Test correction functions 
#-------------------------------------------------

# Generate the data
data <- generate_data(n = 1000, sig_u = 1, binary = FALSE, ba = 1) 

# Get IPW from different approaches
ideal <- srf_naive(data)
naive <- srf_ideal(data)
rc <- srf_rc(data) # regression calibration 
iv <- srf_iv(data) # IV
mime <- srf_mime(data) # multiple imputation
simex <- srf_simex(data) # sim ext.
cv <- srf_cv(data)

# ------------------------------------------------
# Test the get_results() function
# ------------------------------------------------

methods <- c('rc','mime','iv','cv') 
sig_u_grid <- c(0.1,0.3,0.5,0.9) 
aw_grid = c(0.5)
mis_grid = c("out-mis")
ba_grid <- c(1)
n_grid <- c(1000)
rho_grid <- c(0.5)

op_chars <- get_results(methods,
                        sig_u_grid = sig_u_grid,
                        ba_grid = ba_grid,
                        aw_grid = aw_grid,
                        mis_grid = mis_grid,
                        n_grid = n_grid, 
                        rho_grid = rho_grid, 
                        mc.cores = 25,
                        nsim = 100)

# ------------------------------
# Plot 
# ------------------------------

op_chars %>% ggplot(aes(x = sig_u, y = bias, color = method)) + geom_line() + theme_bw()

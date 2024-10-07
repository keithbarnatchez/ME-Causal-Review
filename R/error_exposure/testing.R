rm(list=ls())

library(mvtnorm) # for simulating data
library(tidyverse)
library(mice) # for multiple imputation
library(AER) # for instrumental variables
library(SuperLearner)
library(mgcv)
library(KernSmooth)
library(parallel)

# Code folders
source('~/Github/ME-Causal-Review/R/error_exposure/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/output_functions.R')
source('~/Github/ME-Causal-Review/R/SRF.R')

#-------------------------------------------------
# Test correction functions 
#-------------------------------------------------

# Generate the data
data <- generate_data(n = 1000, sig_u = 0.5, binary = FALSE, ba = 1)

# Get IPW from different approaches
ideal <- srf_ideal(data)
naive <- srf_naive(data)
rc <- srf_rc(data) # propensity score calibration 
iv <- srf_iv(data) # IV
mime <- srf_mime(data) # multiple imputation
simex <- srf_simex(data) # simulation extrapolation
cv <- srf_cv(data) # simulation extrapolation

# ------------------------------------------------
# Test the get_results() function
# ------------------------------------------------

methods <- c('iv', 'cv') 
sig_u_grid <- c(0.1,0.3,0.5,0.9) 
ba_grid <- c(-1)
aw_grid <- 0.5
bw_grid <- -0.5
mis_grid <- c("ps-mis")
n_grid <- c(2000)

op_chars <- get_results(methods,
                        sig_u_grid = sig_u_grid,
                        ba_grid = ba_grid, 
                        aw_grid = aw_grid, 
                        bw_grid = bw_grid,
                        mis_grid = mis_grid,
                        n_grid = n_grid, 
                        mc.cores = 24,
                        nsim = 100)

# ------------------------------
# Plot 
# ------------------------------

op_chars %>% ggplot(aes(x = sig_u, y = bias, color = method)) + geom_line() + theme_bw()

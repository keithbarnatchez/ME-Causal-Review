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
source('~/Github/ME-Causal-Review/R/ATE.R')

#-------------------------------------------------
# Test correction functions 
#-------------------------------------------------

# Generate the data
data <- generate_data(n = 1000, sig_u = 1, binary = FALSE, ba = 1) 

# Get IPW from different approaches
ideal <- ate_naive(data)
naive <- ate_ideal(data)
rc <- ate_rc(data, method = "aipw") # regression calibration 
psc <- ate_psc(data) # propensity score calibration
iv <- ate_iv(data) # IV
mime <- ate_mime(data, method = "aipw")
simex <- ate_simex(data, method = "aipw")

# ------------------------------------------------
# Test the get_results() function
# ------------------------------------------------

methods <- c('psc','iv') 
sig_u_grid <- c(0.1,0.3,0.5,0.9) 
aw_grid = c(0.5)
ba_grid <- c(1)
bw_grid <- c(-1)
n_grid <- c(2000)
bin_grid <- c(FALSE)

op_chars <- get_results(methods,
                        sig_u_grid,
                        aw_grid,
                        ba_grid,
                        bw_grid,
                        n_grid, 
                        bin_grid, 
                        mc.cores = 1,
                        nsim = 100)

# ------------------------------
# Plot 
# ------------------------------

op_chars %>% ggplot(aes(x = sig_u, y = bias, color = method)) + geom_line() + theme_bw()

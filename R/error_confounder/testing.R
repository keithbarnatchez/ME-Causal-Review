rm(list=ls())
library(mvtnorm) # for simulating data
library(simex) # for the indirect SIMEX approach
library(tidyverse)
library(mice) # for multiple imputation
library(ipw)
library(AER) # for instrumental variables
library(betareg) # for using beta regresison in PSC 

# Code folder is root directory in case you're working interactively
source('correction_functions.R')
source('data_functions.R')
source('output_functions.R')

# Simulation code 
n <- 10000 # obs 

# Generate the data
data <- generate_data(n,sig_u = 0.3,binary=0,
                      ba=1) 
#-------------------------------------------------
# Test correction functions 

# Get IPW from different approaches
ate_psc <- psc(data) # propensity score calibration 
ate_iv <- iv_confounder(data) # IV
ate_mime <- mime(data)
ate_simex <- simex_indirect(data)

# ------------------------------------------------
# Test the get_results() function

methods <- c('psc') 
sig_u_grid <- c(0.1,0.2,0.3,0.5,0.9) ; bt_grid <- c(1) ; n_grid <- c(5000)
bin_grid <- c(0)
op_chars <- get_results(methods,sig_u_grid,bt_grid,n_grid, bin_grid, nsim=100)

# ------------------------------
# Plot 

op_chars2 %>% filter(a==1) %>% ggplot(aes(x=u, y=bias, color=method)) + geom_line() +
  theme_bw()

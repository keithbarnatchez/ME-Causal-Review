rm(list=ls())
library(mvtnorm)
library(simex)
library(tidyverse)
library(mice)
library(ipw)
library(AER)

# Code folder is root directory in case you're working interactively
source('correction_functions.R')
source('data_functions.R')
source('output_functions.R')

# Simulation code 
n <- 10000 # obs

# Generate the data
data <- generate_data(n,sig_u = 5,binary=1,
                      ba=0.3)

#-------------------------------------------------
# Test correction functions 

# Get IPW from different approaches
ate_psc <- psc(data) # propensity score calibration
ate_iv <- iv_confounder(data) # IV
ate_mime <- mime(data)

# ------------------------------------------------
# Test the get_results() function

methods <- c('iv')
sig_u_grid <- c(0.1,0.2,0.3) ; bt_grid <- c(1,2) ; n_grid <- c(100,1000,5000)
bin_grid <- c(0,1)
op_chars <- get_results(methods,sig_u_grid,bt_grid,n_grid, bin_grid)


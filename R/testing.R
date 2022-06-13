rm(list=ls())
library(mvtnorm)
library(simex)

# Code folder is root directory in case you're working interactively
source('correction_functions.R')
source('data_functions.R')

# Simulation code 
n <- 10000 # obs
v_share <- 0.1 # share in validation data
v_idx <- rbinom(n,size=1,prob=v_share) # logical noting which obs are in validation data

# Generate the data
data <- generate_data(n,sig_u = 0.2)

# Get IPW from different approaches
iptw_psc <- psc(data,v_idx)
iptw_simex <- simex_indirect(data,v_idx)


# exposure_sim_main.R
#
# Code for implementing exposure ME simulation exercise
rm(list = ls())

# ------------------------------------------------------------------------------
## Preliminaries
library(SuperLearner)
library(parallel)
library(abind)
library(geex)
library(AER)
library(tidyverse)
library(KernSmooth)
library(dplyr)

# Code folders
source('~/Github/ME-Causal-Review/R/error_exposure/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/output_functions.R')
source('~/Github/ME-Causal-Review/R/ERF2.R')

# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- "exposure_results.csv"
simdir <- '~/Github/ME-Causal-Review/output/sim_results/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name 

# ------------------------------------------------------------------------------
# Run It!

methods <- c('rc', 'simex', 'iv', 'mime', 'cv')
sig_u_grid <- c(0.1,0.25,0.5,0.75,0.9) 
ba_grid <- c(0.5, 1)
aw_grid <- c(0.25, 0.75) # coef of w in the ps model
n_grid <- c(1000, 2000) # sample size

op_chars <- get_results(methods,
                        sig_u_grid,
                        ba_grid,
                        aw_grid,
                        n_grid, 
                        mc.cores = 24,
                        nsim = 100)

write.csv(op_chars, file=fullpath)
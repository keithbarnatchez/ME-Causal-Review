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
library(tidyverse)

# Code folders
source('~/Github/ME-Causal-Review/R/error_exposure/correction_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/data_functions.R')
source('~/Github/ME-Causal-Review/R/error_exposure/output_functions.R')
source('~/Github/ME-Causal-Review/R/ERF.R')

# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':','-',Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_',gsub(' ','_',flnm),'.csv',sep = '')

simdir <- '~/Github/ME-Causal-Review/output/sim_results/exposure/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name 

# ------------------------------------------------------------------------------
# Run It!

methods <- c('rc', 'simex', 'iv', 'mime')
sig_u_grid <- c(0.1,0.3,0.5,0.9) 
ba_grid <- c(1)
n_grid <- c(2000)
bin_grid <- c(FALSE)

op_chars <- get_results(methods,
                        sig_u_grid,
                        ba_grid,
                        n_grid, 
                        bin_grid, 
                        mc.cores = 1,
                        nsim = 100)


write.csv(thecheck, file=fullpath)
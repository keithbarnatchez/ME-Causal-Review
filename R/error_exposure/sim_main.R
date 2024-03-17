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

# ------------------------------------------------------------------------------
# Set up relevant paths for output
flnm <- gsub(':','-',Sys.time()) # filename suffixed with date/time
flnm <- paste('sim_results_',gsub(' ','_',flnm),'.csv',sep = '')

simdir <- '../../output/sim_results/exposure/' # directory
fullpath <- paste(simdir,flnm,sep='') # construct path to final file name 

# ------------------------------------------------------------------------------
# Code for generating and fitting data
# Assuming directory is set to this code file's location
source("/../ERF.R")
source("correction_functions.R")
source("output_functions.R")

# ------------------------------------------------------------------------------
# Run It!

methods <- c('simex')
sig_u_grid <- c(0.1,0.5,0.9)
ba_grid <- c(1)
bax1_grid <- c(-0.25) ; bax2_grid <- c(0.5)
n_grid <- c(1500) ; nsim <- 100
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
thecheck <- get_results(methods,
                        sig_u_grid,
                        ba_grid,bax1_grid,bax2_grid,
                        n_grid,
                        nsim=50)


write.csv(thecheck, file=fullpath)
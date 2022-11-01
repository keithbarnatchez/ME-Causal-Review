rm(list=ls())
library(mvtnorm) # for simulating data
library(simex) # for the indirect SIMEX approach
library(tidyverse)
library(mice) # for multiple imputation
library(ipw)
library(AER) # for instrumental variables
library(betareg) # for using beta regression in PSC 
library(SuperLearner)

# Code folder is root directory in case you're working interactively
source('exposure_correction_functions.R')
source('exposure_data_functions.R')

# 
n <- 1e5
data <- generate_data(n)
iv_exposure(data)


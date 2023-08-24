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
source('correction_functions.R')
source('data_functions.R')
source('output_functions.R')

# Simulation code 
n <- 1.5e3 # obs 

# Generate the data
data <- generate_data(n,sig_u = 0.3,binary=0,
                      ba=1) 
#-------------------------------------------------
# Test AIPW implementation

# Fit outcome model
rhs_vars <- data %>% select(X,A,Z)
slhat <- SuperLearner(Y=data$Y,  X=rhs_vars, family=gaussian, SL.library = 'SL.glm')$SL.predict

# Fit PS model
glmhat <- lm(Y ~ X + A + Z,data=data)$fitted.values

# Get muhat|X,A=0, muhat|X,A=1
Y<-data$Y
A<-data$A
W<-data$W
Z<-data$Z
tau2<-0.3
#-------------------------------------------------
# Test correction functions 

# Get IPW from different approaches
ate_psc_reg <- psc(data,nboot=100,iptw=0) # propensity score calibration 
ate_psc <- psc(data,nboot=100,iptw=1)
ate_iv <- iv_confounder(data) # IV
ate_mime <- mime(data)
ate_simex <- simex_indirect(data,nboot=0)
ate_simex_direct <- simex_direct(data$Y,data$A,data$W,data$Z,family=gaussian(),tau2=0.3)
# ------------------------------------------------
# Test the get_results() function

methods <- c('psc_reg') 
sig_u_grid <- c(0.1,.3,0.5,0.9) ; bt_grid <- c(1) ; n_grid <- c(5000)
bin_grid <- c(0)
rho_grid = 0.5; psi_grid = 0.3 ; ax_grid=1/4
op_chars <- get_results(methods,sig_u_grid,bt_grid,n_grid,
                        rho_grid, psi_grid, ax_grid, bin_grid, nsim=10)

# ------------------------------
# Plot 

op_chars %>% ggplot(aes(x=u, y=bias, color=method)) + geom_line() +
  theme_bw()


sig_u_grid <- c(0.1,0.5,0.9) # ME variances
ba_grid <- c(1) # treatment effect
n_grid <- c(5000) # sample size
bin_grid <- c(0) # 0/1 continuous/binary
rho_grid <- c(0.2,0.8) # corr b/w x and z
psi_grid <- c(0.5) # corr b/w v and x
ax_grid <- 0.25 # coeff of x in the ps model

n <- 5000 ; sig_u <- 0.1 ; binary=0 ; rho=0.2;n=5000;sig_e=1


data <- generate_data(n,sig_u = 0.3,binary=0,
                      ba=1) 

temp_data <- generate_data(n,sig_u = 0.9,binary=0,
                           rho = 0.2,
                           ba=1) 
ate_simex <- simex_indirect(temp_data,nboot=0)
naive <- ate_naive(temp_data)
ate_iv <- iv_confounder(temp_data)
ate_mime <- mime(data)

# testing AIPW
aipw_test <- aipw(data)

predicted <- model.matrix(simex_model$model)%*%simex_model$coefficients
pred2 <- exp(predicted)/(1+exp(predicted))



